struct AddCostSpec
    variable_type::Type
    component_type::Type
    has_status_variable::Bool
    has_status_parameter::Bool
    sos_status::SOS_STATUS_VARIABLE
    multiplier::Float64
    variable_cost::Union{Nothing, Function}
    start_up_cost::Union{Nothing, Function}
    shut_down_cost::Union{Nothing, Function}
    fixed_cost::Union{Nothing, Function}
    addtional_linear_terms::Dict{String, Symbol}
end

function AddCostSpec(;
    variable_type,
    component_type,
    has_status_variable = false,
    has_status_parameter = false,
    sos_status = NO_VARIABLE,
    multiplier = OBJECTIVE_FUNCTION_POSITIVE,
    variable_cost = nothing,
    start_up_cost = nothing,
    shut_down_cost = nothing,
    fixed_cost = nothing,
    addtional_linear_terms = Dict{String, Symbol}(),
)
    return AddCostSpec(
        variable_type,
        component_type,
        has_status_variable,
        has_status_parameter,
        sos_status,
        multiplier,
        variable_cost,
        start_up_cost,
        shut_down_cost,
        fixed_cost,
        addtional_linear_terms,
    )
end

function AddCostSpec(
    ::Type{<:T},
    ::Type{<:U},
    ::PSIContainer,
) where {T <: PSY.Component, U <: AbstractDeviceFormulation}
    error("AddCostSpec is not implemented for $T / $U")
end

set_addtional_linear_terms!(spec::AddCostSpec, key, value) =
    spec.addtional_linear_terms[key] = value

function add_service_variables!(spec::AddCostSpec, services)
    for service in services
        name = PSY.get_name(service)
        set_addtional_linear_terms!(spec, name, make_variable_name(name, typeof(service)))
    end
    return
end

"""
Add variables to the PSIContainer for a service.
"""
function cost_function!(
    psi_container::PSIContainer,
    devices::IS.FlattenIteratorWrapper{T},
    ::DeviceModel{T, U},
    ::Type{<:PM.AbstractPowerModel},
    feedforward::Union{Nothing, AbstractAffectFeedForward} = nothing,
) where {T <: PSY.Component, U <: AbstractDeviceFormulation}
    for d in devices
        spec = AddCostSpec(T, U, psi_container)
        @debug T, spec
        services = PSY.get_services(d)
        add_service_variables!(spec, services)
        add_to_cost!(psi_container, spec, PSY.get_operation_cost(d), d)
    end
    return
end

function add_to_cost_expression!(
    psi_container::PSIContainer,
    cost_expression::JuMP.AbstractJuMPScalar,
)
    T_ce = typeof(cost_expression)
    T_cf = typeof(psi_container.cost_function)
    if T_cf <: JuMP.GenericAffExpr && T_ce <: JuMP.GenericQuadExpr
        psi_container.cost_function += cost_expression
    else
        JuMP.add_to_expression!(psi_container.cost_function, cost_expression)
    end
    return
end

function has_on_variable(
    psi_container::PSIContainer,
    ::Type{T};
    variable_type = OnVariable,
) where {T <: PSY.Component}
    # get_variable can't be used because the default behavior is to error if variables is not present
    return !isnothing(get(
        psi_container.variables,
        make_variable_name(variable_type, T),
        nothing,
    ))
end

function has_on_parameter(psi_container::PSIContainer, ::Type{T}) where {T <: PSY.Component}
    if !model_has_parameters(psi_container)
        return false
    end
    return !isnothing(get(
        psi_container.parameters,
        encode_symbol(OnVariable, string(T)),
        nothing,
    ))
end

function _get_pwl_vars_container(psi_container::PSIContainer)
    if !haskey(psi_container.variables, :PWL_cost_vars)
        time_steps = model_time_steps(psi_container)
        contents = Dict{Tuple{String, Int, Int}, Any}()
        container = JuMP.Containers.SparseAxisArray(contents)
        assign_variable!(psi_container, :PWL_cost_vars, container)
    else
        container = get_variable(psi_container, :PWL_cost_vars)
    end
    return container
end

function slope_convexity_check(slopes::Vector{Float64})
    flag = true
    for ix in 1:(length(slopes) - 1)
        if slopes[ix] > slopes[ix + 1]
            @debug slopes
            return flag = false
        end
    end
    return flag
end

@doc raw"""
Returns True/False depending on compatibility of the cost data with the linear implementation method

Returns ```flag```

# Arguments

* cost_ : container for quadratic and linear factors
"""
function pwlparamcheck(cost_)
    slopes = PSY.get_slopes(cost_)
    # First element of the return is the average cost at P_min.
    # Shouldn't be passed for convexity check
    return slope_convexity_check(slopes[2:end])
end

function linear_gen_cost!(
    psi_container::PSIContainer,
    var_name::Symbol,
    component_name::String,
    linear_term::Float64,
    time_period::Int,
)
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    variable = get_variable(psi_container, var_name)[component_name, time_period]
    gen_cost = sum(variable) * linear_term
    add_to_cost_expression!(psi_container, gen_cost * dt)
    return
end

@doc raw"""
Returns piecewise cost expression using SOS Type-2 implementation for psi_container model.

# Equations

``` variable = sum(sos_var[i]*cost_data[2][i])```

``` gen_cost = sum(sos_var[i]*cost_data[1][i]) ```

# LaTeX

`` variable = (sum_{i\in I} c_{2, i} sos_i) ``

`` gen_cost = (sum_{i\in I} c_{1, i} sos_i) ``

Returns ```gen_cost```

# Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* variable::JuMP.Containers.DenseAxisArray{JV} : variable array
* cost_data::PSY.VariableCost{NTuple{2, Float64}} : container for quadratic and linear factors
"""
function pwl_gencost_sos!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component_name::String,
    cost_data::Vector{NTuple{2, Float64}},
    time_period::Int,
)
    var_name = make_variable_name(spec.variable_type, spec.component_type)
    variable = get_variable(psi_container, var_name)[component_name, time_period]
    settings_ext = get_ext(get_settings(psi_container))
    export_pwl_vars = get_export_pwl_vars(psi_container.settings)
    @debug export_pwl_vars
    total_gen_cost = JuMP.AffExpr(0.0)

    if spec.sos_status == NO_VARIABLE
        bin = 1.0
        @debug("Using Piecewise Linear cost function but no variable/parameter ref for ON status is passed. Default status will be set to online (1.0)")
    elseif spec.sos_status == PARAMETER
        param_key = encode_symbol(OnVariable, string(spec.component_type))
        bin =
            get_parameter_container(psi_container, param_key).parameter_array[component_name]
        @debug("Using Piecewise Linear cost function with parameter $(param_key)")
    elseif spec.sos_status == VARIABLE
        var_key = make_variable_name(OnVariable, spec.component_type)
        bin = get_variable(psi_container, var_key)[component_name, time_period]
        @debug("Using Piecewise Linear cost function with variable $(var_key)")
    else
        @assert false
    end

    gen_cost = JuMP.AffExpr(0.0)
    pwlvars = Array{JuMP.VariableRef}(undef, length(cost_data))
    for i in 1:length(cost_data)
        pwlvars[i] = JuMP.@variable(
            psi_container.JuMPmodel,
            base_name = "{$(variable)}_{sos}",
            start = 0.0,
            lower_bound = 0.0,
            upper_bound = 1.0
        )
        if export_pwl_vars
            container = _get_pwl_vars_container(psi_container)
            container[(component_name, time_period, i)] = pwlvars[i]
        end
        JuMP.add_to_expression!(gen_cost, cost_data[i][1] * pwlvars[i])
    end
    JuMP.@constraint(psi_container.JuMPmodel, sum(pwlvars) == bin)
    JuMP.@constraint(
        psi_container.JuMPmodel,
        pwlvars in MOI.SOS2(collect(1:length(pwlvars)))
    )
    JuMP.@constraint(
        psi_container.JuMPmodel,
        variable == sum([var_ * cost_data[ix][2] for (ix, var_) in enumerate(pwlvars)])
    )
    JuMP.add_to_expression!(total_gen_cost, gen_cost)

    return total_gen_cost
end

@doc raw"""
Returns piecewise cost expression using linear implementation for psi_container model.

# Equations

``` 0 <= pwl_var[i] <= (cost_data[2][i] - cost_data[2][i-1])```

``` variable = sum(pwl_var[i])```

``` gen_cost = sum(pwl_var[i]*cost_data[1][i]/cost_data[2][i]) ```

# LaTeX
`` 0 <= pwl_i <= (c_{2, i} - c_{2, i-1})``

`` variable = (sum_{i\in I} pwl_i) ``

`` gen_cost = (sum_{i\in I}  pwl_i) c_{1, i}/c_{2, i} ``

Returns ```gen_cost```

# Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* variable::JuMP.Containers.DenseAxisArray{JV} : variable array
* cost_data::Vector{NTuple{2, Float64}} : container for quadratic and linear factors
"""
function pwl_gencost_linear!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component_name::String,
    cost_data::Vector{NTuple{2, Float64}},
    time_period::Int,
)
    var_name = make_variable_name(spec.variable_type, spec.component_type)
    variable = get_variable(psi_container, var_name)[component_name, time_period]
    settings_ext = get_ext(get_settings(psi_container))
    export_pwl_vars = get_export_pwl_vars(psi_container.settings)
    @debug export_pwl_vars
    total_gen_cost = JuMP.AffExpr(0.0)

    gen_cost = JuMP.AffExpr(0.0)
    total_gen = JuMP.AffExpr(0.0)
    for i in 1:length(cost_data)
        pwlvar = JuMP.@variable(
            psi_container.JuMPmodel,
            base_name = "{$(variable)}_{pwl_$(i)}",
            start = 0.0,
            lower_bound = 0.0,
            upper_bound = PSY.get_breakpoint_upperbounds(cost_data)[i]
        )
        if export_pwl_vars
            container = _get_pwl_vars_container(psi_container)
            container[(component_name, time_period, i)] = pwlvar
        end
        JuMP.add_to_expression!(gen_cost, PSY.get_slopes(cost_data)[i] * pwlvar)
        JuMP.add_to_expression!(total_gen, pwlvar)
    end
    JuMP.@constraint(psi_container.JuMPmodel, variable == total_gen)
    JuMP.add_to_expression!(total_gen_cost, gen_cost)
    return total_gen_cost
end

"""
Adds to the models costs represented by PowerSystems TwoPart costs
"""
function add_to_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    cost_data::PSY.TwoPartCost,
    component::PSY.Component,
)
    component_name = PSY.get_name(component)
    time_steps = model_time_steps(psi_container)
    @debug "TwoPartCost" component_name
    if !(spec.variable_cost === nothing)
        variable_cost = spec.variable_cost(cost_data)
        for t in time_steps
            variable_cost!(psi_container, spec, component_name, variable_cost, t)
        end
    else
        @warn "No variable cost defined for $component_name"
    end

    if !isnothing(spec.fixed_cost) && spec.has_status_variable
        @debug "Fixed cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(OnVariable, spec.component_type),
                component_name,
                spec.fixed_cost,
                t,
            )
        end
    end
    return
end

"""
Adds to the models costs represented by PowerSystems ThreePart costs
"""
function add_to_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    cost_data::PSY.ThreePartCost,
    component::PSY.Component,
)
    component_name = PSY.get_name(component)
    @debug "ThreePartCost" component_name
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    variable_cost = PSY.get_variable(cost_data)
    time_steps = model_time_steps(psi_container)
    for t in time_steps
        variable_cost!(psi_container, spec, component_name, variable_cost, t)
    end

    if !isnothing(spec.start_up_cost)
        @debug "Start up cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(StartVariable, spec.component_type),
                component_name,
                spec.start_up_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    if !isnothing(spec.shut_down_cost)
        @debug "Shut down cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(StopVariable, spec.component_type),
                component_name,
                spec.shut_down_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    if !isnothing(spec.fixed_cost) && spec.has_status_variable
        @debug "Fixed cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(OnVariable, spec.component_type),
                component_name,
                spec.fixed_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    return
end

"""
Adds to the models costs represented by PowerSystems Multi-Start costs.
"""
function add_to_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    cost_data::PSY.MultiStartCost,
    component::PSY.Component,
)
    component_name = PSY.get_name(component)
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    time_steps = model_time_steps(psi_container)

    if !isnothing(spec.fixed_cost) && spec.has_status_variable
        @debug "Fixed cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(OnVariable, spec.component_type),
                component_name,
                spec.fixed_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    if !isnothing(spec.shut_down_cost)
        @debug "Shut down cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(StopVariable, spec.component_type),
                component_name,
                spec.shut_down_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    # Original implementation had SOS by default
    variable_cost_data = PSY.get_cost(PSY.get_variable(cost_data))
    if !all(iszero.(last.(variable_cost_data)))
        for t in time_steps
            gen_cost =
                pwl_gencost_sos!(psi_container, spec, component_name, variable_cost_data, t)
            add_to_cost_expression!(psi_container, spec.multiplier * gen_cost * dt)
        end
    else
        @debug "No Variable Cost associated with $(component_name)"
    end

    # Start-up costs
    start_cost_data = PSY.get_start_up(cost_data)
    for (st, var_type) in enumerate(start_types)
        var_name = make_variable_name(var_type, spec.component_type)
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                var_name,
                component_name,
                start_cost_data[st] * spec.multiplier,
                t,
            )
        end
    end

    return
end

"""
Adds to the models costs represented by PowerSystems Market-Bid costs. Implementation for
devices PSY.ThermalMultiStart
"""
function add_to_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    cost_data::PSY.MarketBidCost,
    component::PSY.ThermalMultiStart,
)
    component_name = PSY.get_name(component)
    @debug "Market Bid" component_name
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    time_steps = model_time_steps(psi_container)
    initial_time = model_initial_time(psi_container)
    variable_cost_forecast = PSY.get_variable_cost(
        component,
        PSY.get_operation_cost(component);
        start_time = initial_time,
        len = length(time_steps),
    )
    variable_cost_forecast_values = TimeSeries.values(variable_cost_forecast)
    for t in time_steps
        variable_cost!(
            psi_container,
            spec,
            component_name,
            variable_cost_forecast_values[t],
            t,
        )
    end

    if !isnothing(spec.start_up_cost)
        start_cost_data = spec.start_up_cost(cost_data)
        for (st, var_type) in enumerate(start_types)
            var_name = make_variable_name(var_type, spec.component_type)
            for t in time_steps
                linear_gen_cost!(
                    psi_container,
                    var_name,
                    component_name,
                    start_cost_data[st] * spec.multiplier,
                    t,
                )
            end
        end
    end

    if spec.has_status_variable
        @debug "no_load cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(OnVariable, spec.component_type),
                component_name,
                PSY.get_no_load(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    if !isnothing(spec.shut_down_cost)
        @debug "Shut down cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(StopVariable, spec.component_type),
                component_name,
                spec.shut_down_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    #Service Cost Bid
    ancillary_services = PSY.get_ancillary_services(cost_data)
    for service in ancillary_services
        add_service_bid_cost!(psi_container, spec, component, service)
    end
    return
end

"""
Adds to the models costs represented by PowerSystems Market-Bid costs. Default implementation for any PSY.Component. Uses by default the cost in of the cold stages for
start up costs.
"""
function add_to_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    cost_data::PSY.MarketBidCost,
    component::PSY.Component,
)
    component_name = PSY.get_name(component)
    @debug "Market Bid" component_name
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    time_steps = model_time_steps(psi_container)
    initial_time = model_initial_time(psi_container)
    variable_cost_forecast = PSY.get_variable_cost(
        component,
        PSY.get_operation_cost(component);
        start_time = initial_time,
        len = length(time_steps),
    )
    variable_cost_forecast_values = TimeSeries.values(variable_cost_forecast)
    for t in time_steps
        variable_cost!(
            psi_container,
            spec,
            component_name,
            variable_cost_forecast_values[t],
            t,
        )
    end

    if !isnothing(spec.start_up_cost)
        start_cost_data = spec.start_up_cost(cost_data)
        var_name = make_variable_name(StartVariable, spec.component_type)
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                var_name,
                component_name,
                start_cost_data.cold * spec.multiplier,
                t,
            )
        end
    end

    if spec.has_status_variable
        @debug "no_load cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(OnVariable, spec.component_type),
                component_name,
                PSY.get_no_load(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    if !isnothing(spec.shut_down_cost)
        @debug "Shut down cost" component_name
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                make_variable_name(StopVariable, spec.component_type),
                component_name,
                spec.shut_down_cost(cost_data) * spec.multiplier,
                t,
            )
        end
    end

    # Service Cost Bid
    ancillary_services = PSY.get_ancillary_services(cost_data)
    for service in ancillary_services
        add_service_bid_cost!(psi_container, spec, component, service)
    end
    return
end

function add_service_bid_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component::PSY.Component,
    service::PSY.Service,
)
    return
end

function add_service_bid_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component::PSY.Component,
    service::PSY.Reserve{T},
) where {T <: PSY.ReserveDirection}
    time_steps = model_time_steps(psi_container)
    initial_time = model_initial_time(psi_container)
    forecast_data = PSY.get_services_bid(
        component,
        PSY.get_operation_cost(component),
        service;
        start_time = initial_time,
        len = length(time_steps),
    )
    forecast_data_values = TimeSeries.values(forecast_data)
    if eltype(forecast_data_values) == PSY.VariableCost{Float64}
        for t in time_steps
            linear_gen_cost!(
                psi_container,
                spec.addtional_linear_terms[PSY.get_name(service)],
                PSY.get_name(component),
                forecast_data_values,
                t,
            )
        end
    else
        error("Current version only supports linear cost bid for services, please change the forecast data for $(PSY.get_name(service))")
    end
    return
end

function add_service_bid_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component::PSY.Component,
    service::PSY.ReserveDemandCurve{T},
) where {T <: PSY.ReserveDirection}
    error("Current version doesn't supports cost bid for ReserveDemandCurve services, please change the forecast data for $(PSY.get_name(service))")
    return
end

@doc raw"""
Adds to the cost function cost terms for sum of variables with common factor to be used for cost expression for psi_container model.

    # Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* var_name::Symbol: The variable name
* component_name::String: The component_name of the variable container
* cost_component::PSY.VariableCost{Float64} : container for cost to be associated with variable
"""
function variable_cost!(
    ::PSIContainer,
    ::AddCostSpec,
    component_name::String,
    ::Nothing,
    ::Int,
)
    @debug "Empty Variable Cost" component_name
    return
end

@doc raw"""
Adds to the cost function cost terms for sum of variables with common factor to be used for cost expression for psi_container model.

    # Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* var_name::Symbol: The variable name
* component_name::String: The component_name of the variable container
* cost_component::PSY.VariableCost{Float64} : container for cost to be associated with variable
"""
function variable_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component_name::String,
    cost_component::PSY.VariableCost{Float64},
    time_period::Int,
)
    @debug "Linear Variable Cost" component_name
    var_name = make_variable_name(spec.variable_type, spec.component_type)
    cost_data = PSY.get_cost(cost_component)
    linear_gen_cost!(
        psi_container,
        var_name,
        component_name,
        cost_data * spec.multiplier,
        time_period,
    )
    return
end

@doc raw"""
Adds to the cost function cost terms for sum of variables with common factor to be used for cost expression for psi_container model.

# Equation

``` gen_cost = dt*sign*(sum(variable.^2)*cost_data[1] + sum(variable)*cost_data[2]) ```

# LaTeX

`` cost = dt\times sign (sum_{i\in I} c_1 v_i^2 + sum_{i\in I} c_2 v_i ) ``

for quadratic factor large enough. If the first term of the quadratic objective is 0.0, adds a
linear cost term `sum(variable)*cost_data[2]`

# Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* var_name::Symbol: The variable name
* component_name::String: The component_name of the variable container
* cost_component::PSY.VariableCost{NTuple{2, Float64}} : container for quadratic and linear factors
"""
function variable_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component_name::String,
    cost_component::PSY.VariableCost{NTuple{2, Float64}},
    time_period::Int,
)
    var_name = make_variable_name(spec.variable_type, spec.component_type)
    cost_data = PSY.get_cost(cost_component)
    if cost_data[1] >= eps()
        @debug "Quadratic Variable Cost" component_name
        resolution = model_resolution(psi_container)
        dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
        variable = get_variable(psi_container, var_name)[component_name, time_period]
        gen_cost = sum(variable .^ 2) * cost_data[1] + sum(variable) * cost_data[2]
        add_to_cost_expression!(psi_container, spec.multiplier * gen_cost * dt)
    else
        @debug "Quadratic Variable Cost with only linear term" component_name
        linear_gen_cost!(
            psi_container,
            var_name,
            component_name,
            cost_data[2] * spec.multiplier,
            time_period,
        )
    end
    return
end

@doc raw"""
Creates piecewise linear cost function using a sum of variables and expression with sign and time step included.

# Expression

```JuMP.add_to_expression!(gen_cost, c)```

Returns sign*gen_cost*dt

# LaTeX

``cost = sign\times dt \sum_{v\in V} c_v``

where ``c_v`` is given by

`` c_v = \sum_{i\in Ix} \frac{y_i - y_{i-1}}{x_i - x_{i-1}} v^{p.w.}_i ``

# Arguments

* psi_container::PSIContainer : the psi_container model built in PowerSimulations
* var_name::Symbol: The variable name
* component_name::String: The component_name of the variable container
* cost_component::PSY.VariableCost{Vector{NTuple{2, Float64}}}
"""
function variable_cost!(
    psi_container::PSIContainer,
    spec::AddCostSpec,
    component_name::String,
    cost_component::PSY.VariableCost{Vector{NTuple{2, Float64}}},
    time_period::Int,
)
    @debug "PWL Variable Cost" component_name
    resolution = model_resolution(psi_container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    # If array is full of tuples with zeros return 0.0
    cost_data = PSY.get_cost(cost_component)
    if all(iszero.(last.(cost_data)))
        @debug "All cost terms for component $(component_name) are 0.0"
        return JuMP.AffExpr(0.0)
    end

    var_name = make_variable_name(spec.variable_type, spec.component_type)
    if !pwlparamcheck(cost_component)
        @warn("The cost function provided for $(var_name) device is not compatible with a linear PWL cost function.
        An SOS-2 formulation will be added to the model.
        This will result in additional binary variables added to the model.")
        gen_cost =
            pwl_gencost_sos!(psi_container, spec, component_name, cost_data, time_period)
    else
        gen_cost =
            pwl_gencost_linear!(psi_container, spec, component_name, cost_data, time_period)
    end
    add_to_cost_expression!(psi_container, spec.multiplier * gen_cost * dt)
    return
end
