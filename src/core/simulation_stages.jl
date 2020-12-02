######## Internal Simulation Object Structs ########
mutable struct StageInternal
    number::Int
    executions::Int
    execution_count::Int
    end_of_interval_step::Int
    # This line keeps track of the executions of a stage relative to other stages.
    # This might be needed in the future to run multiple stages. For now it is disabled
    # synchronized_executions::Dict{Int, Int} # Number of executions per upper level stage step
    psi_container::PSIContainer
    # Caches are stored in set because order isn't relevant and they should be unique
    caches::Set{CacheKey}
    chronolgy_dict::Dict{Int, <:FeedForwardChronology}
    built::BUILD_STATUS
    stage_path::String
    ext::Dict{String, Any}
    function StageInternal(
        number,
        executions,
        execution_count,
        psi_container;
        ext = Dict{String, Any}(),
    )
        new(
            number,
            executions,
            execution_count,
            0,
            psi_container,
            Set{CacheKey}(),
            Dict{Int, FeedForwardChronology}(),
            EMPTY,
            "",
            ext,
        )
    end
end

# TODO: Add DocString
@doc raw"""
    Stage({M<:AbstractOperationsProblem}
        template::OperationsProblemTemplate
        sys::PSY.System
        optimizer::JuMP.MOI.OptimizerWithAttributes
        internal::Union{Nothing, StageInternal}
        )
"""
mutable struct Stage{M <: AbstractOperationsProblem}
    template::OperationsProblemTemplate
    sys::PSY.System
    internal::Union{Nothing, StageInternal}

    function Stage{M}(
        template::OperationsProblemTemplate,
        sys::PSY.System,
        settings::PSISettings,
        jump_model::Union{Nothing, JuMP.AbstractModel} = nothing,
    ) where {M <: AbstractOperationsProblem}
        internal = StageInternal(0, 0, 0, PSIContainer(sys, settings, jump_model))
        new{M}(template, sys, internal)
    end
end

function Stage{M}(
    template::OperationsProblemTemplate,
    sys::PSY.System,
    optimizer::JuMP.MOI.OptimizerWithAttributes,
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
) where {M <: AbstractOperationsProblem}
    check_kwargs(kwargs, STAGE_ACCEPTED_KWARGS, "Stage")
    settings = PSISettings(sys; optimizer = optimizer, use_parameters = true, kwargs...)
    return Stage{M}(template, sys, settings, jump_model)
end

"""
    Stage(::Type{M},
    template::OperationsProblemTemplate,
    sys::PSY.System,
    optimizer::JuMP.MOI.OptimizerWithAttributes,
    jump_model::Union{Nothing, JuMP.AbstractModel}=nothing;
    kwargs...) where {M<:AbstractOperationsProblem}
This builds the optimization problem of type M with the specific system and template for the simulation stage
# Arguments
- `::Type{M} where M<:AbstractOperationsProblem`: The abstract operation model type
- `template::OperationsProblemTemplate`: The model reference made up of transmission, devices,
                                          branches, and services.
- `sys::PSY.System`: the system created using Power Systems
- `jump_model::Union{Nothing, JuMP.AbstractModel}`: Enables passing a custom JuMP model. Use with care
# Output
- `Stage::Stage`: The operation model containing the model type, unbuilt JuMP model, Power
Systems system.
# Example
```julia
template = OperationsProblemTemplate(CopperPlatePowerModel, devices, branches, services)
stage = Stage(MyOpProblemType template, system, optimizer)
```
# Accepted Key Words
- `initial_time::Dates.DateTime`: Initial Time for the model solve
- `PTDF::PTDF`: Passes the PTDF matrix into the optimization model for StandardPTDFModel networks.
- `warm_start::Bool` True will use the current operation point in the system to initialize variable values. False initializes all variables to zero. Default is true
- `balance_slack_variables::Bool` True will add slacks to the system balance constraints
- `services_slack_variables::Bool` True will add slacks to the services requirement constraints
- `export_pwl_vars::Bool` True will write the results of the piece-wise-linear intermediate variables. Slows down the simulation process significantly
- `allow_fails::Bool`  True will allow the simulation to continue if the optimizer can't find a solution. Use with care, can lead to unwanted behaviour or results
"""
function Stage(
    ::Type{M},
    template::OperationsProblemTemplate,
    sys::PSY.System,
    optimizer::JuMP.MOI.OptimizerWithAttributes,
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
) where {M <: AbstractOperationsProblem}
    return Stage{M}(template, sys, optimizer, jump_model; kwargs...)
end

function Stage(
    template::OperationsProblemTemplate,
    sys::PSY.System,
    optimizer::JuMP.MOI.OptimizerWithAttributes,
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
)
    return Stage{GenericOpProblem}(template, sys, optimizer, jump_model; kwargs...)
end

stage_built(s::Stage) = s.internal.built == BUILT
stage_empty(s::Stage) = s.internal.built == EMPTY
get_execution_count(s::Stage) = s.internal.execution_count
get_executions(s::Stage) = s.internal.executions
get_sys(s::Stage) = s.sys
get_template(s::Stage) = s.template
get_number(s::Stage) = s.internal.number
get_psi_container(s::Stage) = s.internal.psi_container
get_end_of_interval_step(s::Stage) = s.internal.end_of_interval_step
warm_start_enabled(s::Stage) = get_warm_start(s.internal.psi_container.settings)
get_initial_time(s::Stage{T}) where {T <: AbstractOperationsProblem} =
    get_initial_time(s.internal.psi_container.settings)
get_resolution(s::Stage) = IS.time_period_conversion(PSY.get_time_series_resolution(s.sys))
get_settings(s::Stage) = get_psi_container(s).settings

function reset!(stage::Stage{M}) where {M <: AbstractOperationsProblem}
    @assert stage_built(stage)
    if stage_built(stage)
        @info("Stage $(stage.internal.number) will be reset by the simulation build call")
    end
    stage.internal.execution_count = 0
    stage.internal.psi_container =
        PSIContainer(stage.sys, stage.internal.psi_container.settings, nothing)
    stage.internal.built = EMPTY
    return
end

function build!(
    stage::Stage{M},
    initial_time::Dates.DateTime,
    horizon::Int,
    stage_interval::Dates.Period,
) where {M <: PowerSimulationsOperationsProblem}
    !stage_empty(stage) && reset!(stage)
    settings = get_settings(get_psi_container(stage))
    # Horizon and initial time are set here because the information is specified in the
    # Simulation Sequence object and not at the stage creation.
    set_horizon!(settings, horizon)
    set_initial_time!(settings, initial_time)
    stage.internal.built = IN_PROGRESS
    psi_container = get_psi_container(stage)
    # TODO: Abstract the code to just require implementation of _build(). The user shouldn't need
    # to re-implement all the code in this function
    _build!(psi_container, stage.template, stage.sys)
    @assert get_horizon(psi_container.settings) == length(psi_container.time_steps)
    stage_resolution = get_resolution(stage)
    stage.internal.end_of_interval_step = Int(stage_interval / stage_resolution)
    stage_path = stage.internal.stage_path
    _write_psi_container(
        stage.internal.psi_container,
        joinpath(stage_path, "Stage$(stage.internal.number)_optimization_model.json"),
    )
    if get_system_to_file(settings)
        PSY.to_json(
            stage.sys,
            joinpath(stage_path, "Stage$(stage.internal.number)_sys_data.json"),
        )
    end
    stage.internal.built = BUILT
    return
end

function run_stage(
    stage::Stage{M},
    start_time::Dates.DateTime,
    results_path::String,
) where {M <: PowerSimulationsOperationsProblem}
    @assert stage.internal.psi_container.JuMPmodel.moi_backend.state != MOIU.NO_OPTIMIZER
    timed_log = Dict{Symbol, Any}()
    model = stage.internal.psi_container.JuMPmodel
    settings = get_settings(stage)
    _, timed_log[:timed_solve_time], timed_log[:solve_bytes_alloc], timed_log[:sec_in_gc] =
        @timed JuMP.optimize!(model)

    @info "JuMP.optimize! completed" timed_log

    model_status = JuMP.primal_status(stage.internal.psi_container.JuMPmodel)
    if model_status != MOI.FEASIBLE_POINT::MOI.ResultStatusCode
        if settings.allow_fails
            @warn("Stage $(stage.internal.number) status is $(model_status)")
        else
            error("Stage $(stage.internal.number) status is $(model_status)")
        end
    end
    # TODO: Add Fallback when optimization fails
    export_model_result(stage, start_time, results_path)
    export_optimizer_log(timed_log, stage.internal.psi_container, results_path)
    stage.internal.execution_count += 1
    # Reset execution count at the end of step
    if stage.internal.execution_count == stage.internal.executions
        stage.internal.execution_count = 0
    end
    return
end

# Here because requires the stage to be defined
# This is a method a user defining a custom cache will have to define. This is the definition
# in PSI for the building the TimeStatusChange
function get_initial_cache(cache::AbstractCache, stage::Stage)
    throw(ArgumentError("Initialization method for cache $(typeof(cache)) not defined"))
end

function get_initial_cache(cache::TimeStatusChange, stage::Stage)
    ini_cond_on = get_initial_conditions(
        stage.internal.psi_container,
        TimeDurationON,
        cache.device_type,
    )

    ini_cond_off = get_initial_conditions(
        stage.internal.psi_container,
        TimeDurationOFF,
        cache.device_type,
    )

    device_axes = Set((
        PSY.get_name(ic.device) for ic in Iterators.Flatten([ini_cond_on, ini_cond_off])
    ),)
    value_array = JuMP.Containers.DenseAxisArray{Dict{Symbol, Any}}(undef, device_axes)

    for ic in ini_cond_on
        device_name = PSY.get_name(ic.device)
        condition = get_condition(ic)
        status = (condition > 0.0) ? 1.0 : 0.0
        value_array[device_name] = Dict(:count => condition, :status => status)
    end

    for ic in ini_cond_off
        device_name = PSY.get_name(ic.device)
        condition = get_condition(ic)
        status = (condition > 0.0) ? 0.0 : 1.0
        if value_array[device_name][:status] != status
            throw(IS.ConflictingInputsError("Initial Conditions for $(device_name) are not compatible. The values provided are invalid"))
        end
    end

    return value_array
end

function get_initial_cache(cache::StoredEnergy, stage::Stage)
    ini_cond_level =
        get_initial_conditions(stage.internal.psi_container, EnergyLevel, cache.device_type)

    device_axes = Set([PSY.get_name(ic.device) for ic in ini_cond_level],)
    value_array = JuMP.Containers.DenseAxisArray{Float64}(undef, device_axes)
    for ic in ini_cond_level
        device_name = PSY.get_name(ic.device)
        condition = get_condition(ic)
        value_array[device_name] = condition
    end
    return value_array
end

function get_timestamps(stage::Stage, start_time::Dates.DateTime)
    resolution = get_resolution(stage)
    horizon = stage.internal.psi_container.time_steps[end]
    range_time = collect(start_time:resolution:(start_time + resolution * horizon))
    time_stamp = DataFrames.DataFrame(Range = range_time[:, 1])

    return time_stamp
end

function write_data(stage::Stage, save_path::AbstractString; kwargs...)
    write_data(stage.internal.psi_container, save_path; kwargs...)
    return
end

# These functions are writing directly to the feather file and skipping printing to memory.
function export_model_result(stage::Stage, start_time::Dates.DateTime, save_path::String)
    duals = Dict()
    if is_milp(stage.internal.psi_container)
        @warn("Stage $(stage.internal.number) is an MILP, duals can't be exported")
    else
        for c in get_constraint_duals(get_psi_container(stage).settings)
            v = get_constraint(get_psi_container(stage), c)
            duals[c] = axis_array_to_dataframe(v)
        end
    end
    write_data(stage, save_path)
    write_data(duals, save_path; duals = true)
    write_data(get_parameters_value(stage.internal.psi_container), save_path; params = true)
    write_data(get_timestamps(stage, start_time), save_path, "time_stamp")
    files = collect(readdir(save_path))
    compute_file_hash(save_path, files)
    return
end

struct StageSerializationWrapper
    template::OperationsProblemTemplate
    sys::String
    settings::PSISettings
    stage_type::DataType
end
