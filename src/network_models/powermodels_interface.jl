#################################################################################
# Comments
#
# - Ideally the net_injection variables would be bounded.  This can be done using an adhoc data model extention
#
#################################################################################
# Model Definitions

""
function instantiate_nip_expr_model(data::Dict{String, Any}, model_constructor; kwargs...)
    return PM.instantiate_model(data, model_constructor, instantiate_nip_expr; kwargs...)
end

""
function instantiate_nip_expr(pm::PM.AbstractPowerModel)
    for (n, network) in PM.nws(pm)
        @assert !PM.ismulticonductor(pm, nw = n)
        PM.variable_bus_voltage(pm, nw = n)
        PM.variable_branch_power(pm, nw = n; bounded = false)
        PM.variable_dcline_power(pm, nw = n)

        PM.constraint_model_voltage(pm, nw = n)

        for i in PM.ids(pm, :ref_buses, nw = n)
            PM.constraint_theta_ref(pm, i, nw = n)
        end

        for i in PM.ids(pm, :bus, nw = n)
            constraint_power_balance_ni_expr(pm, i, nw = n)
        end

        for i in PM.ids(pm, :branch, nw = n)
            PM.constraint_ohms_yt_from(pm, i, nw = n)
            PM.constraint_ohms_yt_to(pm, i, nw = n)

            PM.constraint_voltage_angle_difference(pm, i, nw = n)
        end

        for i in PM.ids(pm, :dcline)
            PM.constraint_dcline_power_losses(pm, i, nw = n)
        end
    end

    return
end

function instantiate_bfp_expr_model(data::Dict{String, Any}, model_constructor; kwargs...)
    return PM.instantiate_model(data, model_constructor, instantiate_bfp_expr; kwargs...)
end

""
function instantiate_bfp_expr(pm::PM.AbstractPowerModel)
    for (n, network) in PM.nws(pm)
        @assert !PM.ismulticonductor(pm, nw = n)
        PM.variable_bus_voltage(pm, nw = n)
        PM.variable_branch_power(pm, nw = n; bounded = false)
        PM.variable_dcline_power(pm, nw = n)

        PM.constraint_model_current(pm, nw = n)

        for i in PM.ids(pm, :ref_buses, nw = n)
            PM.constraint_theta_ref(pm, i, nw = n)
        end

        for i in PM.ids(pm, :bus, nw = n)
            constraint_power_balance_ni_expr(pm, i, nw = n)
        end

        for i in PM.ids(pm, :branch, nw = n)
            PM.constraint_power_losses(pm, i, nw = n)
            PM.constraint_voltage_magnitude_difference(pm, i, nw = n)

            PM.constraint_voltage_angle_difference(pm, i, nw = n)
        end

        for i in PM.ids(pm, :dcline)
            PM.constraint_dcline_power_losses(pm, i, nw = n)
        end
    end

    return
end

function instantiate_vip_expr_model(data::Dict{String, Any}, model_constructor; kwargs...)
    throw(error("VI Models not currently supported"))
end

#################################################################################
# Model Extention Functions

""
function constraint_power_balance_ni_expr(
    pm::PM.AbstractPowerModel,
    i::Int;
    nw::Int = pm.cnw,
)
    if !haskey(PM.con(pm, nw), :power_balance_p)
        PM.con(pm, nw)[:power_balance_p] = Dict{Int, JuMP.ConstraintRef}()
    end
    if !haskey(PM.con(pm, nw), :power_balance_q)
        PM.con(pm, nw)[:power_balance_q] = Dict{Int, JuMP.ConstraintRef}()
    end

    bus = PM.ref(pm, nw, :bus, i)
    bus_arcs = PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = PM.ref(pm, nw, :bus_arcs_dc, i)

    pni_expr = PM.ref(pm, nw, :bus, i, "pni")
    qni_expr = PM.ref(pm, nw, :bus, i, "qni")

    constraint_power_balance_ni_expr(pm, nw, i, bus_arcs, bus_arcs_dc, pni_expr, qni_expr)

    return
end

""
function constraint_power_balance_ni_expr(
    pm::PM.AbstractPowerModel,
    n::Int,
    i::Int,
    bus_arcs,
    bus_arcs_dc,
    pni_expr,
    qni_expr,
)
    p = PM.var(pm, n, :p)
    q = PM.var(pm, n, :q)
    p_dc = PM.var(pm, n, :p_dc)
    q_dc = PM.var(pm, n, :q_dc)

    PM.con(pm, n, :power_balance_p)[i] = JuMP.@constraint(
        pm.model,
        sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == pni_expr
    )
    PM.con(pm, n, :power_balance_q)[i] = JuMP.@constraint(
        pm.model,
        sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == qni_expr
    )

    return
end

""
function constraint_current_balance_ni_expr(
    pm::PM.AbstractPowerModel,
    i::Int;
    nw::Int = pm.cnw,
)
    if !haskey(PM.con(pm, nw), :kcl_cr)
        PM.con(pm, nw)[:kcl_cr] = Dict{Int, JuMP.ConstraintRef}()
    end
    if !haskey(PM.con(pm, nw), :kcl_ci)
        PM.con(pm, nw)[:kcl_ci] = Dict{Int, JuMP.ConstraintRef}()
    end

    bus = PM.ref(pm, nw, :bus, i)
    bus_arcs = PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = PM.ref(pm, nw, :bus_arcs_dc, i)

    pni_expr = PM.ref(pm, nw, :bus, i, "pni")
    qni_expr = PM.ref(pm, nw, :bus, i, "qni")

    constraint_current_balance_ni_expr(pm, nw, i, bus_arcs, bus_arcs_dc, pni_expr, qni_expr)

    return
end

""
function constraint_current_balance_ni_expr(
    pm::PM.AbstractPowerModel,
    n::Int,
    i::Int,
    bus_arcs,
    bus_arcs_dc,
    pni_expr,
    qni_expr,
)
    p = PM.var(pm, n, :p)
    q = PM.var(pm, n, :q)
    p_dc = PM.var(pm, n, :p_dc)
    q_dc = PM.var(pm, n, :q_dc)

    PM.con(pm, n, :power_balance_p)[i] = JuMP.@constraint(
        pm.model,
        sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == pni_expr
    )
    PM.con(pm, n, :power_balance_q)[i] = JuMP.@constraint(
        pm.model,
        sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == qni_expr
    )

    return
end

"active power only models ignore reactive power variables"
function variable_reactive_net_injection(pm::PM.AbstractActivePowerModel; kwargs...)
    return
end

""
function constraint_power_balance_ni_expr(
    pm::PM.AbstractActivePowerModel,
    n::Int,
    i::Int,
    bus_arcs,
    bus_arcs_dc,
    pni_expr,
    qni_expr,
)
    p = PM.var(pm, n, :p)
    p_dc = PM.var(pm, n, :p_dc)

    PM.con(pm, n, :power_balance_p)[i] = JuMP.@constraint(
        pm.model,
        sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == pni_expr
    )

    return
end

""
function powermodels_network!(
    psi_container::PSIContainer,
    system_formulation::Type{S},
    sys::PSY.System,
    instantiate_model = instantiate_nip_expr_model,
) where {S <: PM.AbstractPowerModel}
    time_steps = model_time_steps(psi_container)
    pm_data, PM_map = pass_to_pm(sys, time_steps[end])
    buses = PSY.get_components(PSY.Bus, sys)

    remove_undef!(psi_container.expressions[:nodal_balance_active])
    remove_undef!(psi_container.expressions[:nodal_balance_reactive])

    for t in time_steps, bus in buses
        pm_data["nw"]["$(t)"]["bus"]["$(bus.number)"]["pni"] =
            psi_container.expressions[:nodal_balance_active][bus.number, t]
        pm_data["nw"]["$(t)"]["bus"]["$(bus.number)"]["qni"] =
            psi_container.expressions[:nodal_balance_reactive][bus.number, t]
    end

    psi_container.pm =
        instantiate_model(pm_data, system_formulation, jump_model = psi_container.JuMPmodel)
    psi_container.pm.ext[:PMmap] = PM_map

    return
end

""
function powermodels_network!(
    psi_container::PSIContainer,
    system_formulation::Type{S},
    sys::PSY.System,
    instantiate_model = instantiate_nip_expr_model,
) where {S <: PM.AbstractActivePowerModel}
    time_steps = model_time_steps(psi_container)
    pm_data, PM_map = pass_to_pm(sys, time_steps[end])
    buses = PSY.get_components(PSY.Bus, sys)

    remove_undef!(psi_container.expressions[:nodal_balance_active])

    for t in time_steps, bus in buses
        pm_data["nw"]["$(t)"]["bus"]["$(PSY.get_number(bus))"]["pni"] =
            psi_container.expressions[:nodal_balance_active][PSY.get_number(bus), t]
        # pm_data["nw"]["$(t)"]["bus"]["$(bus.number)"]["qni"] = 0.0
    end

    psi_container.pm =
        instantiate_model(pm_data, system_formulation, jump_model = psi_container.JuMPmodel)
    psi_container.pm.ext[:PMmap] = PM_map

    return
end

#### PM accessor functions ########

function PMvarmap(system_formulation::Type{S}) where {S <: PM.AbstractDCPModel}
    pm_var_map = Dict{Type, Dict{Symbol, Union{String, NamedTuple}}}()

    pm_var_map[PSY.Bus] = Dict(:va => THETA)
    pm_var_map[PSY.ACBranch] = Dict(:p => (from_to = FLOW_ACTIVE_POWER, to_from = nothing))
    pm_var_map[PSY.DCBranch] =
        Dict(:p_dc => (from_to = FLOW_ACTIVE_POWER, to_from = nothing))

    return pm_var_map
end

function PMvarmap(system_formulation::Type{S}) where {S <: PM.AbstractActivePowerModel}
    pm_var_map = Dict{Type, Dict{Symbol, Union{String, NamedTuple}}}()

    pm_var_map[PSY.Bus] = Dict(:va => THETA)
    pm_var_map[PSY.ACBranch] = Dict(
        :p =>
            (from_to = FLOW_ACTIVE_POWER_FROM_TO, to_from = FLOW_ACTIVE_POWER_TO_FROM),
    )
    pm_var_map[PSY.DCBranch] = Dict(
        :p_dc =>
            (from_to = FLOW_ACTIVE_POWER_FROM_TO, to_from = FLOW_ACTIVE_POWER_TO_FROM),
    )

    return pm_var_map
end

function PMvarmap(system_formulation::Type{S}) where {S <: PM.AbstractPowerModel}
    pm_var_map = Dict{Type, Dict{Symbol, Union{String, NamedTuple}}}()

    pm_var_map[PSY.Bus] = Dict(:va => THETA, :vm => VM)
    pm_var_map[PSY.ACBranch] = Dict(
        :p =>
            (from_to = FLOW_ACTIVE_POWER_FROM_TO, to_from = FLOW_ACTIVE_POWER_TO_FROM),
        :q => (
            from_to = FLOW_REACTIVE_POWER_FROM_TO,
            to_from = FLOW_REACTIVE_POWER_TO_FROM,
        ),
    )
    pm_var_map[PSY.DCBranch] = Dict(
        :p_dc =>
            (from_to = FLOW_ACTIVE_POWER_FROM_TO, to_from = FLOW_ACTIVE_POWER_TO_FROM),
        :q_dc => (
            from_to = FLOW_REACTIVE_POWER_FROM_TO,
            to_from = FLOW_REACTIVE_POWER_TO_FROM,
        ),
    )

    return pm_var_map
end

function PMconmap(system_formulation::Type{S}) where {S <: PM.AbstractActivePowerModel}
    pm_con_map = Dict{Type, Dict{Symbol, String}}()

    pm_con_map[PSY.Bus] = Dict(:power_balance_p => NODAL_BALANCE_ACTIVE)
    return pm_con_map
end

function PMconmap(system_formulation::Type{S}) where {S <: PM.AbstractPowerModel}
    pm_con_map = Dict{Type, Dict{Symbol, String}}()

    pm_con_map[PSY.Bus] = Dict(
        :power_balance_p => NODAL_BALANCE_ACTIVE,
        :power_balance_q => NODAL_BALANCE_REACTIVE,
    )
    return pm_con_map
end

function add_pm_var_refs!(
    psi_container::PSIContainer,
    system_formulation::Type{S},
    sys::PSY.System,
) where {S <: PM.AbstractPowerModel}
    time_steps = model_time_steps(psi_container)
    bus_dict = psi_container.pm.ext[:PMmap].bus
    ACbranch_dict = psi_container.pm.ext[:PMmap].arcs
    ACbranch_types = typeof.(values(ACbranch_dict))
    DCbranch_dict = psi_container.pm.ext[:PMmap].arcs_dc
    DCbranch_types = typeof.(values(DCbranch_dict))

    pm_var_names = keys(psi_container.pm.var[:nw][1])

    pm_var_map = PMvarmap(system_formulation)

    for (pm_v, ps_v) in pm_var_map[PSY.Bus]
        if pm_v in pm_var_names
            container = PSI.container_spec(
                psi_container.JuMPmodel,
                [PSY.get_name(b) for b in values(bus_dict)],
                time_steps,
            )
            assign_variable!(psi_container, ps_v, PSY.Bus, container)
            for t in time_steps, (pm_bus, bus) in bus_dict
                name = PSY.get_name(bus)
                container[name, t] = PM.var(psi_container.pm, t, pm_v)[pm_bus] # pm_vars[pm_v][pm_bus]
            end
        end
    end

    add_pm_var_refs!(
        psi_container,
        PSY.ACBranch,
        ACbranch_types,
        ACbranch_dict,
        pm_var_map,
        pm_var_names,
        time_steps,
    )
    add_pm_var_refs!(
        psi_container,
        PSY.DCBranch,
        DCbranch_types,
        DCbranch_dict,
        pm_var_map,
        pm_var_names,
        time_steps,
    )
end

function add_pm_var_refs!(
    psi_container::PSIContainer,
    d_class::Type,
    device_types::Vector,
    pm_map::Dict,
    pm_var_map::Dict,
    pm_var_names::Base.KeySet,
    time_steps::UnitRange{Int},
)
    for d_type in Set(device_types)
        devices = [d for d in pm_map if typeof(d[2]) == d_type]
        for (pm_v, ps_v) in pm_var_map[d_class]
            if pm_v in pm_var_names
                for dir in fieldnames(typeof(ps_v))
                    isnothing(getfield(ps_v, dir)) && continue
                    # TODO: make a better mapping with the var names in the definitions file
                    var_name = getfield(ps_v, dir)
                    container = PSI.container_spec(
                        psi_container.JuMPmodel,
                        [PSY.get_name(d[2]) for d in devices],
                        time_steps,
                    )
                    assign_variable!(psi_container, var_name, d_type, container)
                    for t in time_steps, (pm_d, d) in devices
                        container[PSY.get_name(d), t] =
                            PM.var(psi_container.pm, t, pm_v, getfield(pm_d, dir))
                    end
                end
            end
        end
    end
end

function add_pm_con_refs!(
    psi_container::PSIContainer,
    system_formulation::Type{S},
    sys::PSY.System,
) where {S <: PM.AbstractPowerModel}
    time_steps = model_time_steps(psi_container)
    bus_dict = psi_container.pm.ext[:PMmap].bus

    pm_con_names = [
        k
        for
        k in keys(psi_container.pm.con[:nw][1]) if !isempty(PM.con(psi_container.pm, 1, k))
    ]

    pm_con_map = PMconmap(system_formulation)
    for (pm_v, ps_v) in pm_con_map[PSY.Bus]
        if pm_v in pm_con_names
            container = PSI.add_cons_container!(
                psi_container,
                make_constraint_name(ps_v, PSY.Bus),
                [PSY.get_name(b) for b in values(bus_dict)],
                time_steps,
            )
            for t in time_steps, (pm_bus, bus) in bus_dict
                name = PSY.get_name(bus)
                container[name, t] = PM.con(psi_container.pm, t, pm_v)[pm_bus]
            end
        end
    end
end
