abstract type AbstractDCLineFormulation <: AbstractBranchFormulation end
struct HVDCUnbounded <: AbstractDCLineFormulation end
struct HVDCLossless <: AbstractDCLineFormulation end
struct HVDCDispatch <: AbstractDCLineFormulation end
struct VoltageSourceDC <: AbstractDCLineFormulation end

const UnboundedBranches =
    Union{HVDCUnbounded, StaticLineUnbounded, StaticTransformerUnbounded}

#################################### Branch Variables ##################################################
flow_variables!(
    psi_container::PSIContainer,
    ::Type{<:PM.AbstractPowerModel},
    devices::IS.FlattenIteratorWrapper{<:PSY.DCBranch},
) = nothing

function add_variables!(
    psi_container::PSIContainer,
    ::StandardPTDFModel,
    devices::IS.FlattenIteratorWrapper{B},
) where {B <: PSY.DCBranch}
    time_steps = model_time_steps(psi_container)
    var_name = make_variable_name(FLOW_ACTIVE_POWER, B)
    container = container_spec(
        psi_container.JuMPmodel,
        [PSY.get_name(d) for d in devices],
        time_steps,
    )
    assign_variable!(psi_container, FLOW_ACTIVE_POWER, B, container)
    for (ix, d) in enumerate(devices)
        bus_fr = PSY.get_number(PSY.get_arc(d).from)
        bus_to = PSY.get_number(PSY.get_arc(d).to)
        for t in time_steps
            jvariable = JuMP.@variable(
                psi_container.JuMPmodel,
                base_name = "$(bus_fr), $(bus_to)_{$(PSY.get_name(d)), $(t)}",
            )
            container[PSY.get_name(d), t] = jvariable
            add_to_expression!(
                psi_container.expressions[:nodal_balance_active],
                PSY.get_number(PSY.get_arc(d).from),
                t,
                jvariable,
                -1.0,
            )
            add_to_expression!(
                psi_container.expressions[:nodal_balance_active],
                PSY.get_number(PSY.get_arc(d).to),
                t,
                jvariable,
                1.0,
            )
        end
    end
    return
end

#################################### Flow Variable Bounds ##################################################
#################################### Rate Limits Constraints ##################################################
function branch_rate_constraints!(
    psi_container::PSIContainer,
    devices::IS.FlattenIteratorWrapper{B},
    model::DeviceModel{B, <:AbstractDCLineFormulation},
    ::Type{<:PM.AbstractDCPModel},
    feedforward::Union{Nothing, AbstractAffectFeedForward},
) where {B <: PSY.DCBranch}
    var = get_variable(psi_container, FLOW_ACTIVE_POWER, B)
    time_steps = model_time_steps(psi_container)
    constraint_val =
        JuMPConstraintArray(undef, [PSY.get_name(d) for d in devices], time_steps)
    assign_constraint!(psi_container, FLOW_ACTIVE_POWER, B, constraint_val)
    for t in time_steps, d in devices
        min_rate = max(
            PSY.get_active_power_limits_from(d).min,
            PSY.get_active_power_limits_to(d).min,
        )
        max_rate = min(
            PSY.get_active_power_limits_from(d).max,
            PSY.get_active_power_limits_to(d).max,
        )
        constraint_val[PSY.get_name(d), t] = JuMP.@constraint(
            psi_container.JuMPmodel,
            min_rate <= var[PSY.get_name(d), t] <= max_rate
        )
    end
    return
end

function branch_rate_constraints!(
    psi_container::PSIContainer,
    devices::IS.FlattenIteratorWrapper{B},
    model::DeviceModel{B, HVDCLossless},
    ::Type{<:PM.AbstractPowerModel},
    feedforward::Union{Nothing, AbstractAffectFeedForward},
) where {B <: PSY.DCBranch}
    for (var_type, cons_type) in zip(
        (FLOW_ACTIVE_POWER_FROM_TO, FLOW_ACTIVE_POWER_TO_FROM),
        (RATE_LIMIT_FT, RATE_LIMIT_TF),
    )
        var = get_variable(psi_container, var_type, B)
        time_steps = model_time_steps(psi_container)
        constraint_val =
            JuMPConstraintArray(undef, [PSY.get_name(d) for d in devices], time_steps)
        assign_constraint!(psi_container, cons_type, B, constraint_val)
        time_steps = model_time_steps(psi_container)

        for t in time_steps, d in devices
            min_rate = max(
                PSY.get_active_power_limits_from(d).min,
                PSY.get_active_power_limits_to(d).min,
            )
            max_rate = min(
                PSY.get_active_power_limits_from(d).max,
                PSY.get_active_power_limits_to(d).max,
            )
            name = PSY.get_name(d)
            constraint_val[name, t] = JuMP.@constraint(
                psi_container.JuMPmodel,
                min_rate <= var[name, t] <= max_rate
            )
        end
    end
    return
end

function branch_rate_constraints!(
    psi_container::PSIContainer,
    devices::IS.FlattenIteratorWrapper{B},
    model::DeviceModel{B, <:AbstractDCLineFormulation},
    ::Union{
        Type{<:PM.AbstractActivePowerModel},
        Type{<:PM.AbstractPowerModel},
    },
    feedforward::Union{Nothing, AbstractAffectFeedForward},
) where {B <: PSY.DCBranch}
    time_steps = model_time_steps(psi_container)
    for (var_type, cons_type) in zip(
        (FLOW_ACTIVE_POWER_FROM_TO, FLOW_ACTIVE_POWER_TO_FROM),
        (RATE_LIMIT_FT, RATE_LIMIT_TF),
    )
        var = get_variable(psi_container, var_type, B)
        constraint_val =
            JuMPConstraintArray(undef, [PSY.get_name(d) for d in devices], time_steps)
        assign_constraint!(psi_container, cons_type, B, constraint_val)

        for t in time_steps, d in devices
            min_rate = max(
                PSY.get_active_power_limits_from(d).min,
                PSY.get_active_power_limits_to(d).min,
            )
            max_rate = min(
                PSY.get_active_power_limits_from(d).max,
                PSY.get_active_power_limits_to(d).max,
            )
            constraint_val[PSY.get_name(d), t] = JuMP.@constraint(
                psi_container.JuMPmodel,
                min_rate <= var[PSY.get_name(d), t] <= max_rate
            )
            add_to_expression!(
                psi_container.expressions[:nodal_balance_active],
                PSY.get_number(PSY.get_arc(d).to),
                t,
                var[PSY.get_name(d), t],
                -PSY.get_loss(d).l1,
                -PSY.get_loss(d).l0,
            )
        end
    end
    return
end
