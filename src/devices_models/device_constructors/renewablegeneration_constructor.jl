function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{R, D},
    ::Type{S},
) where {
    R <: PSY.RenewableGen,
    D <: AbstractRenewableDispatchFormulation,
    S <: PM.AbstractPowerModel,
}
    devices = get_available_components(R, sys)

    if !validate_available_devices(R, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)

    # Constraints
    add_constraints!(
        psi_container,
        RangeConstraint,
        ActivePowerVariable,
        devices,
        model,
        S,
        get_feedforward(model),
    )
    add_constraints!(
        psi_container,
        RangeConstraint,
        ReactivePowerVariable,
        devices,
        model,
        S,
        get_feedforward(model),
    )
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S)

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{R, D},
    ::Type{S},
) where {
    R <: PSY.RenewableGen,
    D <: AbstractRenewableDispatchFormulation,
    S <: PM.AbstractActivePowerModel,
}
    devices = get_available_components(R, sys)

    if !validate_available_devices(R, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)

    # Constraints
    add_constraints!(
        psi_container,
        RangeConstraint,
        ActivePowerVariable,
        devices,
        model,
        S,
        get_feedforward(model),
    )
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S)

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    ::DeviceModel{R, FixedOutput},
    ::Type{S},
) where {R <: PSY.RenewableGen, S <: PM.AbstractPowerModel}
    devices = get_available_components(R, sys)

    if !validate_available_devices(R, devices)
        return
    end

    nodal_expression!(psi_container, devices, S)

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    ::DeviceModel{PSY.RenewableFix, D},
    ::Type{S},
) where {D <: AbstractRenewableDispatchFormulation, S <: PM.AbstractPowerModel}
    @warn("The Formulation $(D) only applies to FormulationControllable Renewable Resources, \n Consider Changing the Device Formulation to FixedOutput")

    construct_device!(psi_container, sys, DeviceModel(PSY.RenewableFix, FixedOutput), S)

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    ::DeviceModel{PSY.RenewableFix, FixedOutput},
    ::Type{S},
) where {S <: PM.AbstractPowerModel}
    devices = get_available_components(PSY.RenewableFix, sys)

    if !validate_available_devices(PSY.RenewableFix, devices)
        return
    end

    nodal_expression!(psi_container, devices, S)

    return
end
