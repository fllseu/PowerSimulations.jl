"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, D},
    ::Type{S},
) where {T <: PSY.ThermalGen, D <: AbstractThermalFormulation, S <: PM.AbstractPowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)
    add_variables!(psi_container, OnVariable, devices)
    add_variables!(psi_container, StartVariable, devices)
    add_variables!(psi_container, StopVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, D)

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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, D},
    ::Type{S},
) where {
    T <: PSY.ThermalGen,
    D <: AbstractThermalFormulation,
    S <: PM.AbstractActivePowerModel,
}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, OnVariable, devices)
    add_variables!(psi_container, StartVariable, devices)
    add_variables!(psi_container, StopVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, D)

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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, ThermalBasicUnitCommitment},
    ::Type{S},
) where {T <: PSY.ThermalGen, S <: PM.AbstractPowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)
    add_variables!(psi_container, OnVariable, devices)
    add_variables!(psi_container, StartVariable, devices)
    add_variables!(psi_container, StopVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

    # Constraints
    # TODO: active_power_constraints
    # TODO: refactor constraints such that ALL variables for all devices are added first, and then the constraint creation is trigged
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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, ThermalBasicUnitCommitment},
    ::Type{S},
) where {T <: PSY.ThermalGen, S <: PM.AbstractActivePowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, OnVariable, devices)
    add_variables!(psi_container, StartVariable, devices)
    add_variables!(psi_container, StopVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, ThermalRampLimited},
    ::Type{S},
) where {T <: PSY.ThermalGen, S <: PM.AbstractPowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

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
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

"""
This function creates the model for a full thermal dispatch formulation depending on combination of devices, device_formulation and system_formulation
"""
function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, ThermalRampLimited},
    ::Type{S},
) where {T <: PSY.ThermalGen, S <: PM.AbstractActivePowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

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
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))

    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, D},
    ::Type{S},
) where {
    T <: PSY.ThermalGen,
    D <: AbstractThermalDispatchFormulation,
    S <: PM.AbstractPowerModel,
}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)

    # Initial Conditions

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
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, D},
    ::Type{S},
) where {
    T <: PSY.ThermalGen,
    D <: AbstractThermalDispatchFormulation,
    S <: PM.AbstractActivePowerModel,
}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)

    # Initial Conditions

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
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{T, FixedOutput},
    ::Type{S},
) where {T <: PSY.ThermalGen, S <: PM.AbstractActivePowerModel}
    devices = get_available_components(T, sys)

    if !validate_available_devices(T, devices)
        return
    end

    nodal_expression!(psi_container, devices, S)

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{PSY.ThermalMultiStart, ThermalMultiStartUnitCommitment},
    ::Type{S};
    kwargs...,
) where {S <: PM.AbstractPowerModel}
    devices = PSY.get_components(PSY.ThermalMultiStart, sys)

    if !validate_available_devices(PSY.ThermalMultiStart, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    add_variables!(psi_container, ReactivePowerVariable, devices)
    commitment_variables!(psi_container, devices)
    add_variables!(psi_container, ColdStartVariable, devices)
    add_variables!(psi_container, WarmStartVariable, devices)
    add_variables!(psi_container, HotStartVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_type_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_initial_condition_constraints!(
        psi_container,
        devices,
        model,
        S,
        get_feedforward(model),
    )
    must_run_constraints!(psi_container, devices, model, S, get_feedforward(model))
    initial_range_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))
    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end

function construct_device!(
    psi_container::PSIContainer,
    sys::PSY.System,
    model::DeviceModel{PSY.ThermalMultiStart, ThermalMultiStartUnitCommitment},
    ::Type{S};
    kwargs...,
) where {S <: PM.AbstractActivePowerModel}
    devices = PSY.get_components(PSY.ThermalMultiStart, sys)

    if !validate_available_devices(PSY.ThermalMultiStart, devices)
        return
    end

    # Variables
    add_variables!(psi_container, ActivePowerVariable, devices)
    commitment_variables!(psi_container, devices)
    add_variables!(psi_container, ColdStartVariable, devices)
    add_variables!(psi_container, WarmStartVariable, devices)
    add_variables!(psi_container, HotStartVariable, devices)

    # Initial Conditions
    initial_conditions!(psi_container, devices, model.formulation)

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
    commitment_constraints!(psi_container, devices, model, S, get_feedforward(model))
    ramp_constraints!(psi_container, devices, model, S, get_feedforward(model))
    time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_time_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_type_constraints!(psi_container, devices, model, S, get_feedforward(model))
    startup_initial_condition_constraints!(
        psi_container,
        devices,
        model,
        S,
        get_feedforward(model),
    )
    must_run_constraints!(psi_container, devices, model, S, get_feedforward(model))
    initial_range_constraints!(psi_container, devices, model, S, get_feedforward(model))
    feedforward!(psi_container, devices, model, get_feedforward(model))
    # Cost Function
    cost_function!(psi_container, devices, model, S, get_feedforward(model))

    return
end
