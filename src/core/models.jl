
### Services

abstract type AbstractServiceFormulation end

mutable struct ServiceModel{D <: PSY.Service, B <: AbstractServiceFormulation}
    service_type::Type{D}
    formulation::Type{B}
end


abstract type AbstractDeviceFormulation end
struct FixedOutput <: AbstractDeviceFormulation end

function _check_device_formulation(
    ::Type{D},
) where {D <: Union{AbstractDeviceFormulation, PSY.Device}}
    if !isconcretetype(D)
        throw(ArgumentError("The device model must contain only concrete types, $(D) is an Abstract Type"))
    end
end

### Device models

"""
    DeviceModel(::Type{D}, ::Type{B}) where {D<:PSY.Device,
                                       B<:AbstractDeviceFormulation}
This validates the device formulation for the Power System Device and the
abstract device formulation and returns  Power System Device and the
abstract device formulation if the power system device is a concrete type.

# Arguments
-`::Type{D}`: Power System Device
-`::Type{B}`: Abstract Device Formulation

# Outputs
`DeviceModel(D, B, nothing)`: D::PSY.Device, B::AbstractDeviceFormulation

# Example
```julia
branches = Dict{Symbol, DeviceModel}
    (:L => DeviceModel(PSY.Line, StaticLine),
    :T => DeviceModel(PSY.Transformer2W, StaticTransformer),
    :TT => DeviceModel(PSY.TapTransformer , StaticTransformer),
    :dc_line => DeviceModel(PSY.HVDCLine, HVDCDispatch))
```
"""
mutable struct DeviceModel{D <: PSY.Device, B <: AbstractDeviceFormulation}
    device_type::Type{D}
    formulation::Type{B}
    # TODO: Needs to be made into an array if more than one feedforward is desired
    feedforward::Union{Nothing, AbstractAffectFeedForward}
    services::Vector{ServiceModel}

    function DeviceModel(
        ::Type{D},
        ::Type{B},
        FF = nothing,
    ) where {D <: PSY.Device, B <: AbstractDeviceFormulation}
        _check_device_formulation(D)
        _check_device_formulation(B)
        new{D, B}(D, B, FF, Vector{ServiceModel}())
    end
end

get_device_type(m::DeviceModel) = m.device_type
get_formulation(m::DeviceModel) = m.formulation
get_feedforward(m::DeviceModel) = m.feedforward
get_services(m::Union{DeviceModel, Nothing}) = isnothing(m) ? nothing : m.services

### Network models

struct CopperPlatePowerModel <: PM.AbstractActivePowerModel end
struct AreaBalancePowerModel <: PM.AbstractActivePowerModel end
struct StandardPTDFModel <: PM.AbstractDCPModel end

#================================================
    # exact non-convex models
    ACPPowerModel, ACRPowerModel, ACTPowerModel

    # linear approximations
    DCPPowerModel, NFAPowerModel

    # quadratic approximations
    DCPLLPowerModel, LPACCPowerModel

    # quadratic relaxations
    SOCWRPowerModel, SOCWRConicPowerModel,
    SOCBFPowerModel, SOCBFConicPowerModel,
    QCRMPowerModel, QCLSPowerModel,

    # sdp relaxations
    SDPWRMPowerModel, SparseSDPWRMPowerModel
================================================#

##### Exact Non-Convex Models #####
import PowerModels: ACPPowerModel

import PowerModels: ACRPowerModel

import PowerModels: ACTPowerModel

##### Linear Approximations #####
import PowerModels: DCPPowerModel

import PowerModels: NFAPowerModel

##### Quadratic Approximations #####
import PowerModels: DCPLLPowerModel

import PowerModels: LPACCPowerModel

##### Quadratic Relaxations #####
import PowerModels: SOCWRPowerModel

import PowerModels: SOCWRConicPowerModel

import PowerModels: QCRMPowerModel

import PowerModels: QCLSPowerModel
