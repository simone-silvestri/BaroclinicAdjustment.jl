using Oceananigans
using KernelAbstractions: @index, @kernel

using Oceananigans.TurbulenceClosures:
        tapering_factorᶠᶜᶜ,
        tapering_factorᶜᶠᶜ,
        tapering_factorᶜᶜᶠ,
        SmallSlopeIsopycnalTensor,
        AbstractScalarDiffusivity,
        ExplicitTimeDiscretization,
        FluxTapering,
        isopycnal_rotation_tensor_xz_ccf,
        isopycnal_rotation_tensor_yz_ccf,
        isopycnal_rotation_tensor_zz_ccf

import Oceananigans.TurbulenceClosures:
        calculate_diffusivities!,
        DiffusivityFields,
        viscosity, 
        diffusivity,
        diffusive_flux_x,
        diffusive_flux_y, 
        diffusive_flux_z

using Oceananigans.Utils: launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Operators
using Oceananigans.BuoyancyModels: ∂x_b, ∂y_b, ∂z_b 

using Oceananigans.Operators: ℑxyzᶜᶜᶠ, ℑyzᵃᶜᶠ, ℑxyzᶜᶜᶠ, ℑxzᶜᵃᶠ, Δxᶜᶜᶜ, Δyᶜᶜᶜ

struct QGLeithViscosity{A, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation}
    C :: A
    isopycnal_tensor :: M
    slope_limiter :: S
end

QGLeithViscosity(; C=1.0, isopycnal_model=SmallSlopeIsopycnalTensor(), slope_limiter=FluxTapering(1e-2)) =
    QGLeithViscosity(C, isopycnal_model, slope_limiter) 

DiffusivityFields(grid, tracer_names, bcs, ::QGLeithViscosity) = 
                (; νₑ = CenterField(grid), )

@inline function abs²_∇h_ζ(i, j, k, grid, fields)

    ∂xζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)

    return (∂xζ^2 + ∂yζ^2)
end

@inline function abs²_∇h_δ(i, j, k, grid, fields)

    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    
    return (∂xδ^2 + ∂yδ^2)
end

@inline ∂yb_times_f2_div_N2(i, j, k, grid, coriolis, buoyancy, tracers) = ℑxyzᶜᶜᶠ(i, j, k, grid, fᶠᶠᵃ, coriolis)^2 / 
                                                                          max(1e-10, ∂z_b(i, j, k, grid, buoyancy, tracers)) *
                                                                          ℑyzᵃᶜᶠ(i, j, k, grid, ∂y_b, buoyancy, tracers)

@inline ∂xb_times_f2_div_N2(i, j, k, grid, coriolis, buoyancy, tracers) = ℑxyzᶜᶜᶠ(i, j, k, grid, fᶠᶠᵃ, coriolis)^2 / 
                                                                          max(1e-10, ∂z_b(i, j, k, grid, buoyancy, tracers))  *
                                                                          ℑxzᶜᵃᶠ(i, j, k, grid, ∂x_b, buoyancy, tracers)

@inline function abs²_∇h_q(i, j, k, grid, coriolis, buoyancy, tracers)

    ∂zqx = ∂zᶜᶜᶜ(i, j, k, grid, ∂xb_times_f2_div_N2, coriolis, buoyancy, tracers)
    ∂zqy = ∂zᶜᶜᶜ(i, j, k, grid, ∂yb_times_f2_div_N2, coriolis, buoyancy, tracers)
    
    return (∂zqx^2 + ∂zqy^2)
end

@inline Δᶠ(i, j, k, grid) = sqrt(Δxᶜᶜᶜ(i, j, k, grid) * Δyᶜᶜᶜ(i, j, k, grid)) 

@kernel function calculate_qgleith_viscosity!(ν, grid, closure, coriolis, buoyancy, velocities, tracers)
    i, j, k = @index(Global, NTuple)

    ∂ζ² =  abs²_∇h_ζ(i, j, k, grid, velocities)
    ∂δ² =  abs²_∇h_δ(i, j, k, grid, velocities)
    ∂q² =  abs²_∇h_q(i, j, k, grid, coriolis, buoyancy, tracers)

    C = closure.C

    @inbounds ν[i, j, k] = (C * Δᶠ(i, j, k, grid) / π)^3 * sqrt(∂ζ² + ∂δ² + ∂q²) 
end

function calculate_diffusivities!(diffusivity_fields, closure::QGLeithViscosity, model)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    coriolis = model.coriolis

    launch!(arch, grid, :xyz,
            calculate_qgleith_viscosity!,
            diffusivity_fields.νₑ, grid, closure, coriolis, buoyancy, velocities, tracers)

    return nothing
end

"Return the filter width for a Leith Diffusivity on a general grid."
@inline Δᶠ(i, j, k, grid, ::QGLeithViscosity) = sqrt(Δxᶜᶜᶜ(i, j, k, grid) * Δyᶜᶜᶜ(i, j, k, grid)) 


@inline viscosity(::QGLeithViscosity, K) = K.νₑ
@inline diffusivity(::QGLeithViscosity, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

# Diffusive fluxes for Leith diffusivities

@inline diffusive_flux_x(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)

#     νₑ = diffusivities.νₑ

#     νₑⁱʲᵏ = ℑxᶠᵃᵃ(i, j, k, grid, νₑ)
#     ∂x_c  = ∂xᶠᶜᶜ(i, j, k, grid, c)

#     ϵ = tapering_factorᶠᶜᶜ(i, j, k, grid, closure, fields, buoyancy)

#     return - ϵ * νₑⁱʲᵏ * ∂x_c
# end

@inline diffusive_flux_y(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)

#     νₑ = diffusivities.νₑ

#     νₑⁱʲᵏ = ℑyᵃᶠᵃ(i, j, k, grid, νₑ)
#     ∂y_c  = ∂yᶜᶠᶜ(i, j, k, grid, c)

#     ϵ = tapering_factorᶜᶠᶜ(i, j, k, grid, closure, fields, buoyancy)

#     return - ϵ *νₑⁱʲᵏ * ∂y_c
# end

@inline diffusive_flux_z(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)

#     νₑ = diffusivities.νₑ

#     νₑⁱʲᵏ = ℑzᵃᵃᶠ(i, j, k, grid, νₑ)

#     ∂x_c = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, c)
#     ∂y_c = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, c)
#     ∂z_c = ∂zᶜᶜᶠ(i, j, k, grid, c)

#     R₃₁ = isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)
#     R₃₂ = isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)
#     R₃₃ = isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)

#     ϵ = tapering_factorᶜᶜᶠ(i, j, k, grid, closure, fields, buoyancy)

#     return - ϵ * νₑⁱʲᵏ * (
#           2 * R₃₁ * ∂x_c
#         + 2 * R₃₂ * ∂y_c
#             + R₃₃ * ∂z_c)
# end

                
                
