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

using Oceananigans.Operators: ℑxyzᶜᶜᶠ, ℑyzᵃᶜᶠ, ℑxzᶜᵃᶠ, Δxᶜᶜᶜ, Δyᶜᶜᶜ

struct QGLeithViscosity{A, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation}
    C :: A
    isopycnal_tensor :: M
    slope_limiter :: S
end

QGLeithViscosity(; C=1.0, isopycnal_model=SmallSlopeIsopycnalTensor(), slope_limiter=FluxTapering(1e-2)) =
    QGLeithViscosity(C, isopycnal_model, slope_limiter) 

DiffusivityFields(grid, tracer_names, bcs, ::QGLeithViscosity) = 
                (; νₑ = CenterField(grid),
                   Ld = Field{Center, Center, Nothing}(grid))

@inline function abs²_∇h_ζ(i, j, k, grid, coriolis, fields)

    ∂xζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)

    ∂xf = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, fᶠᶠᵃ, coriolis)
    ∂yf = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, fᶠᶠᵃ, coriolis)
    
    return (∂xζ^2 + ∂yζ^2 + ∂xf^2 + ∂yf^2)
end

@inline function abs²_∇h_δ(i, j, k, grid, fields)

    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    
    return (∂xδ^2 + ∂yδ^2)
end

@inline ∂yb_times_f2_div_N2(i, j, k, grid, coriolis, buoyancy, tracers) = ℑxyzᶜᶜᶠ(i, j, k, grid, fᶠᶠᵃ, coriolis)^2 / 
                                                                          max(1e-20, ∂z_b(i, j, k, grid, buoyancy, tracers)) *
                                                                          ℑyzᵃᶜᶠ(i, j, k, grid, ∂y_b, buoyancy, tracers)

@inline ∂xb_times_f2_div_N2(i, j, k, grid, coriolis, buoyancy, tracers) = ℑxyzᶜᶜᶠ(i, j, k, grid, fᶠᶠᵃ, coriolis)^2 / 
                                                                          max(1e-20, ∂z_b(i, j, k, grid, buoyancy, tracers))  *
                                                                          ℑxzᶜᵃᶠ(i, j, k, grid, ∂x_b, buoyancy, tracers)

@inline function abs²_∇h_q(i, j, k, grid, coriolis, buoyancy, tracers)

    ∂zqx = ∂zᶜᶜᶜ(i, j, k, grid, ∂xb_times_f2_div_N2, coriolis, buoyancy, tracers)
    ∂zqy = ∂zᶜᶜᶜ(i, j, k, grid, ∂yb_times_f2_div_N2, coriolis, buoyancy, tracers)

    return (∂zqx^2 + ∂zqy^2)
end

"Return the filter width for a Leith Diffusivity on a general grid."
@inline Δ²ᶜᶜᶜ(i, j, k, grid) =  2 * (1 / (1 / Δxᶜᶜᶜ(i, j, k, grid)^2 + 1 / Δyᶜᶜᶜ(i, j, k, grid)^2))

@kernel function calculate_qgleith_viscosity!(ν, grid, Ld, closure, velocities, tracers, buoyancy, coriolis)
    i, j, k = @index(Global, NTuple)

    ∂ζ² =  abs²_∇h_ζ(i, j, k, grid, coriolis, velocities)
    ∂δ² =  abs²_∇h_δ(i, j, k, grid, velocities)
    ∂q² =  abs²_∇h_q(i, j, k, grid, coriolis, buoyancy, tracers)

    A  = Δ²ᶜᶜᶜ(i, j, k, grid)

    Bu  = A / Ld[i, j, 1]^2
    ∂Q² = min(∂ζ² + ∂q², ∂ζ² * (1 + 1 / Bu))

    C = closure.C

    @inbounds ν[i, j, k] = (C *  A/ π)^(3 / 2) * sqrt(∂Q² + ∂δ²) 
end

@inline _deformation_radius(i, j, k, grid, C, buoyancy, coriolis) = sqrt(max(0, ∂z_b(i, j, k, grid, buoyancy, C))) / π /
                                                                         abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis))

@kernel function calculate_deformation_radius!(Ld, grid, tracers, buoyancy, coriolis)
    i, j = @index(Global, NTuple)

    @inbounds begin
        Ld[i, j, 1] = 0

        @unroll for k in 1:grid.Nz
            Ld[i, j, 1] += Δzᶜᶜᶠ(i, j, k, grid) * _deformation_radius(i, j, k, grid, tracers, buoyancy, coriolis)
        end
    end
end

function calculate_diffusivities!(diffusivity_fields, closure::QGLeithViscosity, model)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    coriolis = model.coriolis

    launch!(arch, grid, :xy, 
            calculate_deformation_radius!, diffusivity_fields.Ld, grid, tracers, buoyancy, coriolis)

    launch!(arch, grid, :xyz,
            calculate_qgleith_viscosity!,
            diffusivity_fields.νₑ, diffusivity_fields.Ld, grid, closure, velocities, tracers, buoyancy, coriolis)

    return nothing
end

@inline viscosity(::QGLeithViscosity, K) = K.νₑ
@inline diffusivity(::QGLeithViscosity, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

# Diffusive fluxes for Leith diffusivities
@inline diffusive_flux_x(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::QGLeithViscosity, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
