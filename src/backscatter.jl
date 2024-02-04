using Oceananigans
using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.TurbulenceClosures:
        tapering_factorᶠᶜᶜ,
        tapering_factorᶜᶠᶜ,
        tapering_factorᶜᶜᶠ,
        tapering_factor,
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

struct BackScatter{FT, P} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    Cᴮˢ :: FT
    ϕ   :: P
end

DiffusivityFields(grid, tracer_names, bcs, ::BackScatter) = (; νₑ = CenterField(grid))


"Return the filter width for a Leith Diffusivity on a general grid."
@inline Δᶜᶜᶜ(i, j, k, grid) =  sqrt(2 * (1 / (1 / Δxᶜᶜᶜ(i, j, k, grid)^2 + 1 / Δyᶜᶜᶜ(i, j, k, grid)^2)))

function calculate_diffusivities!(diffusivity_fields, closure::BackScatter, model; kwargs...)
    arch = model.architecture
    grid = model.grid
    tracers = model.tracers
    coriolis = model.coriolis

    launch!(arch, grid, :xyz,
            _compute_backscatter_viscosity!,
            diffusivity_fields.νₑ, grid, closure, coriolis, tracers)

    return nothing
end 

@kernel function _compute_backscatter_viscosity!(ν, grid, closure, coriolis, tracers)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ϕ = closure.ϕ[k]
        e = tracers.e[i, j, k]
        Δ = Δᶜᶜᶜ(i, j, k, grid) 
        β = abs(∂yᶜᶜᶜ(i, j, k, grid, ℑxᶜᵃᵃ, fᶠᶠᵃ, coriolis))
        β⁻¹ = min(1e20, 1 / β)
        R = sqrt(sqrt(2e) * β⁻¹)
        L = min(Δ, R)
        ν[i, j, k] = closure.C * ϕ * sqrt(e) * L
    end
end

ImplicitBackScatter()