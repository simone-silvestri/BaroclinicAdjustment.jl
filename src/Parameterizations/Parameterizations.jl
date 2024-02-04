module Parameterizations

export QGLeith
export Smagorinsky
export OMp25Closure
export GeometricBilaplacian
export Leith
export BiharmonicLeith
export EnergyBackScatter

using Oceananigans
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
        compute_diffusivities!,
        DiffusivityFields,
        viscosity, 
        diffusivity,
        diffusive_flux_x,
        diffusive_flux_y, 
        diffusive_flux_z,
        viscous_flux_ux,
        viscous_flux_vx,
        viscous_flux_uy,
        viscous_flux_vy

using Oceananigans.Utils: launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Operators
using Oceananigans.BuoyancyModels: ∂x_b, ∂y_b, ∂z_b 

using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalFormulation
using Oceananigans.Operators
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz

using Oceananigans.Operators: ℑxyzᶜᶜᶠ, ℑyzᵃᶜᶠ, ℑxzᶜᵃᶠ, Δxᶜᶜᶜ, Δyᶜᶜᶜ

"Return the filter width for a Leith Diffusivity on a general grid."
@inline Δ²ᶜᶜᶜ(i, j, k, grid) =  2 * (1 / (1 / Δxᶜᶜᶜ(i, j, k, grid)^2 + 1 / Δyᶜᶜᶜ(i, j, k, grid)^2))

"The averaged filter width"
@inline Δ̃ᶜᶜᶜ(i, j, k, grid) = sqrt((Δxᶜᶜᶜ(i, j, k, grid)^2 + Δyᶜᶜᶜ(i, j, k, grid)^2)/2)

include("qg_leith_viscosity.jl")
include("geometric_bilaplacian.jl")
include("smagorinsky_laplacian_viscosity.jl")
include("leith_laplacian_viscosity.jl")
include("biharmonic_leith_viscosity.jl")
include("omp25_lateral_friction.jl")
include("energy_backscatter.jl")

end