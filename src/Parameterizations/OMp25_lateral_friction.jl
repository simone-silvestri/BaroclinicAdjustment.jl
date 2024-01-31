using Oceananigans.Operators
struct OMp25Closure{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C₄ :: FT
    u₄ :: FT
    C₂ :: FT
    u₂ :: FT
end

const OMP = OMp25Closure

OMp25Closure(FT::DataType = Float64; C₄ = FT(0.06), u₄ = FT(0.01), C₂ = FT(0.15), u₂ = FT(0.01)) = 
        OMp25Closure(C₄, u₄, C₂, u₂) 

DiffusivityFields(grid, tracer_names, bcs, ::OMp25Closure) = 
                (; ν₂ = CenterField(grid),
                   ν₄ = CenterField(grid), 
                   Ld = Field{Center, Center, Nothing}(grid))

@inline Δ̃ᶜᶜᶜ(i, j, k, grid) = sqrt((Δxᶜᶜᶜ(i, j, k, grid)^2 + Δyᶜᶜᶜ(i, j, k, grid)^2)/2)

@kernel function _calculate_omp25_viscosities!(ν₂, ν₄, Ld, grid, closure, velocities)
    i, j, k = @index(Global, NTuple)
    u, v, _ = velocities
    Δ = min(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

    C₂ = closure.C₂
    u₂ = closure.u₂

    δ₁ = Dₛ(i, j, k, grid, u, v)    
    δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, u, v)    

    D_abs = sqrt(δ₁^2 + δ₂^2)

    ν₂_smag = C₂ * Δ^2 * D_abs
    ν₂_stat = u₂ * Δ

    @inbounds Rh = Ld[i, j, 1] / Δ̃ᶜᶜᶜ(i, j, k, grid)
    F₂ = 1 / (1 + 0.25 * Rh)
    C₄ = closure.C₄
    u₄ = closure.u₄

    ν₄_smag = C₄ * Δ^4 * D_abs
    ν₄_stat = u₄ * Δ^3
    
    @inbounds ν₂[i, j, k] = max(ν₂_smag, ν₂_stat) * F₂
    @inbounds ν₄[i, j, k] = max(ν₄_smag, ν₄_stat)
end

function compute_diffusivities!(diffusivity_fields, closure::OMp25Closure, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    coriolis = model.coriolis

    launch!(arch, grid, :xy, 
            calculate_deformation_radius!, diffusivity_fields.Ld, grid, tracers, buoyancy, coriolis)

    launch!(arch, grid, parameters,
            _calculate_omp25_viscosities!,
            diffusivity_fields.ν₂, diffusivity_fields.ν₄, 
            diffusivity_fields.Ld, grid, closure, velocities)

    return nothing
end

@inline viscosity(::OMp25Closure, K) = K.ν₂
@inline diffusivity(::OMp25Closure, K, ::Val{id}) where id = K.ν₂

#######
####### Fluxes! 
#######

using Oceananigans.TurbulenceClosures: δ★ᶜᶜᶜ, ζ★ᶠᶠᶜ
using Oceananigans.Operators: div_xyᶜᶜᶜ, ζ₃ᶠᶠᶜ

@inline function viscous_flux_ux(i, j, k, grid, ::OMP, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * div_xyᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * δ★ᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_vx(i, j, k, grid, ::OMP, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * ζ₃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * ζ★ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_uy(i, j, k, grid, ::OMP, K, clk, fields, b) 
    Δ  = @inbounds + K.ν₂[i, j, k] * ζ₃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds - K.ν₄[i, j, k] * ζ★ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_vy(i, j, k, grid, ::OMP, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * div_xyᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * δ★ᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end
 
@inline diffusive_flux_x(i, j, k, grid, cl::OMP, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, cl::OMP, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, cl::OMP, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
