using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Advection: div_Uc
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.TurbulenceClosures: 
                    biharmonic_mask_x, 
                    biharmonic_mask_y, 
                    ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ,
                    ∂x_∇²h_cᶠᶜᶜ, 
                    ∂y_∇²h_cᶜᶠᶜ

struct EnergyBackScatter{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C₄ :: FT
    C₂ :: FT
    Cᴰ :: FT
    Uᵇ :: FT
    τ  :: FT
end

const EBS = EnergyBackScatter

EnergyBackScatter(FT::DataType = Float64; 
                  C₄ = FT(0.06), 
                  C₂ = FT(-0.6), 
                  Cᴰ = FT(3e-2),
                  Uᵇ = FT(0.1),
                  τ  = FT(1 / 45days)) = EnergyBackScatter(C₄, C₂, Cᴰ, Uᵇ, τ) 

DiffusivityFields(grid, tracer_names, bcs, ::EBS) = 
                    (; ν₂ = CenterField(grid),
                       ν₄ = CenterField(grid), 
                       e  = Field((Center, Center, Nothing), grid),
                       G⁻ = Field((Center, Center, Nothing), grid),
    previous_compute_time = Ref(zero(grid)))

# For the moment just one!
@inline shape_function(i, j, k, grid) = one(grid)

@kernel function _calculate_ebs_viscosities!(ν₂, ν₄, e, grid, closure, velocities, coriolis)
    i, j, k = @index(Global, NTuple)
    u, v, _ = velocities
    Δ = min(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

    ϕ = shape_function(i, j, k, grid)
    C₂ = closure.C₂

    δ₁ = Dₛ(i, j, k, grid, u, v)    
    δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, u, v)    

    D_abs = sqrt(δ₁^2 + δ₂^2)
    
    C₄ = closure.C₄
    τ  = closure.τ
    u⁺ = @inbounds sqrt(2 * max(0, e[i, j, 1]))
    Lᵦ = sqrt(u⁺) * sqrt(1 / abs(∂yᶠᶜᶜ(i, j, k, grid, fᶠᶠᵃ, coriolis)))
    
    Lᵐⁱˣ = min(Δ, Lᵦ)

    ν₂_unbounded = C₂ * u⁺ * Lᵐⁱˣ * ϕ^2 

    ∂ˣu = ∂xᶜᶜᶜ(i, j, k, grid, u)
    ∂ʸv = ∂yᶜᶜᶜ(i, j, k, grid, v)
    
    ∂ˣv = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)
    ∂ʸu = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᶠᶠᶜ, u)

    bound = @inbounds 2 * e[i, j, 1] / sqrt(∂ˣu^2 + 0.5 * (∂ˣv + ∂ʸu)^2 + ∂ʸv^2)

    @inbounds ν₂[i, j, k] = max(ν₂_unbounded, - bound)
    @inbounds ν₄[i, j, k] = (C₄ * D_abs + τ) * Δ^4 
end

function compute_diffusivities!(diffusivity_fields, closure::EBS, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    advection = model.advection.b
    coriolis = model.coriolis
    Δt = model.clock.time - diffusivity_fields.previous_compute_time[]
    diffusivity_fields.previous_compute_time[] = model.clock.time

    launch!(arch, grid, :xy, 
            _advance_tke!, diffusivity_fields,
                           grid, model.clock, velocities, advection, 
                           closure, Δt, Val(grid.Nz))

    launch!(arch, grid, parameters,
            _calculate_ebs_viscosities!,
            diffusivity_fields.ν₂, diffusivity_fields.ν₄, 
            diffusivity_fields.e, grid, closure, velocities, coriolis)

    return nothing
end

using Oceananigans.Operators

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2

@inline ∂u_τˣ₂(i, j, k, grid, u, args...) = ∂yᶠᶠᶜ(i, j, k, grid, u) * viscous_flux_uy(i, j, k, grid, args...)
@inline ∂v_τʸ₁(i, j, k, grid, v, args...) = ∂xᶠᶠᶜ(i, j, k, grid, v) * viscous_flux_vx(i, j, k, grid, args...)

@inline ∂u_τˣ₁(i, j, k, grid, u, args...) = ∂xᶜᶜᶜ(i, j, k, grid, u) * viscous_flux_ux(i, j, k, grid, args...)
@inline ∂v_τʸ₂(i, j, k, grid, v, args...) = ∂yᶜᶜᶜ(i, j, k, grid, v) * viscous_flux_vy(i, j, k, grid, args...)

@kernel function _advance_tke!(diffusivities, grid, clock, velocities, advection, closure, Δt, ::Val{Nz}) where Nz
    i, j  = @index(Global, NTuple)
    u, v, _ = velocities

    εˣ = zero(grid)
    εʸ = zero(grid)

    @unroll for k in 1:Nz
        
        Δz = Δzᶜᶜᶜ(i, j, k, grid) / grid.Lz
        
        args = (closure, diffusivities, clock, velocities, nothing)

        εˣ += Δz * (∂u_τˣ₁(i, j, k, grid, u, args...) + ℑxyᶜᶜᵃ(i, j, k, grid, ∂u_τˣ₂, u, args...))
        εʸ += Δz * (∂v_τʸ₂(i, j, k, grid, v, args...) + ℑxyᶜᶜᵃ(i, j, k, grid, ∂v_τʸ₁, v, args...))
    end

    FT = eltype(grid)
    α = convert(FT, 1.6)
    β = convert(FT, 0.6)

    Cᴰ = closure.Cᴰ

    e  = diffusivities.e
    G⁻ = diffusivities.G⁻

    s² = ℑxᶜᵃᵃ(i, j, 1, grid, ϕ², u) + ℑyᵃᶜᵃ(i, j, 1, grid, ϕ², v)

    implicit_bottom_friction = @inbounds Cᴰ * sqrt(s² + e[i, j, 1] + closure.Uᵇ^2) 

    Gⁿ = - div_Uc(i, j, 1, grid, advection, velocities, e) - (εˣ + εʸ)

    #  - bottom_friction + 
    @inbounds e[i, j, 1]  = (e[i, j, 1] + Δt * (α * Gⁿ - β * G⁻[i, j, 1])) / (1 + Δt * implicit_bottom_friction)
    @inbounds G⁻[i, j, 1] = Gⁿ
end

@inline viscosity(::EBS, K) = K.ν₂
@inline diffusivity(::EBS, K, ::Val{id}) where id = K.ν₂

#######
####### Fluxes! 
#######

using Oceananigans.TurbulenceClosures: δ★ᶜᶜᶜ, ζ★ᶠᶠᶜ
using Oceananigans.Operators: div_xyᶜᶜᶜ, ζ₃ᶠᶠᶜ

@inline function viscous_flux_ux(i, j, k, grid, ::EBS, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * div_xyᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * δ★ᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_vx(i, j, k, grid, ::EBS, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * ζ₃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * ζ★ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_uy(i, j, k, grid, ::EBS, K, clk, fields, b) 
    Δ  = @inbounds + K.ν₂[i, j, k] * ζ₃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds - K.ν₄[i, j, k] * ζ★ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end

@inline function viscous_flux_vy(i, j, k, grid, ::EBS, K, clk, fields, b)
    Δ  = @inbounds - K.ν₂[i, j, k] * div_xyᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    Δ² = @inbounds + K.ν₄[i, j, k] * δ★ᶜᶜᶜ(i, j, k, grid, fields.u, fields.v)
    return Δ + Δ²
end
 
@inline diffusive_flux_x(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
