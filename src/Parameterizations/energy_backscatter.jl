using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Advection: div_Uc, U_dot_∇u, U_dot_∇v 
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.TurbulenceClosures: 
                    biharmonic_mask_x, 
                    biharmonic_mask_y, 
                    ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ,
                    ∂x_∇²h_cᶠᶜᶜ, 
                    ∂y_∇²h_cᶜᶠᶜ
"""
    struct EnergyBackScatter{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}

The `EnergyBackScatter` struct represents a parameterization for energy backscatter in a baroclinic adjustment model.

Fields
=======

- `C₄`: Coefficient C₄ for the biharmonic smagorinky dissipation
- `C₂`: Coefficient C₂ for the backscattering antiviscosity
- `Cᴰ`: Coefficient Cᴰ for the bottom drag dissipation
- `Uᵇ`: Background velocity Uᵇ 
- `τ`: dissipation timescale τ

Reference
=========
Jansen, M., Adcroft, A., Khani, S., & Kong, H. (2019). Toward an energetically consistent, resolution aware 
parameterization of ocean mesoscale eddies. Journal of Advances in Modeling Earth Systems, 
11 (8), 2844-2860. doi: https://doi.org/10.1029/2019MS001750
"""
struct EnergyBackScatter{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C₄ :: FT
    C₂ :: FT
    Cᴰ :: FT
    Uᵇ :: FT
    τ  :: FT
    implicit_dissipation :: Bool
end

const EBS = EnergyBackScatter

EnergyBackScatter(FT::DataType = Float64; 
                  C₄ = FT(0.06), 
                  C₂ = FT(-0.6), 
                  Cᴰ = FT(3e-2),
                  Uᵇ = FT(0.1),
                  τ  = FT(1 / 45days),
                  implicit_dissipation = true) = EnergyBackScatter(C₄, C₂, Cᴰ, Uᵇ, τ, implicit_dissipation) 

DiffusivityFields(grid, tracer_names, bcs, ::EBS) = 
                    (; ν₂ = CenterField(grid),
                       ν₄ = CenterField(grid), 
                       e  = Field((Center, Center, Nothing), grid),
                       G⁻ = Field((Center, Center, Nothing), grid),
    previous_compute_time = Ref(zero(grid)))

# For the moment just one!
@inline shape_function(i, j, k, grid) = one(grid)

@inline biharmonic_diffusion(i, j, k, grid, ::Val{true}, args...) = zero(grid)

@inline function biharmonic_diffusion(i, j, k, grid, implicit_dissipation, u, v, closure) 
    δ₁ = Dₛ(i, j, k, grid, u, v)    
    δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, u, v)    

    D_abs = sqrt(δ₁^2 + δ₂^2)
    
    C₄ = closure.C₄
    τ  = closure.τ

    return (C₄ * D_abs + τ) 
end

@kernel function _calculate_ebs_viscosities!(ν₂, ν₄, e, grid, closure, velocities, coriolis)
    i, j, k = @index(Global, NTuple)
    u, v, _ = velocities
    Δ = min(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

    ϕ = shape_function(i, j, k, grid)
    C₂ = closure.C₂

    @inbounds ν₄[i, j, k] = Δ^4 * biharmonic_diffusion(i, j, k, grid, Val(closure.implicit_dissipation), u, v, closure)

    Lᵦ = sqrt(u⁺) * sqrt(1 / abs(∂yᶠᶜᶜ(i, j, k, grid, fᶠᶠᵃ, coriolis)))
        
    Lᵐⁱˣ = min(Δ, Lᵦ)

    u⁺ = @inbounds sqrt(2 * max(0, e[i, j, 1]))

    ν₂_unbounded = C₂ * u⁺ * Lᵐⁱˣ * ϕ^2 

    ∂ˣu = ∂xᶜᶜᶜ(i, j, k, grid, u)
    ∂ʸv = ∂yᶜᶜᶜ(i, j, k, grid, v)
    
    ∂ˣv = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)
    ∂ʸu = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᶠᶠᶜ, u)

    bound = @inbounds 2 * e[i, j, 1] / sqrt(∂ˣu^2 + 0.5 * (∂ˣv + ∂ʸu)^2 + ∂ʸv^2)
    @inbounds ν₂[i, j, k] = max(ν₂_unbounded, - bound)
end

function compute_diffusivities!(diffusivity_fields, closure::EBS, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    barotropic_velocities = (; u = model.free_surface.state.U̅, v = model.free_surface.state.V̅) 
    advection = model.advection.b
    coriolis = model.coriolis
    Δt = model.clock.time - diffusivity_fields.previous_compute_time[]
    diffusivity_fields.previous_compute_time[] = model.clock.time

    launch!(arch, grid, :xy, 
            _advance_tke!, diffusivity_fields,
                           grid, model.clock, barotropic_velocities, velocities, advection, 
                           closure, Δt, Val(grid.Nz))

    launch!(arch, grid, parameters,
            _calculate_ebs_viscosities!,
            diffusivity_fields.ν₂, diffusivity_fields.ν₄, 
            diffusivity_fields.e, grid, closure, velocities, coriolis)

    return nothing
end

using Oceananigans.Operators

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2

@inline u_∂τˣ₂(i, j, k, grid, u, args...) = ∂yᶠᶠᶜ(i, j, k, grid, u) * viscous_flux_uy(i, j, k, grid, args...)
@inline v_∂τʸ₁(i, j, k, grid, v, args...) = ∂xᶠᶠᶜ(i, j, k, grid, v) * viscous_flux_vx(i, j, k, grid, args...)

@inline u_∂τˣ₁(i, j, k, grid, u, args...) = ∂xᶜᶜᶜ(i, j, k, grid, u) * viscous_flux_ux(i, j, k, grid, args...)
@inline v_∂τʸ₂(i, j, k, grid, v, args...) = ∂yᶜᶜᶜ(i, j, k, grid, v) * viscous_flux_vy(i, j, k, grid, args...)

@inline u_times_U_dot_∇u(i, j, k, grid, scheme, U) = @inbounds velocities.u[i, j, k] * U_dot_∇u(i, j, k, grid, scheme, U)
@inline v_times_U_dot_∇v(i, j, k, grid, scheme, U) = @inbounds velocities.v[i, j, k] * U_dot_∇v(i, j, k, grid, scheme, U)

@inline implicit_dissipation(i, j, k, grid, scheme, U, ::Val{false}) = zero(grid)

@inline function implicit_dissipation(i, j, k, grid, scheme, U, args...) 
    εˣ = ℑxᶜᵃᵃ(i, j, k, grid, u_times_U_dot_∇u, scheme, U)
    εʸ = ℑyᵃᶜᵃ(i, j, k, grid, v_times_U_dot_∇v, scheme, U)

    return εˣ + εʸ
end

@kernel function _advance_tke!(diffusivities, grid, clock, barotropic_velocities, velocities, advection, closure, Δt, ::Val{Nz}) where Nz
    i, j  = @index(Global, NTuple)
    u, v, _ = velocities

    εˣ = zero(grid)
    εʸ = zero(grid)

    εⁱ = zero(grid)

    @unroll for k in 1:Nz
        Δz = Δzᶜᶜᶜ(i, j, k, grid) / grid.Lz
        
        args = (closure, diffusivities, clock, velocities, nothing)

        εˣ += Δz * (u_∂τˣ₁(i, j, k, grid, u, args...) + ℑxyᶜᶜᵃ(i, j, k, grid, u_∂τˣ₂, u, args...))
        εʸ += Δz * (v_∂τʸ₂(i, j, k, grid, v, args...) + ℑxyᶜᶜᵃ(i, j, k, grid, v_∂τʸ₁, v, args...))

        args = (velocities, nothing)

        εⁱ += Δz * implicit_dissipation(i, j, k, grid, Val(closure.implicit_dissipation), advection, velocities)
    end

    FT = eltype(grid)
    α = convert(FT, 1.6)
    β = convert(FT, 0.6)

    Cᴰ = closure.Cᴰ

    e  = diffusivities.e
    G⁻ = diffusivities.G⁻

    s² = ℑxᶜᵃᵃ(i, j, 1, grid, ϕ², u) + ℑyᵃᶜᵃ(i, j, 1, grid, ϕ², v)

    implicit_bottom_friction = @inbounds Cᴰ * sqrt(s² + e[i, j, 1] + closure.Uᵇ^2) 

    Gⁿ = - div_Uc(i, j, 1, grid, advection, barotropic_velocities, e) - (εˣ + εʸ + εⁱ)

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
