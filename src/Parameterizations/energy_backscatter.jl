using Oceananigans.Operators
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.TurbulenceClosures: biharmonic_mask_x, biharmonic_mask_y, ∂x_∇²h_cᶠᶜᶜ, ∂y_∇²h_cᶜᶠᶜ

struct EnergyBackScatter{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C₄ :: FT
    C₂ :: FT
    Cᴰ :: FT
end

const EBS = EnergyBackScatter

EnergyBackScatter(FT::DataType = Float64; C₄ = FT(0.2), C₂ = FT(-0.15), Cᴰ = FT(3e-3)) = EnergyBackScatter(C₄, C₂, Cᴰ) 

DiffusivityFields(grid, tracer_names, bcs, ::EBS) = 
                    (; ν₂ = CenterField(grid),
                       ν₄ = CenterField(grid), 
                       e  = CenterField(grid),
                       G⁻ = CenterField(grid),
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
    u⁺ = @inbounds sqrt(2 * max(0, e[i, j, k]))
    Lᵦ = sqrt(u⁺) * sqrt(1 / abs(∂yᶠᶜᶜ(i, j, k, grid, fᶠᶠᵃ, coriolis)))
    
    Lᵐⁱˣ = min(Δ, Lᵦ)

    @inbounds ν₂[i, j, k] = C₂ * u⁺ * Lᵐⁱˣ * ϕ^2 
    @inbounds ν₄[i, j, k] = C₄ * Δ^4 * D_abs
end

function compute_diffusivities!(diffusivity_fields, closure::EBS, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    advection = model.advection.b
    coriolis = model.coriolis
    Δt = model.clock.time - diffusivities.previous_compute_time[]
    diffusivities.previous_compute_time[] = model.clock.time

    launch!(arch, grid, parameters, 
            _advance_tke!, diffusivity_fields,
                           grid, velocities, advection, 
                           closure, Δt)

    launch!(arch, grid, parameters,
            _calculate_ebs_viscosities!,
            diffusivity_fields.ν₂, diffusivity_fields.ν₄, 
            diffusivity_fields.e, grid, closure, velocities, coriolis)

    return nothing
end

using Oceananigans.Operators

@kernel function _advance_tke!(diffusivities, grid, velocities, advection, closure, Δt)
    i, j, k = @index(Global, NTuple)

    FT = eltype(grid)
    α = convert(FT, 1.6)
    β = convert(FT, 0.6)

    Cᴰ = closure.Cᴰ

    e  = diffusivities.e
    G⁻ = diffusivities.G⁻
    Gⁿ = - div_Uc(i, j, k, grid, advection, velocities, e) 
         - ∇_dot_qᵉ(i, j, k, grid, closure, diffusivities, e)
         - Cᴰ * e[i, j, k]

    @inbounds e[i, j, k] += Δt * (α * G - β * G⁻[i, j, k])
    @inbounds G⁻[i, j, k] = Gⁿ
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
 
@inline function ∇_dot_qᵉ(i, j, k, grid, ::EnergyBackScatter, diffusivities, e)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Ax_qᶠᶜᶜ, _backscatter_diffusive_flux_x, diffusivities, e) +
                                    δyᵃᶜᵃ(i, j, k, grid, Ay_qᶜᶠᶜ, _backscatter_diffusive_flux_y, diffusivities, e))
end

@inline _backscatter_diffusive_flux_x(args...) = backscatter_diffusive_flux_x(args...)
@inline _backscatter_diffusive_flux_y(args...) = backscatter_diffusive_flux_y(args...)

@inline backscatter_diffusive_flux_x(i, j, k, grid, K, e) = @inbounds - K.ν₂[i, j, k] * ∂xᶠᶜᶜ(i, j, k, grid, e) + K.ν₄[i, j, k] * biharmonic_mask_x(i, j, k, grid, ∂x_∇²h_cᶠᶜᶜ, c)
@inline backscatter_diffusive_flux_y(i, j, k, grid, K, e) = @inbounds - K.ν₂[i, j, k] * ∂yᶜᶠᶜ(i, j, k, grid, e) + K.ν₄[i, j, k] * biharmonic_mask_y(i, j, k, grid, ∂y_∇²h_cᶜᶠᶜ, c)

@inline diffusive_flux_x(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, ::EBS, K, ::Val{id}, c, clk, fields, b) where id = zero(grid)
