
"""
    struct QGLeith{FT, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}

The `QGLeith` struct represents a parameterization for the Leith scalar diffusivity in a quasi-geostrophic model.

Fields
=======
- `C`: The coefficient of the Leith scalar diffusivity.
- `min_N²`: The minimum value of the buoyancy frequency squared.
- `isopycnal_tensor`: The isopycnal tensor used in the parameterization.
- `slope_limiter`: The slope limiter used in the parameterization.

Reference
=========
Bachman, S., Fox-Kemper, B., & Pearson, B. (2017). A scale-aware subgrid model for
quasi-geostrophic turbulence. Journal of Geophysical Research: Oceans. 
doi: 10.1002/2016JC012265
"""
struct QGLeith{FT, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C :: FT
    min_N² :: FT
    isopycnal_tensor :: M
    slope_limiter :: S
end
QGLeith(FT::DataType = Float64; C=FT(1.0), min_N² = FT(1e-20), isopycnal_model=SmallSlopeIsopycnalTensor(), slope_limiter=FluxTapering(1e-2)) =
    QGLeith(C, min_N², isopycnal_model, slope_limiter) 

build_diffusivity_fields(grid, clock, tracer_names, bcs, ::QGLeith) = 
                (; νₑ = CenterField(grid),
                   Ld = Field{Center, Center, Nothing}(grid))

@inline function abs²_∇h_ζ(i, j, k, grid, coriolis, fields)

    ∂xζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)

    ∂xf = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, fᶠᶠᵃ, coriolis)
    ∂yf = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, fᶠᶠᵃ, coriolis)
    
    return ∂xζ + ∂xf, ∂yζ + ∂yf
end

@inline function abs²_∇h_δ(i, j, k, grid, fields)

    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    
    return (∂xδ^2 + ∂yδ^2)
end

@inline ∂yb_times_f2_div_N2(i, j, k, grid, clo, coriolis, buoyancy, tracers) = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) / 
                                                                               max(clo.min_N², ∂z_b(i, j, k, grid, buoyancy, tracers)) *
                                                                               ℑyzᵃᶜᶠ(i, j, k, grid, ∂y_b, buoyancy, tracers)

@inline ∂xb_times_f2_div_N2(i, j, k, grid, clo, coriolis, buoyancy, tracers) = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) / 
                                                                               max(clo.min_N², ∂z_b(i, j, k, grid, buoyancy, tracers))  *
                                                                               ℑxzᶜᵃᶠ(i, j, k, grid, ∂x_b, buoyancy, tracers)

@inline function abs²_∇h_q(i, j, k, grid, closure, coriolis, buoyancy, tracers)

    ∂zqx = ∂zᶜᶜᶜ(i, j, k, grid, ∂xb_times_f2_div_N2, closure, coriolis, buoyancy, tracers)
    ∂zqy = ∂zᶜᶜᶜ(i, j, k, grid, ∂yb_times_f2_div_N2, closure, coriolis, buoyancy, tracers)

    return ∂zqx, ∂zqy
end

@kernel function _calculate_qgleith_viscosity!(ν, Ld, grid, closure, velocities, tracers, buoyancy, coriolis)
    i, j, k = @index(Global, NTuple)

    ∂ζx, ∂ζy =  abs²_∇h_ζ(i, j, k, grid, coriolis, velocities)
    ∂qx, ∂qy =  abs²_∇h_q(i, j, k, grid, closure, coriolis, buoyancy, tracers)

    ∂δ² =  abs²_∇h_δ(i, j, k, grid, velocities)

    fᶜᶜᶜ = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis)

    ∂ζ² = ∂ζx^2 + ∂ζy^2
    ∂q² = (∂qx + ∂ζx)^2 + (∂qy + ∂ζy)^2

    A  = Δ²ᶜᶜᶜ(i, j, k, grid)
    Δs = A^0.5

    Bu  = Ld[i, j, 1]^2 / A
    Ro  = 1 / (fᶜᶜᶜ * Δs)
    
    ∂Q² = min(∂q², ∂ζ² * (1 + 1 / Bu)^2)
    ∂Q² = min(∂Q², ∂ζ² * (1 + 1 / Ro^2)^2)

    C = closure.C

    @inbounds ν[i, j, k] = (C * Δs / π)^(3) * sqrt(∂Q² + ∂δ²) 
end

@inline _deformation_radius(i, j, k, grid, C, buoyancy, coriolis) = sqrt(max(0, ∂z_b(i, j, k, grid, buoyancy, C))) / π /
                                                                         abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis))

@kernel function _calculate_deformation_radius!(Ld, grid, tracers, buoyancy, coriolis)
    i, j = @index(Global, NTuple)

    @inbounds begin
        Ld[i, j, 1] = 0
        @unroll for k in 1:grid.Nz
            Ld[i, j, 1] += Δzᶜᶜᶠ(i, j, k, grid) * _deformation_radius(i, j, k, grid, tracers, buoyancy, coriolis)
        end
    end
end

function compute_diffusivities!(diffusivity_fields, closure::QGLeith, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    coriolis = model.coriolis

    launch!(arch, grid, :xy, 
            _calculate_deformation_radius!, diffusivity_fields.Ld, grid, tracers, buoyancy, coriolis)

    launch!(arch, grid, parameters,
            _calculate_qgleith_viscosity!,
            diffusivity_fields.νₑ, diffusivity_fields.Ld, grid, closure, velocities, tracers, buoyancy, coriolis)

    return nothing
end

@inline viscosity(::QGLeith, K) = K.νₑ
@inline diffusivity(::QGLeith, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

@inline diffusive_flux_x(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)

#=
@inline function diffusive_flux_x(i, j, k, grid, closure::QGLeith, diffusivities, 
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ    = diffusivities.νₑ
    νₑⁱʲᵏ = ℑxᶠᵃᵃ(i, j, k, grid, νₑ)
    ∂x_c  = ∂xᶠᶜᶜ(i, j, k, grid, c)

    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)

    return - νₑⁱʲᵏ * ϵ * ∂x_c
end

@inline function diffusive_flux_y(i, j, k, grid, closure::QGLeith, diffusivities,
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ    = diffusivities.νₑ
    νₑⁱʲᵏ = ℑyᶠᵃᵃ(i, j, k, grid, νₑ)
    ∂y_c  = ∂yᶜᶠᶜ(i, j, k, grid, c)
    
    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)

    return - νₑⁱʲᵏ * ϵ * ∂y_c
end

@inline function diffusive_flux_z(i, j, k, grid, closure::QGLeith, diffusivities, 
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ = diffusivities.νₑ

    νₑⁱʲᵏ = ℑzᵃᵃᶠ(i, j, k, grid, νₑ)

    ∂x_c = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, c)
    ∂y_c = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, c)
    ∂z_c = ∂zᶜᶜᶠ(i, j, k, grid, c)

    R₃₁ = isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_model)
    R₃₂ = isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_model)
    R₃₃ = isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_model)

    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)

    return - νₑⁱʲᵏ * ϵ * (
        2 * R₃₁ * ∂x_c +
        2 * R₃₂ * ∂y_c + 
            R₃₃ * ∂z_c)
end
=#