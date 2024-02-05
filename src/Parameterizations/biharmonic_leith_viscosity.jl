"""
    struct BiharmonicLeith{FT} <: AbstractScalarBiharmonicDiffusivity{HorizontalFormulation, 3}

The `BiharmonicLeith` struct represents a biharmonic Leith viscosity parameterization.

Fields
======
- `C::FT`: The coefficient of the biharmonic Leith viscosity.

Reference
=========
Fox-Kemper, B., & Menemenlis, D. (2004). Can large eddy simulation techniques improve
mesoscale rich ocean models? In M. Hecht & H. Hasumi (Eds.), Ocean modeling in an eddying regime 
(pp. 319 - 337). doi: 10.1029/177GM19
"""
struct BiharmonicLeith{FT} <: AbstractScalarBiharmonicDiffusivity{HorizontalFormulation, 3}
    C :: FT
end

BiharmonicLeith(FT::DataType = Float64; C=FT(2.0)) = BiharmonicLeith(C) 

DiffusivityFields(grid, tracer_names, bcs, ::BiharmonicLeith) = 
                (; νₑ = CenterField(grid))

@kernel function _calculate_biharmonicleith_viscosity!(ν, grid, closure, fields)
    i, j, k = @index(Global, NTuple)

    ∂ζx = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂ζy = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    
    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)


    ∂ζ² = ∂ζx^2 + ∂ζy^2
    ∂δ² = ∂xδ^2 + ∂yδ^2

    A  = Δ²ᶜᶜᶜ(i, j, k, grid)
    Δs = A^0.5
    
    C = closure.C

    @inbounds ν[i, j, k] = (C * Δs / π)^(3) * sqrt(∂ζ² + ∂δ²) * Δs^2 / 8
end

function compute_diffusivities!(diffusivity_fields, closure::BiharmonicLeith, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities

    launch!(arch, grid, parameters,
            _calculate_biharmonicleith_viscosity!,
            diffusivity_fields.νₑ, grid, closure, velocities)

    return nothing
end

@inline viscosity(::BiharmonicLeith, K) = K.νₑ
@inline diffusivity(::BiharmonicLeith, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

@inline diffusive_flux_x(i, j, k, grid, closure::BiharmonicLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::BiharmonicLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::BiharmonicLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
