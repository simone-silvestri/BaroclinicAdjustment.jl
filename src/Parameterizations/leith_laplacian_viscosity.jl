struct Leith{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C :: FT
end

Leith(FT::DataType = Float64; C=FT(1.0)) = Leith(C) 

DiffusivityFields(grid, tracer_names, bcs, ::Leith) = 
                (; νₑ = CenterField(grid))


@inline function abs²_∇h_δ(i, j, k, grid, fields)

    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    
    return (∂xδ^2 + ∂yδ^2)
end

@kernel function _calculate_leith_viscosity!(ν, grid, closure, fields)
    i, j, k = @index(Global, NTuple)

    ∂ζx = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    
    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)


    ∂ζ² = ∂ζx^2 + ∂ζy^2
    ∂δ² = ∂xδ^2 + ∂yδ^2

    A  = Δ²ᶜᶜᶜ(i, j, k, grid)
    Δs = A^0.5
    
    C = closure.C

    @inbounds ν[i, j, k] = (C * Δs / π)^(3) * sqrt(∂ζ² + ∂δ²) 
end

function compute_diffusivities!(diffusivity_fields, closure::Leith, model)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities

    launch!(arch, grid, :xyz,
            _calculate_leith_viscosity!,
            diffusivity_fields.νₑ, grid, closure, velocities)

    return nothing
end

@inline viscosity(::Leith, K) = K.νₑ
@inline diffusivity(::Leith, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

@inline diffusive_flux_x(i, j, k, grid, closure::Leith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::Leith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::Leith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
