
struct Smagorinsky{FT} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    Cₛ :: FT
end

Smagorisky(FT::DataType = Float64; C = FT(0.15)) = Smagorinsky(C)

DiffusivityFields(grid, tracer_names, bcs, ::Smagorinsky) = 
                (; νₑ = CenterField(grid))


@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)

@kenrel function _calculate_smagorinsky_viscosity!(νₑ, grid, closure, velocities)
    i, j, k = @index(Global, NTuple)
    u, v, w = velocities

    δ₁ = Dₛ(i, j, k, grid, u, v)    
    δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, u, v)    
    A  = Δ²ᶜᶜᶜ(i, j, k, grid)
 
    @inbounds νₑ[i, j, k] = A * closure.Cₛ * sqrt(δ₁^2 + δ₂^2)
end

function compute_diffusivities!(diffusivity_fields, closure::QGLeith, model; parameters)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities

    launch!(arch, grid, parameters,
            _calculate_smagorinsky_viscosity!,
            diffusivity_fields.νₑ, grid, closure, velocities)

    return nothing
end

@inline viscosity(::Smagorinsky, K) = K.νₑ
@inline diffusivity(::Smagorinsky, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

@inline diffusive_flux_x(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)

