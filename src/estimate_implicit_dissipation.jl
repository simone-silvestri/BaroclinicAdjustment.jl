function compute_χ_values(simulation)
    model = simulation.model
    advection = model.advection.b
    grid = model.grid
    arch = architecture(grid)
    b = model.tracers.b
    χ, bⁿ⁻¹, Uⁿ⁻¹ = simulation.model.auxiliary_fields

    launch!(arch, grid, :xyz, _compute_dissipation!, χ, grid, advection, 
            Uⁿ⁻¹, b, bⁿ⁻¹)

    return nothing
end

function update_previous_values(simulation)
    u, v, w = simulation.model.velocities
    b = simulation.model.tracers

    bⁿ⁻¹, Uⁿ⁻¹ = simulation.model.auxiliary_fields
    set!(bⁿ⁻¹,   b)
    set!(Uⁿ⁻¹.u, u)
    set!(Uⁿ⁻¹.v, v)
    set!(Uⁿ⁻¹.w, w)

    return nothing
end

@kernel function _compute_dissipation!(χᵁ, χⱽ, χᵂ, grid, advection, Uⁿ⁻¹, b, bⁿ⁻¹)
    i, j, k = @index(Global, NTuple)

    @inbounds χ.u[i, j, k] = compute_χᵁ(i, j, k, grid, advection, Uⁿ⁻¹.u, b, bⁿ⁻¹)
    @inbounds χ.v[i, j, k] = compute_χⱽ(i, j, k, grid, advection, Uⁿ⁻¹.v, b, bⁿ⁻¹)
    @inbounds χ.w[i, j, k] = compute_χᵂ(i, j, k, grid, advection, Uⁿ⁻¹.w, b, bⁿ⁻¹)
end

@inline b★(i, j, k, grid, bⁿ, bⁿ⁻¹) = @inbounds (bⁿ[i, j, k] + bⁿ⁻¹[i, j, k]) / 2
@inline b²(i, j, k, grid, b₁, b₂)   = @inbounds (b₁[i, j, k] * b₂[i, j, k])

@inline function compute_χᵁ(i, j, k, grid, advection, U, bⁿ, bⁿ⁻¹)
   
    δˣb★ = δxᶠᶜᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δˣb² = δxᶠᶜᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    𝒜x = _advective_tracer_flux_x(i, j, k, grid, advection, U, bⁿ⁻¹)
    𝒟x = @inbounds Axᶠᶜᶜ(i, j, k, grid) * U[i, j, k] * δˣb²

    return (𝒜x * 2 * δˣb★ - 𝒟x) / Vᶠᶜᶜ(i, j, k, grid)
end

@inline function compute_χⱽ(i, j, k, grid, advection, V, bⁿ, bⁿ⁻¹)
   
    δʸb★ = δyᶜᶠᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δʸb² = δyᶜᶠᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    𝒜y = _advective_tracer_flux_y(i, j, k, grid, advection, V, bⁿ⁻¹)
    𝒟y = @inbounds Ayᶜᶠᶜ(i, j, k, grid) * V[i, j, k] * δʸb²

    return (𝒜y * 2 * δʸb★ - 𝒟y) / Vᶜᶠᶜ(i, j, k, grid)
end

@inline function compute_χᵂ(i, j, k, grid, advection, W, bⁿ, bⁿ⁻¹)
   
    δᶻb★ = δzᶜᶜᶠ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δᶻb² = δzᶜᶜᶠ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    𝒜z = _advective_tracer_flux_z(i, j, k, grid, advection, W, bⁿ⁻¹)
    𝒟z = @inbounds Azᶜᶜᶠ(i, j, k, grid) * W[i, j, k] * δᶻb²

    return (𝒜z * 2 * δᶻb★ - 𝒟z) / Vᶜᶜᶠ(i, j, k, grid)
end