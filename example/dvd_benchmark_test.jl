using BaroclinicAdjustment
using BaroclinicAdjustment: baroclinic_adjustment_latlong
using Oceananigans
using Oceananigans.Units

using Oceananigans.Utils: launch!
using Oceananigans.Grids: architecture
using KernelAbstractions: @kernel, @index

using Oceananigans.Operators
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y, _advective_tracer_flux_z


# architecture and resolution
arch = CPU()
resolution = 1/2 # resolution in degrees
filename = "benchmark_dvd_sixteen"

# Hypothesis: we can change the order to reduce the variance dissipation
tracer_advection = WENO(order = 5)

# Momentum advection: the best one we have yet
momentum_advection = VectorInvariant(vorticity_scheme = WENO(order = 9), 
                                      vertical_scheme = WENO(order = 5))

# Horizontal dissipation: we don't need it when using WENO, otherwise uncomment below
horizontal_closure = nothing 

# @inline my_ν(i, j, k, grid, clock, fields) = fields.u[i, j, k]
# horizontal_closure = HorizontalScalarBiharmonicDiffusivity(ν = my_v, κ = 1e2, discrete_form = true)

# Simulation parameters
buoyancy_forcing_timescale = 50days
stop_time                  = 1000days

simulation = baroclinic_adjustment_latlong(resolution, filename; arch)

# This is the grid
grid = simulation.model.grid

# Define the 3 additional dissipation fields
bⁿ⁻¹ = CenterField(grid)
χᵁ = XFaceField(grid)
χⱽ = YFaceField(grid)
χᵂ = ZFaceField(grid)

# Construct the simulation
simulation = baroclinic_adjustment_latlong(resolution, filename; 
                                           arch,
                                           tracer_advection, 
                                           momentum_advection,
                                           horizontal_closure,
                                           auxiliary_fields = (; bⁿ⁻¹, χᵁ, χⱽ, χᵂ),
                                           buoyancy_forcing_timescale,
                                           stop_time)

@inline function store_previous_b(simulation) 

    bⁿ   = simulation.model.tracers.b
    bⁿ⁻¹ = simulation.model.auxiliary_fields.bⁿ⁻¹

    parent(bⁿ⁻¹) .= parent(bⁿ)

    return nothing
end
                
@inline b★(i, j, k, grid, b, bⁿ⁻¹) = @inbounds (b[i, j, k] + bⁿ⁻¹[i, j, k]) / 2

@kernel function _compute_χ(χᵁ, χⱽ, χᵂ, U, V, W, b, bⁿ⁻¹, grid, advection)
    i, j, k = @index(Global, NTuple)

    δˣb★ = δxᶠᶜᶜ(i, j, k, grid, b★, b, bⁿ⁻¹)
    δʸb★ = δyᶜᶠᶜ(i, j, k, grid, b★, b, bⁿ⁻¹)
    δᶻb★ = δzᶜᶜᶠ(i, j, k, grid, b★, b, bⁿ⁻¹)

    𝒜x = _advective_tracer_flux_x(i, j, k, grid, advection, U, b) 
    𝒜y = _advective_tracer_flux_y(i, j, k, grid, advection, V, b) 
    𝒜z = _advective_tracer_flux_z(i, j, k, grid, advection, W, b) 

    @inbounds begin
        χᵁ[i, j, k] = 𝒜x * 2 * δˣb★ / Vᶠᶜᶜ(i, j, k, grid)
        χⱽ[i, j, k] = 𝒜y * 2 * δʸb★ / Vᶜᶠᶜ(i, j, k, grid)
        χᵂ[i, j, k] = 𝒜z * 2 * δᶻb★ / Vᶜᶜᶠ(i, j, k, grid)
    end
end

@inline function compute_χ(simulation)
    model = simulation.model
    grid  = model.grid
    
    (; bⁿ⁻¹, χᵁ, χⱽ, χᵂ) = model.auxiliary_fields
    (; u, v, w) = model.velocities
    b = model.tracers.b

    launch!(architecture(grid), grid, :xyz, _compute_χ, χᵁ, χⱽ, χᵂ, u, v, w, b, bⁿ⁻¹, grid, model.advection.b)

    return nothing
end

simulation.callbacks[:compute_χ] = Callback(compute_χ, IterationInterval(1))
simulation.callbacks[:store_previous_b] = Callback(store_previous_b, IterationInterval(1))

simulation.output_writers[:dissipation] = JLD2OutputWriter(simulation.model, (; χᵁ, χⱽ, χᵂ);
                                                           filename = filename * "_chi", 
                                                           schedule = TimeInterval(5days),
                                                           overwrite_existing = true)

# Run it!!!
run!(simulation)