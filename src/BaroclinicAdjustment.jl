module BaroclinicAdjustment

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing, architecture
using Oceananigans.Operators
using Oceananigans.Utils: getnamewrapper, launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Advection: VelocityStencil, DefaultStencil
using KernelAbstractions: @kernel, @index
using JLD2

include("horizontal_visc.jl")
include("qg_leith_viscosity.jl")
include("outputs.jl")

@inline ∂z_uᴳ(i, j, k, grid, b, coriolis) = 1 / fᶠᶠᵃ(i, j, k, grid, coriolis) * ℑxyzᶠᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, b)

@kernel function _geostrophic_velocity(u, b, grid, coriolis)
    i, j = @index(Global, NTuple)

    @inbounds begin
        u[i, j, 1] = Δzᶠᶜᶠ(i, j, 1, grid) * ∂z_uᴳ(i, j, 1, grid, b, coriolis)
        for k in 2:grid.Nz
            u[i, j, k] = u[i, j, k-1] + Δzᶠᶜᶠ(i, j, k, grid) * ∂z_uᴳ(i, j, k, grid, b, coriolis)
        end
    end
end

function set_geostrophic_velocity!(u, b, coriolis)
    grid = u.grid
    launch!(architecture(grid), grid, :xy, _geostrophic_velocity, u, b, grid, coriolis)

    return nothing
end

function barotropic_substeps(Δt, grid, gravitational_acceleration)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δ  = 1 / sqrt(1 / Δx^2 + 1 / Δy^2)

    return  Base.Int(ceil(2 * Δt / (0.7 / wave_speed * Δ)))
end

function baroclinic_adjustment_rectilinear(resolution, filename; arch = GPU(), 
                                           horizontal_closure = nothing,
                                           momentum_advection = VectorInvariant(), 
                                           xdirection = true)

    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(2000kilometers ÷ (resolution * 100kilometers))
    Nz = 100
    stop_time = 200days
    Δt = 2.5minutes

    if xdirection 
        grid = RectilinearGrid(arch;
                               topology = (Periodic, Bounded, Bounded),
                               size = (Ny, Ny, Nz), 
                               x = (-1000kilometers, 1000kilometers),
                               y = (-1000kilometers, 1000kilometers),
                               z = (-Lz, 0),
                               halo = (6, 6, 6))
    else
        grid = RectilinearGrid(arch;
                               topology = (Flat, Bounded, Bounded),
                               size = (Ny, Nz), 
                               y = (-1000kilometers, 1000kilometers),
                               z = (-Lz, 0),
                               halo = (6, 6))
    end

    coriolis = BetaPlane(latitude = -45)

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = one(grid),
                                                               convective_νz = zero(grid),
                                                               background_κz = 1e-4,
                                                               background_νz = 1e-4)

    
    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    @info "Building a model..."

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          buoyancy = BuoyancyTracer(),
                                          closure = closures,
                                          tracers = :b,
                                          momentum_advection,
                                          tracer_advection = WENO(),
                                          free_surface = ImplicitFreeSurface())

    @info "Built $model."

    ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)

    # Parameters
    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification

    Δb = 0.06
    Δy = 200kilometers
        
    bᵢ(x, y, z) = N² * z + Δb * ramp(y, Δy)

    set!(model, b=bᵢ)
    # set_geostrophic_velocity!(model.velocities.u, model.tracers.b, model.coriolis)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

    # add progress callback
    wall_clock = [time_ns()]

    function print_progress(sim)
        @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
        100 * (sim.model.clock.time / sim.stop_time),
        sim.model.clock.iteration,
        prettytime(sim.model.clock.time),
        prettytime(1e-9 * (time_ns() - wall_clock[1])),
        maximum(abs, sim.model.velocities.u),
        maximum(abs, sim.model.velocities.v),
        maximum(abs, sim.model.velocities.w),
        prettytime(sim.Δt))

        wall_clock[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(20))

    reduced_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end
    
function baroclinic_adjustment_latlong(resolution, filename; arch = GPU(), 
                                                     horizontal_closure = nothing,
                                                     momentum_advection = VectorInvariant(), 
                                                     xdirection = true)
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 100
    stop_time = 200days
    Δt = 2.5minutes

    if xdirection 
        grid = LatitudeLongitudeGrid(arch;
                                    topology = (Periodic, Bounded, Bounded),
                                    size = (Ny, Ny, Nz), 
                                    longitude = (-10, 10),
                                    latitude = (-60, -40),
                                    z = (-Lz, 0),
                                    halo = (6, 6, 6))
    else
        grid = LatitudeLongitudeGrid(arch;
                                     topology = (Flat, Bounded, Bounded),
                                     size = (Ny, Nz), 
                                     latitude = (-60, -40),
                                     z = (-Lz, 0),
                                     halo = (6, 6))
    end

    coriolis = HydrostaticSphericalCoriolis()

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = one(grid),
                                                               convective_νz = zero(grid),
                                                               background_κz = 1e-4,
                                                               background_νz = 1e-4)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    @info "Building a model..."

    model = HydrostaticFreeSurfaceModel(; grid,
                                        coriolis,
                                        buoyancy = BuoyancyTracer(),
                                        closure = closures,
                                        tracers = :b,
                                        momentum_advection,
                                        tracer_advection = WENO(),
                                        free_surface = ImplicitFreeSurface())

    @info "Built $model."

    gradient = "y"

    function ramp(λ, y, Δ)
        gradient == "x" && return min(max(0, λ / Δ + 1/2), 1)
        gradient == "y" && return min(max(0, (50 + y) / Δ + 1/2), 1)
    end

    # Parameters
    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification

    Δy = 5 # degree
    Δb = 0.06

    bᵢ(λ, y, z) = N² * z + Δb * ramp(λ, y, Δy)

    set!(model, b=bᵢ)
    set_geostrophic_velocity!(model.velocities.u, model.tracers.b, model.coriolis)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

    # add progress callback
    wall_clock = [time_ns()]

    function print_progress(sim)
        @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
                100 * (sim.model.clock.time / sim.stop_time),
                sim.model.clock.iteration,
                prettytime(sim.model.clock.time),
                prettytime(1e-9 * (time_ns() - wall_clock[1])),
                maximum(abs, sim.model.velocities.u),
                maximum(abs, sim.model.velocities.v),
                maximum(abs, sim.model.velocities.w),
                prettytime(sim.Δt))

        wall_clock[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(20))

    reduced_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end

viscosity_name(clo::ScalarBiharmonicDiffusivity) = typeof(clo.ν)
viscosity_name(clo::ScalarDiffusivity)           = typeof(clo.ν)
viscosity_name(clo) = typeof(clo)
advection_name(adv) = getnamewrapper(adv.vorticity_scheme)

add_trailing_characters(name, trailing_character = "_weaker") = name * trailing_character

function run_eight_degree_simulations(trailing_character = "_weaker")

    vi1 = VectorInvariant()
    vi2 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO())
    vi3 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO(), vorticity_stencil  = DefaultStencil())
    vi4 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO(), divergence_stencil = VelocityStencil())
    vi5 = VectorInvariant(vorticity_scheme = WENO(order = 9), divergence_scheme = WENO(), vertical_scheme = WENO())
    vi6 = VectorInvariant(vorticity_scheme = WENO(order = 9), divergence_scheme = WENO(), vertical_scheme = WENO(), vorticity_stencil  = DefaultStencil())

    hi1 = nothing
    hi2 = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)
    hi3 = leith_viscosity(HorizontalFormulation())
    hi4 = leith_laplacian_viscosity(HorizontalFormulation())
    hi5 = smagorinski_viscosity(HorizontalFormulation())
    hi6 = QGLeithViscosity()

    advection_schemes   = [vi1, vi2, vi1, vi1, vi1, vi3, vi4, vi5, vi1, vi6]
    horizontal_closures = [hi2, hi1, hi3, hi4, hi5, hi1, hi1, hi1, hi6, hi1]

    names = ["bilap", "weno5vd", "leith", "lapleith", "smag", "weno5dd", "weno5vv", "weno9", "qgleith", "weno9dd"]
    names = add_trailing_characters.(names, Ref(trailing_character))
    
    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        baroclinic_adjustment_latlong(1/8, name; momentum_advection, horizontal_closure)
    end

    return nothing
end

function run_high_res_simulation(trailing_character = "_weaker")

    vi1 = VectorInvariant()

    hi4 = leith_laplacian_viscosity(HorizontalFormulation(), C_vort = 2.0, C_div = 2.0)

    advection_schemes   = [vi1]
    horizontal_closures = [hi4]

    names = ["highres"]
    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        baroclinic_adjustment_latlong(1/50, name; momentum_advection, horizontal_closure)
    end

    return nothing
end

run_2d_flat_simulation(trailing_character = "_weaker") = 
    baroclinic_adjustment_latlong(1/8, add_trailing_characters("2dsim", trailing_character);
                          horizontal_closure = leith_viscosity(HorizontalFormulation()),
                          xdirection = false)

function run_all(trailing_character = "_weaker")
    run_2d_flat_simulation(trailing_character)
    run_eight_degree_simulations(trailing_character)
    run_high_res_simulation(trailing_character)
end

include("Diagnostics/Diagnostics.jl")

using .Diagnostics

include("calculate_diagnostics.jl")

end
