module BaroclinicAdjustment

using Random
using Oceananigans: getnamewrapper
using Oceananigans.Units
using GLMakie
using Oceananigans

include("horizontal_visc.jl")
include("outputs.jl")

function baroclinic_adjustment(resolution; horizontal_closure = nothing, momentum_advection = VectorInvariant())

    # Architecture
    arch = GPU()
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Int(10 / resolution)
    Nz = 100
    stop_time = 200days
    Δt = 10minutes

    grid = LatitudeLongitudeGrid(arch;
                                topology = (Periodic, Bounded, Bounded),
                                size = (Ny, Ny, Nz), 
                                longitude = (-5,   5),
                                latitude = (-60, -50),
                                z = (-Lz, 0),
                                halo = (3, 3, 3))

    coriolis = HydrostaticSphericalCoriolis()

    Δy = 1000kilometers / Ny
    vertical_closure = VerticalScalarDiffusivity(ν=1e-4, κ=1e-5)

    closures = (vertical_closure, horizontal_closure)

    @info "Building a model..."

    model = HydrostaticFreeSurfaceModel(grid = grid,
                                        coriolis = coriolis,
                                        buoyancy = BuoyancyTracer(),
                                        closure = closures,
                                        tracers = :b,
                                        momentum_advection = VectorInvariant(),
                                        tracer_advection = WENO(),
                                        free_surface = ImplicitFreeSurface())

    @info "Built $model."

    """
    Linear ramp from 0 to 1 between -Δy/2 and +Δy/2.

    For example:

    y < y₀           => ramp = 0
    y₀ < y < y₀ + Δy => ramp = y / Δy
    y > y₀ + Δy      => ramp = 1
    """
    function ramp(λ, y, Δ)
        gradient == "x" && return min(max(0, λ / Δ + 1/2), 1)
        gradient == "y" && return min(max(0, (y - 45) / Δ + 1/2), 1)
    end

    # Parameters
    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification
    M² = 8e-8 # [s⁻²] horizontal buoyancy gradient

    Δy = 1 # degree
    Δb = 100kilometers * Δy * M²

    bᵢ(λ, y, z) = N² * z + Δb * ramp(λ, y, Δy)
    
    set!(model, b=bᵢ)

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

    filename = "experiment_$(viscosity_name(horizontal_closure))_$(advection_name(momentum_advection))"
    
    standard_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end

viscosity_name(clo::ScalarBiharmonicDiffusivity) = typeof(clo.ν)
viscosity_name(clo::ScalarDiffusivity)           = typeof(clo.ν)
viscosity_name(clo) = typeof(clo)
advection_name(adv) = getnamewrapper(adv.vorticity_scheme)

function run_quarter_degree_simulations()

    vi1 = VectorInvariant()
    vi2 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO())
    vi3 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO(), vorticity_stencil  = DefaultStencil())
    vi4 = VectorInvariant(vorticity_scheme = WENO(), divergence_scheme = WENO(), vertical_scheme = WENO(), divergence_stencil = VelocityStencil())
    vi5 = VectorInvariant(vorticity_scheme = WENO(order = 9), divergence_scheme = WENO(), vertical_scheme = WENO())

    hi1 = nothing
    hi2 = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)
    hi3 = leith_viscosity(HorizontalFormulation())
    hi3 = leith_laplacian_viscosity(HorizontalFormulation())
    hi5 = smagorinski_viscosity(HorizontalFormulation())

    advection_schemes   = [vi1, vi1, vi1, vi1, vi2, vi3, vi4, vi5]
    horizontal_closures = [hi2, hi3, hi4, hi5, hi1, hi1, hi1, hi1]


    for (momentum_advection, horizontal_closure) in zip(advection_schemes, horizontal_closures)
        baroclinic_adjustment(1/4; momentum_advection, horizontal_closure)
    end

    return nothing
end


end
