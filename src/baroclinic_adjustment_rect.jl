function baroclinic_adjustment_rectilinear(resolution, filename; arch = GPU(), 
                                           horizontal_closure = nothing,
                                           momentum_advection = VectorInvariant(), 
                                           xdirection = true)

    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(2000kilometers ÷ (resolution * 100kilometers))
    Nz = 50
    stop_time = 200days
    Δt = 5.0minutes

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

    coriolis = BetaPlane(latitude = -50)

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = one(grid),
                                                               convective_νz = zero(grid),
                                                               background_κz = 1e-5,
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

    Δb = 0.006
    Δy = 100kilometers
    ϵb = 1e-2 * Δb # noise amplitude

    Random.seed!(1234)
    bᵢ(x, y, z) = N² * z + Δb * ramp(y, Δy) + ϵb * randn()
        
    set!(model, b=bᵢ)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=10minutes)
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
