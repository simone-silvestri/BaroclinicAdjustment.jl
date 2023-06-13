using Oceananigans.Grids: φnode

function baroclinic_double_gyre(resolution, filename, FT::DataType = Float64; arch = GPU(), 
                                horizontal_closure = nothing,
                                momentum_advection = VectorInvariant())

    # Domain
    stop_time = 5000days
    Δt = 2.5minutes

    λ_west = -30 # [°] longitude of west boundary
    λ_east = +30 # [°] longitude of east boundary
    φ_south = 15 # [°] latitude of south boundary
    φ_north = 75 # [°] latitude of north boundary

    φ₀ = 0.5(φ_north + φ_south) # [°] latitude of the center of the domain

    Lλ = λ_east - λ_west   # [°] longitude extent of the domain
    Lφ = φ_north - φ_south # [°] latitude extent of the domain

    Lz = 1.8kilometers # depth [m]

    Nλ = Base.Int(Lλ ÷ resolution)
    Nφ = Base.Int(Lφ ÷ resolution)
    Nz = 50

    σ = 1.2 # stretching factor
    hyperbolically_spaced_faces(k) = - Lz * (1 - tanh(σ * (k - 1) / Nz) / tanh(σ))

    grid = LatitudeLongitudeGrid(arch, FT;
                                 size = (Nλ, Nφ, Nz),
                                 longitude = (λ_west, λ_east),
                                 latitude = (φ_south, φ_north),
                                 z = hyperbolically_spaced_faces,
                                 topology = (Bounded, Bounded, Bounded),
                                 halo = (6, 6, 6))

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(FT; convective_κz = 0.1,
                                    convective_νz = 1e-2,
                                    background_κz = 1e-5,
                                    background_νz = 5e-4)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    gravity = FT(Oceananigans.BuoyancyModels.g_Earth)
    max_Δt  = 20minutes

    substeps     = barotropic_substeps(max_Δt, grid, gravity)    
    free_surface = SplitExplicitFreeSurface(FT; substeps)

    @info "running with substeps $substeps"
    @info "Building a model..."

    α  = 2e-4 # [K⁻¹] thermal expansion coefficient
    g  = 9.81 # [m s⁻²] gravitational constant
    ρ₀ = 1028 # [kg m⁻³] reference density

    Δzₛ = minimum_zspacing(grid) # vertical spacing at the surface [m] 

    parameters = (Lφ = Lφ,
                  Lz = Lz,
                  φ₀ = φ₀,                # latitude of the center of the domain [°]
                  φₛ = φ_south,           # latitude of the southern edge of the domain [°]
                  τ = 0.1 / ρ₀,           # surface kinematic wind stress [m² s⁻²]
                  μ = 0.001,              # bottom drag damping parameter [s⁻¹]
                  Δb = 30 * α * g,        # surface vertical buoyancy gradient [s⁻²]
                  timescale = 30days,     # relaxation time scale [s]  
                  vˢ = Δzₛ/30days)        # buoyancy pumping velocity [ms⁻¹]

    @inline u_stress(λ, φ, t, p) = p.τ * sin(2π * (φ - p.φ₀) / p.Lφ)

    # ### Surface buoyancy relaxation
    @inline function b_relax(i, j, grid, clock, fields, p) 
        bˢ = fields.b[i, j, grid.Nz]
        φ  = φnode(j, grid, Center())
        bᴿ = p.Δb * (φ - p.φₛ) / p.Lφ
        return p.vˢ * (bˢ - bᴿ)
    end

    # ### Bottom drag
    @inline u_drag(i, j, grid, clock, fields, p) = - p.μ * Δzᶠᶜᶜ(i, j, 1, grid) * fields.u[i, j, 1]
    @inline v_drag(i, j, grid, clock, fields, p) = - p.μ * Δzᶜᶠᶜ(i, j, 1, grid) * fields.v[i, j, 1]

    u_stress_bc = FluxBoundaryCondition(u_stress; parameters)
    b_relax_bc  = FluxBoundaryCondition(nothing) #b_relax;  parameters, discrete_form = true)
    u_drag_bc   = FluxBoundaryCondition(u_drag;   parameters, discrete_form = true)
    v_drag_bc   = FluxBoundaryCondition(v_drag;   parameters, discrete_form = true)

    u_bcs = FieldBoundaryConditions(bottom = u_drag_bc, top = u_stress_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_drag_bc)
    b_bcs = FieldBoundaryConditions(top = b_relax_bc)

    model = HydrostaticFreeSurfaceModel(; grid,
                                        coriolis = HydrostaticSphericalCoriolis(FT),
                                        buoyancy = BuoyancyTracer(),
                                        closure = closures,
                                        tracers = :b,
                                        boundary_conditions = (u = u_bcs, v = v_bcs, b = b_bcs),
                                        momentum_advection,
                                        tracer_advection = WENO(FT),
                                        free_surface)

    @info "Built $model."

    ## Initial conditions
    bᵢ(λ, φ, z) = parameters.Δb * ( 1 + z / grid.Lz )

    set!(model, b = bᵢ)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.2; max_change=1.1, max_Δt, min_Δt = 15)
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

    standard_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end