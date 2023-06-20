using Statistics: mean
using Oceananigans.Grids: node

function restoring_baroclinic_adjustment_latlong(resolution, filename, FT::DataType = Float64; arch = GPU(), 
                                                 horizontal_closure = nothing,
                                                 momentum_advection = VectorInvariant(), 
                                                 background_νz = 1e-4,
                                                 φ₀ = -50)
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 50
    stop_time = 400days
    Δt = 2.5minutes

    grid = LatitudeLongitudeGrid(arch, FT;
                                topology = (Periodic, Bounded, Bounded),
                                size = (Ny, Ny, Nz), 
                                longitude = (-10, 10),
                                latitude = (φ₀-10, φ₀+10),
                                z = (-Lz, 0),
                                halo = (6, 6, 6))

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(FT; convective_κz = 0.1,
                                                                   convective_νz = 0.0,
                                                                   background_κz = 1e-5,
                                                                   background_νz)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    gravity = FT(Oceananigans.BuoyancyModels.g_Earth)
    max_Δt  = 20minutes

    substeps = barotropic_substeps(max_Δt, grid, gravity)

    @inline ramp(λ, y, Δ) = min(max(0, (φ₀ + y) / Δ + 1/2), 1)

    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification
    Δy = 1.0 # degree
    Δb = 0.006

    # Parameters
    param = (; N², Δy, Δb, ramp)

    @inline bᵢ(x, y, z, p) = p.N² * z + p.Δb * p.ramp(x, y, p.Δy) 

    B = Field{Nothing, Center, Center}(grid)

    @inline function B_forcing(i, j, k, grid, clock, fields, p)
        λ, φ, z = node(i, j, k, grid, Center(), Center(), Center())
        @inbounds Bᴹ = fields.B[1, j, k]
        @inbounds bᴿ = p.bᵢ(λ, φ, z, p)
        return - p.λ * (Bᴹ - bᴿ)
    end

    Fb = Forcing(B_forcing; discrete_form = true, parameters = (; ramp, bᵢ, N², Δy, Δb, λ = 1 / 50days))

    @inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

    @inline speedᶠᶜᶜ(i, j, k, grid, fields) = sqrt(fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))
    @inline speedᶜᶠᶜ(i, j, k, grid, fields) = sqrt(ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u) + fields.v[i, j, k]^2)
 
    @inline u_drag(i, j, grid, clock, fields, τ) = - τ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
    @inline v_drag(i, j, grid, clock, fields, τ) = - τ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

    u_drag_bcs = FluxBoundaryCondition(u_drag; parameters = 0.001)    
    v_drag_bcs = FluxBoundaryCondition(v_drag; parameters = 0.001)

    u_bcs = FieldBoundaryConditions(bottom = u_drag_bcs)
    v_bcs = FieldBoundaryConditions(bottom = v_drag_bcs)

    free_surface = SplitExplicitFreeSurface(FT; substeps)
    @info "running with substeps $substeps"
    @info "Building a model..."

    model = HydrostaticFreeSurfaceModel(; grid,
                                        coriolis = HydrostaticSphericalCoriolis(FT),
                                        buoyancy = BuoyancyTracer(),
                                        closure = closures,
                                        tracers = :b,
                                        momentum_advection,
                                        auxiliary_fields = (; B),
                                        forcing = (; b = Fb),
                                        boundary_conditions = (; u = u_bcs, v = v_bcs),
                                        tracer_advection = WENO(FT),
                                        free_surface)

    @info "Built $model."

    ϵb = 1e-2 * Δb # noise amplitude
    Random.seed!(1234)
    bᵢᵣ(x, y, z) = bᵢ(x, y, z, param) + ϵb * randn()

    set!(model, b=bᵢᵣ)

    model.auxiliary_fields.B .= mean(model.tracers.b, dims = 1)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    @show simulation.model.auxiliary_fields.B

    function set_average_B!(sim)
        sim.model.auxiliary_fields.B .= mean(sim.model.tracers.b, dims=1)
        return nothing
    end

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.1; max_change=1.1, max_Δt, min_Δt = 15)
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
    simulation.callbacks[:set_average_B!] = Callback(set_average_B!, IterationInterval(1))

    reduced_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end