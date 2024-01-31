using Oceananigans.Coriolis: hack_sind
using Statistics

@inline function buoyancy_forcing(i, j, k, grid, clock, fields, p)
    @inbounds B  = p.B[1, j, k]
    @inbounds B★ = p.B★[1, j, k]
    return 1 / p.τ * (B★ - B)
end

@inline function u_velocity_forcing(i, j, k, grid, clock, fields, p)
    @inbounds U  = p.U[1, j, k]
    @inbounds U★ = p.U★[1, j, k]
    return 1 / p.τ * (U★ - U)
end

@inline function v_velocity_forcing(i, j, k, grid, clock, fields, p)
    @inbounds V = p.V[1, j, k]
    return - V / p.τ 
end

@inline transform(φ, p) = (p.φinit + φ) / p.Lφ * 2π - π/2

@inline function bᵢ(λ, φ, z, p) 
    x = transform(φ, p)
    b = ifelse(x < 0, 0, ifelse(x > π, 1, 1 - (π - x - sin(π - x) * cos(π - x)) / π))
    return p.N² * z + p.Δb * b
end

@inline function uᵢ(λ, φ, z, p) 
    f   = 2 * p.coriolis.rotation_rate * hack_sind(φ)
    x   = transform(φ, p)
    ∂b∂x = p.Δb * ifelse(x < 0, 0, ifelse(x > π, 0, (sin(x)^2 - cos(x)^2 + 1) / π))
    ∂x∂φ = 1 / p.Lφ * 2π
    ∂φ∂y = 1 / p.R * 180 / π
    return - 1 / f * (p.Lz + z) * ∂b∂x * ∂x∂φ * ∂φ∂y
end

function baroclinic_adjustment_latlong(testcase::TestCase, resolution, trailing; kwargs...) 
    name = testcase.n * trailing
    momentum_advection = testcase.a
    horizontal_closure = testcase.h

    return baroclinic_adjustment_latlong(resolution, name; momentum_advection, horizontal_closure, kwargs...)
end

function baroclinic_adjustment_latlong(resolution, filename, FT::DataType = Float64; 
                                       arch = GPU(), 
                                       horizontal_closure = nothing,
                                       momentum_advection = VectorInvariant(), 
                                       tracer_advection = WENO(FT; order = 7),
                                       buoyancy_forcing_timescale = 50days,
                                       background_νz = 1e-4,
                                       φ₀ = - 50,
                                       stop_time = 1000days)
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 50
    Δt = 2.5minutes

    grid = LatitudeLongitudeGrid(arch, FT;
                                topology = (Periodic, Bounded, Bounded),
                                size = (Ny, Ny, Nz), 
                                longitude = (-10, 10),
                                latitude = (φ₀-10, φ₀+10),
                                z = (-Lz, 0),
                                halo = (6, 6, 6))
    
    vertical_closure = VerticalScalarDiffusivity(FT; κ = 1e-5, ν = background_νz)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification
    Δb = 0.005 # [m/s²] buoyancy difference

    coriolis = HydrostaticSphericalCoriolis(FT)

    # Parameters
    param = (; N², Δb, Lz, Lφ = grid.Ly, 
               φinit = - (φ₀-10),
               τ  = buoyancy_forcing_timescale, 
               B  = Field{Nothing, Center, Center}(grid), 
               U  = Field{Nothing, Center, Center}(grid), 
               V  = Field{Nothing, Face,   Center}(grid),
               U★ = Field{Nothing, Center, Center}(grid),
               B★ = Field{Nothing, Center, Center}(grid),
               coriolis,
               R = grid.radius)

    free_surface = SplitExplicitFreeSurface(FT; grid, cfl = 0.75)
    @info "Building a model..."

    if !isnothing(buoyancy_forcing_timescale)
        set!(param.U★, (y, z) -> uᵢ(1, y, z, param))
        set!(param.B★, (y, z) -> bᵢ(1, y, z, param))
        bf = Forcing(  buoyancy_forcing; discrete_form=true, parameters=param)
        uf = Forcing(u_velocity_forcing; discrete_form=true, parameters=param)
        vf = Forcing(v_velocity_forcing; discrete_form=true, parameters=param)
        forcing = (; b = bf, u = uf, v = vf)
    else 
        forcing = NamedTuple()
    end

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          buoyancy = BuoyancyTracer(),
                                          closure = closures,
                                          tracers = :b,
                                          momentum_advection,
                                          tracer_advection,
                                          forcing,
                                          free_surface)

    ϵb = 1e-2 * Δb # noise amplitude
    Random.seed!(1234)

    bᵢᵣ(x, y, z) = bᵢ(x, y, z, param) + ϵb * randn()
    uᵢᵣ(x, y, z) = uᵢ(x, y, z, param)

    set!(model, b=bᵢᵣ, u=uᵢᵣ)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

    # add timestep wizard callback
    wizard = TimeStepWizard(cfl=0.1; max_change=1.1, max_Δt = 20minutes, min_Δt = 15)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

    function update_mean_values(sim)
        sim.model.forcing.b.parameters.B .= mean(sim.model.tracers.b,    dims = 1)
        sim.model.forcing.u.parameters.U .= mean(sim.model.velocities.u, dims = 1)
        sim.model.forcing.v.parameters.V .= mean(sim.model.velocities.v, dims = 1)
        return nothing
    end

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

    simulation.callbacks[:print_progress]     = Callback(print_progress,     IterationInterval(20))
    simulation.callbacks[:update_mean_values] = Callback(update_mean_values, IterationInterval(10))

    reduced_outputs!(simulation, filename)

    return simulation
end
