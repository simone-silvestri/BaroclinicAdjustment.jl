using Oceananigans.Coriolis: hack_sind
using Oceananigans.Grids
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ZStar
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

function baroclinic_adjustment_simulation(testcase::TestCase, resolution, trailing = ""; kwargs...) 
    name = testcase.n * trailing
    momentum_advection = testcase.a
    horizontal_closure = testcase.h
    timestepper = testcase.t

    @info "Running case $name with" momentum_advection horizontal_closure timestepper

    return baroclinic_adjustment_simulation(resolution, name; momentum_advection, timestepper, horizontal_closure, kwargs...)
end

@inline u_drag_bc(i, j, grid, clock, fields, C) = - C * fields.u[i, j, 1]
@inline v_drag_bc(i, j, grid, clock, fields, C) = - C * fields.v[i, j, 1]

@inline u_immersed_drag_bc(i, j, k, grid, clock, fields, C) = - C * fields.u[i, j, k]
@inline v_immersed_drag_bc(i, j, k, grid, clock, fields, C) = - C * fields.v[i, j, k]

function u_bottom_drag()
    u_drag = FluxBoundaryCondition(u_drag_bc, discrete_form=true, parameters=5e-3)
    u_immersed_drag = FluxBoundaryCondition(u_immersed_drag_bc, discrete_form=true, parameters=5e-3)
    return FieldBoundaryConditions(bottom = u_drag, immersed = u_immersed_drag)
end

function v_bottom_drag()
    v_drag = FluxBoundaryCondition(v_drag_bc, discrete_form=true, parameters=5e-3)
    v_immersed_drag = FluxBoundaryCondition(v_immersed_drag_bc, discrete_form=true, parameters=5e-3)
    return FieldBoundaryConditions(bottom = v_drag, immersed = v_immersed_drag)
end

function gaussian_ridge(Lz, φ₀)
    function immersed_bottom(x, y)
        return - Lz + exp(- ((y - φ₀)^2) / 2.5^2) * Lz / 2
    end

    return immersed_bottom
end

function add_immersed_boundary(grid, Lz, φ₀)
    immersed_boundary = GridFittedBottom(gaussian_ridge(Lz, φ₀))
    return ImmersedBoundaryGrid(grid, immersed_boundary; active_cells_map = true)
end

function baroclinic_adjustment_simulation(resolution, filename, FT::DataType = Float64; 
                                          arch = GPU(), 
                                          horizontal_closure = nothing,
                                          momentum_advection = VectorInvariant(), 
                                          immersed = false,
                                          tracer_advection = WENO(FT; order = 7),
                                          auxiliary_fields = NamedTuple(),
                                          buoyancy_forcing_timescale = 50days,
                                          timestepper = :QuasiAdamsBashforth2,
                                          background_νz = 1e-4,
                                          φ₀ = - 50,
                                          stop_time = 1000days)
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 50
    Δt = if timestepper == :QuasiAdamsBashforth2
        5minutes
    else
        15minutes
    end

    grid = LatitudeLongitudeGrid(arch, FT;
                                 topology = (Periodic, Bounded, Bounded),
                                 size = (Ny, Ny, Nz), 
                                 longitude = (-10, 10),
                                 latitude = (φ₀-10, φ₀+10),
                                 z = MutableVerticalDiscretization((-Lz, 0)),
                                 halo = (6, 6, 6))
    
    if immersed
        grid = add_immersed_boundary(grid, Lz, φ₀)
    end

    vertical_closure = VerticalScalarDiffusivity(FT; κ=1e-5, ν=background_νz)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    N² = 4e-6  # [s⁻²] buoyancy frequency / stratification
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

    free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=Δt)
    @info "Building a model..."

    if !isnothing(buoyancy_forcing_timescale)
        set!(param.U★, (y, z) -> uᵢ(1, y, z, param))
        set!(param.B★, (y, z) -> bᵢ(1, y, z, param))
        bf = Forcing(  buoyancy_forcing; discrete_form=true, parameters=param)
        uf = Forcing(u_velocity_forcing; discrete_form=true, parameters=param)
        vf = Forcing(v_velocity_forcing; discrete_form=true, parameters=param)
        forcing = (; b = bf, u = uf, v = vf)
        boundary_conditions = NamedTuple()
    else 
        forcing = NamedTuple()
        boundary_conditions=(u=u_bottom_drag(), v=v_bottom_drag())
    end

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          buoyancy = BuoyancyTracer(),
                                          closure = closures,
                                          tracers = :b,
                                          momentum_advection,
                                          tracer_advection,
                                          auxiliary_fields,
                                          timestepper,
                                          forcing,
                                          boundary_conditions,
                                          vertical_coordinate = ZStar(),
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
    
    if !isnothing(buoyancy_forcing_timescale)
        simulation.callbacks[:update_mean_values] = Callback(update_mean_values, IterationInterval(10))
    end
    
    reduced_outputs!(simulation, filename)

    return simulation
end
