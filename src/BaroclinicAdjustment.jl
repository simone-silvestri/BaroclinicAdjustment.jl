module BaroclinicAdjustment

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing, architecture
using Oceananigans.Operators
using Oceananigans.Utils: getnamewrapper, launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Advection: VelocityStencil, DefaultStencil, EnergyConservingScheme

using Oceananigans.Advection: FunctionStencil, divergence_smoothness
using Oceananigans.Advection: CrossUpwinding, SelfUpwinding, VelocityUpwinding

using Oceananigans.Advection: VectorInvariantCrossVerticalUpwinding, VectorInvariantSelfVerticalUpwinding, VectorInvariantVelocityVerticalUpwinding

using KernelAbstractions: @kernel, @index
using JLD2
using Random

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

function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.5)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δ  = 1 / sqrt(1 / Δx^2 + 1 / Δy^2)

    return  Base.Int(ceil(2 * Δt / (CFL / wave_speed * Δ)))
end

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
    
function baroclinic_adjustment_latlong(resolution, filename, FT::DataType = Float64; arch = GPU(), 
                                                     horizontal_closure = nothing,
                                                     momentum_advection = VectorInvariant(), 
                                                     background_νz = 1e-4,
                                                     xdirection = true)
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 50
    stop_time = 200days
    Δt = 2.5minutes

    if xdirection 
        grid = LatitudeLongitudeGrid(arch, FT;
                                    topology = (Periodic, Bounded, Bounded),
                                    size = (Ny, Ny, Nz), 
                                    longitude = (-10, 10),
                                    latitude = (-60, -40),
                                    z = (-Lz, 0),
                                    halo = (6, 6, 6))
    else
        grid = LatitudeLongitudeGrid(arch, FT;
                                     topology = (Flat, Bounded, Bounded),
                                     size = (Ny, Nz), 
                                     latitude = (-60, -40),
                                     z = (-Lz, 0),
                                     halo = (6, 6))
    end

    vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(FT; convective_κz = one(grid),
                                                                   convective_νz = zero(grid),
                                                                   background_κz = 1e-5,
                                                                   background_νz)

    closures = isnothing(horizontal_closure) ? vertical_closure : (vertical_closure, horizontal_closure)

    gravity = FT(Oceananigans.BuoyancyModels.g_Earth)
    max_Δt = 10minutes

    substeps = barotropic_substeps(max_Δt, grid, gravity)

    free_surface = SplitExplicitFreeSurface(FT; substeps)
    @info "running with substeps $substeps"
    @info "Building a model..."

    model = HydrostaticFreeSurfaceModel(; grid,
                                        coriolis = HydrostaticSphericalCoriolis(FT),
                                        buoyancy = BuoyancyTracer(),
                                        closure = closures,
                                        tracers = :b,
                                        momentum_advection,
                                        tracer_advection = WENO(FT),
                                        free_surface)

    @info "Built $model."

    gradient = "y"

    function ramp(λ, y, Δ)
        gradient == "x" && return min(max(0, λ / Δ + 1/2), 1)
        gradient == "y" && return min(max(0, (50 + y) / Δ + 1/2), 1)
    end

    # Parameters
    N² = 4e-6 # [s⁻²] buoyancy frequency / stratification

    Δy = 1.0 # degree
    Δb = 0.006
    ϵb = 1e-2 * Δb # noise amplitude

    Random.seed!(1234)
    bᵢ(x, y, z) = N² * z + Δb * ramp(x, y, Δy) + ϵb * randn()
    
    set!(model, b=bᵢ)

    #####
    ##### Simulation building
    #####

    simulation = Simulation(model; Δt, stop_time)

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

    reduced_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end

viscosity_name(clo::ScalarBiharmonicDiffusivity) = typeof(clo.ν)
viscosity_name(clo::ScalarDiffusivity)           = typeof(clo.ν)
viscosity_name(clo) = typeof(clo)
advection_name(adv) = getnamewrapper(adv.vorticity_scheme)

add_trailing_characters(name, trailing_character = "_weaker") = name * trailing_character

getupwindingscheme(::VectorInvariantCrossVerticalUpwinding)    = "f"
getupwindingscheme(::VectorInvariantSelfVerticalUpwinding)     = "p"
getupwindingscheme(::VectorInvariantVelocityVerticalUpwinding) = "s"

getname(s::VectorInvariant{N}) where N = "weno" * string(N * 2 - 1) * getupwindingscheme(s) * "$(s.vorticity_stencil isa VelocityStencil ? "V" : "D")"

function run_simulations(resolution, FT::DataType = Float64; trailing_character = "_weaker")

    advection_schemes   = []
    horizontal_closures = []
    names               = []

    # First five are the "Explicit" LES closures
    for i in 1:5
        push!(advection_schemes, VectorInvariant(vorticity_scheme = EnergyConservingScheme(), vertical_scheme = EnergyConservingScheme()))
    end

    hi2 = HorizontalScalarBiharmonicDiffusivity(FT; ν = geometric_νhb, discrete_form = true, parameters = 5days)
    hi3 = leith_viscosity(HorizontalFormulation(), FT)
    hi4 = leith_laplacian_viscosity(HorizontalFormulation(), FT)
    hi5 = smagorinski_viscosity(HorizontalFormulation(), FT)
    hi6 = QGLeith(FT)

    push!(horizontal_closures, hi2, hi3, hi4, hi5, hi6)
    push!(names, "bilap", "leith", "lapleith", "smag", "qgleith")

    # Next twelve are the "Implicit" LES closures
    for order in [5, 9]
        for upwinding_treatment in (CrossUpwinding(), SelfUpwinding(), VelocityUpwinding())
            for vorticity_stencil in (VelocityStencil(), DefaultStencil())
                push!(advection_schemes, VectorInvariant(; vorticity_scheme = WENO(FT; order), 
                                                           vorticity_stencil,
                                                           vertical_scheme = WENO(FT), 
                                                           upwinding_treatment))
                push!(horizontal_closures, nothing)
                push!(names, getname(advection_schemes[end]))
            end
        end
    end

    # Incorrect stencil usage
    for order in [5, 9]
        push!(advection_schemes, VectorInvariant(; vorticity_scheme = WENO(FT; order),  
                                                   vorticity_stencil = DefaultStencil(),
                                                   δU_stencil  = DefaultStencil(),
                                                   δV_stencil  = DefaultStencil(),
                                                   δu²_stencil = DefaultStencil(),
                                                   δv²_stencil = DefaultStencil(),
                                                   vertical_scheme = WENO(FT)))
                
        push!(horizontal_closures, nothing)
    end

    push!(names, "weno5pAllD", "weno9pAllD")
    
    # Flux form WENO schemes
    push!(advection_schemes, WENO(FT))
    push!(advection_schemes, WENO(FT; order = 9))
    push!(horizontal_closures, nothing, nothing)

    push!(names, "weno5Fl", "weno9Fl")

    push!(advection_schemes, VectorInvariant(vorticity_scheme = WENO(FT),            vertical_scheme = WENO(FT), multi_dimensional_stencil = true))
    push!(advection_schemes, VectorInvariant(vorticity_scheme = WENO(FT; order = 9), vertical_scheme = WENO(FT), multi_dimensional_stencil = true))
    push!(horizontal_closures, nothing, nothing)
    push!(names, "weno5MD", "weno9MD")

    @info "Running simulations with resolution $resolution" names

    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        @info "running simulation" name momentum_advection horizontal_closure
        baroclinic_adjustment_latlong(resolution, name, FT; momentum_advection, horizontal_closure, arch = GPU())
    end

    return nothing
end

function run_high_res_simulation(resolution; trailing_character = "_weaker")

    hi1 = nothing
    hi2 = leith_viscosity(HorizontalFormulation())
    
    vi1 = VectorInvariant()
    vi2 = VectorInvariant(vorticity_scheme = WENO(), vertical_scheme = WENO())

    advection_schemes   = [vi1, vi2]
    horizontal_closures = [hi1, hi2]

    names = ["leith", "weno5v"]
    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        @show name, momentum_advection, horizontal_closure
        baroclinic_adjustment_latlong(resolution, name; momentum_advection, horizontal_closure)
    end

    return nothing
end

run_2d_flat_simulation(resolution; trailing_character = "_weaker") = 
    baroclinic_adjustment_latlong(resolution, add_trailing_characters("2dsim", trailing_character);
                          horizontal_closure = leith_viscosity(HorizontalFormulation()),
                          xdirection = false)

function run_all(resolutions, FT::DataType = Float64; trailing_character = ["_weaker"])
    for (res, char) in zip(resolutions, trailing_character)
        run_simulations(res, FT; trailing_character = char)
    end
    run_high_res_simulation(1/50; trailing_character = "_fifty")
end

function run_all_rectilinear(resolutions; trailing_character = ["_weaker"])
    for (res, char) in zip(resolutions, trailing_character)
        run_simulations(res; trailing_character = char)
    end
    run_high_res_simulation(1/50; trailing_character = "_fifty")
end

include("Diagnostics/Diagnostics.jl")

using .Diagnostics

end
