module BaroclinicAdjustment

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.Utils: getnamewrapper
using Oceananigans.Advection: VelocityStencil, DefaultStencil
using JLD2

include("horizontal_visc.jl")
include("qg_leith_viscosity.jl")
include("outputs.jl")

function barotropic_substeps(Δt, grid, gravitational_acceleration)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δ  = 1 / sqrt(1 / Δx^2 + 1 / Δy^2)

    return  Base.Int(ceil(2 * Δt / (0.7 / wave_speed * Δ)))
end
    
function baroclinic_adjustment(resolution, filename; horizontal_closure = nothing, momentum_advection = VectorInvariant())

    # Architecture
    arch = GPU()
    
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = Base.Int(20 / resolution)
    Nz = 50
    stop_time = 200days
    Δt = 10minutes

    grid = LatitudeLongitudeGrid(arch;
                                topology = (Periodic, Bounded, Bounded),
                                size = (Ny, Ny, Nz), 
                                longitude = (-10, 10),
                                latitude = (-60, -40),
                                z = (-Lz, 0),
                                halo = (6, 6, 6))

    coriolis = HydrostaticSphericalCoriolis()

    Δy = 1000kilometers / Ny
    vertical_closure = VerticalScalarDiffusivity(ν=1e-5, κ=1e-6)

    closures = (vertical_closure, horizontal_closure)

    substeps = barotropic_substeps(10minutes, grid, Oceananigans.BuoyancyModels.g_Earth)

    @info "running with $substeps substeps"
    free_surface = SplitExplicitFreeSurface(; substeps)

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

    """
    Linear ramp from 0 to 1 between -Δy/2 and +Δy/2.

    For example:

    y < y₀           => ramp = 0
    y₀ < y < y₀ + Δy => ramp = y / Δy
    y > y₀ + Δy      => ramp = 1
    """

    gradient = "y"

    function ramp(λ, y, Δ)
        gradient == "x" && return min(max(0, λ / Δ + 1/2), 1)
        gradient == "y" && return min(max(0, (50 + y) / Δ + 1/2), 1)
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

    reduced_outputs!(simulation, filename)

    run!(simulation)

    return nothing
end

viscosity_name(clo::ScalarBiharmonicDiffusivity) = typeof(clo.ν)
viscosity_name(clo::ScalarDiffusivity)           = typeof(clo.ν)
viscosity_name(clo) = typeof(clo)
advection_name(adv) = getnamewrapper(adv.vorticity_scheme)

add_trailing_characters(name) = name * "_larger"

function run_eight_degree_simulations()

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
    names = add_trailing_characters.(names)
    
    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        baroclinic_adjustment(1/8, name; momentum_advection, horizontal_closure)
    end

    return nothing
end

function run_high_res_simulation()

    vi1 = VectorInvariant()

    hi4 = leith_laplacian_viscosity(HorizontalFormulation(), C_vort = 2.0, C_div = 2.0)

    advection_schemes   = [vi1]
    horizontal_closures = [hi4]

    names = ["highres"]
    names = add_trailing_characters.(names)

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        baroclinic_adjustment(1/50, name; momentum_advection, horizontal_closure)
    end

    return nothing
end

include("Diagnostics/Diagnostics.jl")

using .Diagnostics

add_trailing_name(name) = name * "_snapshots.jld2"

function calculate_diagnostics()
    file_prefix = ["bilap", "weno5vd", "leith", "lapleith", 
                   "smag", "weno5dd", "weno5vv", "weno9", "highres"]
    filenames = add_trailing_characters.(file_prefix)
    filenames = add_trailing_name.(filenames)

    energies   = Dict()
    spectras   = Dict()
    zonalmeans = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        fields    = all_fieldtimeseries(filename)
        energy    = compute_spurious_mixing(fields)
        zonalmean = compute_zonal_mean(fields)
        # spectra = compute_spectra(fields)

        energies[Symbol(prefix)] = energy
        # spectras[Symbol(prefix)] = spectra
        zonalmeans[Symbol(prefix)] = zonalmean
    end

    jldopen("energies.jld2","w") do f
        for (key, value) in energies
            f[string(key)] = value
        end
    end

    jldopen("zonalmeans.jld2","w") do f
        for (key, value) in zonalmeans
            f[string(key)] = value
        end
    end

    return nothing
end

end
