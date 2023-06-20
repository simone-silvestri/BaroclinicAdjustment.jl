function run_simulations(resolution, FT::DataType = Float64; φ₀ = -50, simulation = baroclinic_adjustment_latlong, trailing_character = "_weaker")

    advection_schemes, horizontal_closures, names = testcases(FT)

    @info "Running simulations with resolution $resolution" names

    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        @info "running simulation" name momentum_advection horizontal_closure
        simulation(resolution, name, FT; momentum_advection, φ₀, horizontal_closure)
    end

    return nothing
end

function run_high_res_simulation(resolution, FT::DataType = Float64; φ₀ = -50, simulation = baroclinic_adjustment_latlong, trailing_character = "_weaker")

    hi1 = nothing
    hi2 = QGLeith(FT)
    
    vi1 = VectorInvariant(vorticity_scheme = EnergyConservingScheme(FT), vertical_scheme = EnergyConservingScheme(FT))
    vi2 = VectorInvariant(vorticity_scheme = WENO(FT), vertical_scheme = WENO(FT))

    advection_schemes   = [vi2, vi1]
    horizontal_closures = [hi1, hi2]

    names = ["leith", "weno5pV"]
    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        @show name, momentum_advection, horizontal_closure
        simulation(resolution, name, FT; momentum_advection, φ₀, horizontal_closure)
    end

    return nothing
end

function run_all(resolutions, FT::DataType = Float64; φ₀ = -50, simulation = baroclinic_adjustment_latlong, trailing_character = ["_weaker"])
    for (res, char) in zip(resolutions, trailing_character)
        run_simulations(res, FT; simulation, φ₀, trailing_character = char)
    end
    run_high_res_simulation(1/50, FT; simulation, φ₀, trailing_character = "_fifty")
end

run_all_resolutions(simulation, FT::DataType = Float64; φ₀ = -50) = run_all([1/4, 1/8, 1/16], FT; φ₀, simulation, trailing_character = ["_quarter", "_eight", "_sixteen"])