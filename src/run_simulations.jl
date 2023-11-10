function run_simulations(resolution, FT::DataType = Float64; φ₀ = -50, simulation = baroclinic_adjustment_latlong, trailing_character = "")

    advection_schemes, horizontal_closures, names = testcases(FT)

    @info "Running simulations with resolution $resolution" names

    names = add_trailing_characters.(names, Ref(trailing_character))

    for (momentum_advection, horizontal_closure, name) in zip(advection_schemes, horizontal_closures, names)
        @info "running simulation" name momentum_advection horizontal_closure
        simulation(resolution, name, FT; momentum_advection, φ₀, horizontal_closure)
    end

    return nothing
end

function run_all_simulations(; simulation = baroclinic_adjustment_latlong, 
                               FT::DataType = Float64; 
                               φ₀ = -50) 
                                                          
    run_simulations(1/8,  FT; φ₀, simulation; trailing_character = "_eight")
    run_simulations(1/16, FT; φ₀, simulation; trailing_character = "_sixteen")
    run_simulations(1/32, FT; φ₀, simulation; trailing_character = "_thirtytwo")
    
    return nothing
end