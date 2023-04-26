using Statistics: mean

function compute_rpe_density(var)
    ze = calculate_z★_diagnostics(var[:b])

    εe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)
    αe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)

    zfield = HeightField(ze.grid)

    @info "computing resting and available potential energy density..."
    for t in 1:length(ze.times)
        @info "doing time $t"
        ρ = DensityField(var[:b][t])
        set!(εe[t], ze[t] * ρ)
        set!(αe[t], (- zfield - ze[t]) * ρ)
    end

    εeavg = mean(εe[length(ze.times)] - εe[1], dims = 1)

    return (; ze, εe, αe, εeavg)
end

function calculate_RPE(st)
    RPE = Float64[]

    for t in 1:length(st.ze.times)
        @info "doing time $t"
        push!(RPE, mean(st.εe[t]))
    end

    return RPE
end

function calculate_APE(st)
    APE = Float64[]

    for t in 1:length(st.ze.times)
        @info "doing time $t"
        push!(APE, mean(st.αe[t]))
    end

    return APE
end

function calculate_KE(var)
    KE  = Float64[]

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        v = var[:v][t]
        u = var[:u][t]
        ke = @at (Center, Center, Center) u^2 + v^2

        push!(KE, mean(ke))
    end

    return KE
end
