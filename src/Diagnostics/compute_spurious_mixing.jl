using Statistics: mean
using Oceananigans.Fields: location

function compute_rpe_density(var)
    ze = calculate_z★_diagnostics(var[:b])

    εe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)
    αe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)

    zfield = HeightField(ze.grid)

    @info "computing resting and available potential energy density..."
    for t in 1:length(ze.times)
        @info "doing time $t"
        ρ = DensityOperation(var[:b][t])
        set!(εe[t], ze[t] * ρ)
        set!(αe[t], (zfield - ze[t]) * ρ)
    end

    return (; ze, εe, αe)
end

function calculate_RPE(st)
    RPE = Float64[]

    vol = VolumeField(st.εe[1].grid, location(st.εe[1]))

    for t in 1:length(st.ze.times)
        @info "doing time $t"
        push!(RPE, sum(interior(compute!(Field(st.εe[t] * vol)))))
    end

    return RPE
end

function calculate_APE(st)
    APE = Float64[]

    vol = VolumeField(st.αe[1].grid, location(st.αe[1]))

    for t in 1:length(st.ze.times)
        @info "doing time $t"
        push!(APE, sum(interior(compute!(Field(st.αe[t] * vol)))))
    end

    return APE
end

function calculate_KE(var)
    KE  = Float64[]

    vol = VolumeField(var[:u].grid)

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        v = var[:v][t]
        u = var[:u][t]
        ke = compute!(Field(@at (Center, Center, Center) (u^2 + v^2) * vol))

        push!(KE, 0.5 * sum(interior(ke)))
    end

    return KE
end

function calculate_Ω(var)
    Ω = Float64[]

    vol = VolumeField(var[:u].grid)

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        ζ = compute!(Field(VerticalVorticityOperation(var, t)^2 * vol))

        push!(Ω, sum(interior(ζ)))
    end

    return Ω
end

function calculate_N²(var)    
    N²avg = Float64[]

    vol = VolumeField(var[:u].grid)
    mean_vol = mean(vol)

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        N² = compute!(Field(StratificationOperation(var, t) * vol))

        push!(N²avg, mean(interior(N²)) / mean_vol)
    end

    return N²avg
end

function calculate_b_dissipation(var)    
    b_diss = Float64[]

    vol = VolumeField(var[:u].grid)
    mean_vol = mean(vol)

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        diss = compute!(Field(VerticalDissipationOperation(var, t) * vol))

        push!(b_diss, mean(interior(diss)) / mean_vol)
    end

    return b_diss
end

function calculate_deformation_radius(var)    
    LR = Float64[]

    area = AreaField(var[:u].grid)
    mean_area = mean(area)

    @info "computing resting and available potential energy density..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        R = compute!(Field(DeformationRadius(var, t) * area))

        push!(LR, mean(interior(R)) / mean_area)
    end

    return LR
end
