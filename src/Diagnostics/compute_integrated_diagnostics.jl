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

function calculate_KE(var; compute_projection = true)
    KE  = Float64[]
    
    BTKE = []
    BCKE = []

    vol = VolumeField(var[:u].grid)

    @info "computing kinetic energy..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        v = var[:v][t]
        u = var[:u][t]
        
        ke = compute!(Field(@at (Center, Center, Center) (u^2 + v^2) * vol))

        push!(KE, 0.5 * sum(interior(ke)))

        if compute_projection && t == length(var[:u].times)
            btke = sum(interior(ke),         dims = 3) 
            bcke = sum(interior(ke) .- btke, dims = 3)
    
            push!(BTKE, btke)
            push!(BCKE, bcke)
        end
    end

    if compute_projection
        return KE, BTKE, BCKE
    else
        return KE
    end
end

function calculate_Ω(var)
    Ω = Float64[]

    vol = VolumeField(var[:u].grid)

    @info "computing enstrophy..."
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

    @info "computing stratification..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        N² = compute!(Field(StratificationOperation(var, t) * vol))

        push!(N²avg, mean(interior(N²)) / mean_vol)
    end

    return N²avg
end

using Oceananigans.BoundaryConditions

function calculate_b_budget(var)    
    dissipation     = Float64[]
    time_derivative = Float64[]
    
    vol = VolumeField(var[:u].grid)
    fill_halo_regions!(vol)

    b1  = CenterField(var[:b].grid)
    b2  = CenterField(var[:b].grid)

    @info "computing buoyancy budget..."
    for t in 1:2:length(var[:b].times)-1
        @info "doing time $t"
        dt    = (var[:b].times[t+1] - var[:b].times[t]) 

        set!(b1, var[:b][t])
        set!(b2, var[:b][t+1])

        fill_halo_regions!((b1, b2))
        diss1 = compute!(Field(∂z(b1)^2 * vol))
        diss2 = compute!(Field(∂z(b2)^2 * vol))

        b12 = sum(interior(compute!(Field(b1^2 * vol))))
        b22 = sum(interior(compute!(Field(b2^2 * vol))))
        db2dt = (b22 - b12) / dt
        push!(dissipation, 1e-5 * (sum(interior(diss1)) + sum(interior(diss2))) / 2)
        push!(time_derivative,     db2dt)
    end

    return (; dissipation, time_derivative)
end

function calculate_deformation_radius(var)    
    LR = Float64[]

    area = AreaField(var[:u].grid)
    mean_area = mean(area)

    @info "computing deformation radius..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        R = compute!(Field(DeformationRadius(var, t) * area))

        push!(LR, mean(interior(R)) / mean_area)
    end

    return LR
end

@inline function _slope_operation(i, j, k, grid, B)

    ∂yb = ℑzᵃᵃᶠ(i, j, k, grid, ∂yᶜᶠᶜ, B)
    ∂zb = max(1e-20, ℑyᵃᶠᵃ(i, j, k, grid, ∂zᶜᶜᶠ, B))

    return (∂yb / ∂zb)
end

function calculate_slope(var)    
    slope = Float64[]
    b = CenterField(var[:b].grid)

    vol = VolumeField(var[:u].grid)
    mean_vol = mean(vol)

    @info "computing deformation radius..."
    for t in 1:length(var[:u].times)
        @info "doing time $t"
        set!(b, var[:b][t])
        fill_halo_regions!(b)
        B = mean(b, dims = 1)
        S = KernelFunctionOperation{Nothing, Face, Face}(_slope_operation, grid, B)
        push!(slope, mean(filter(x -> x != 0, interior(S))) / mean_vol)
    end

    return slope
end
