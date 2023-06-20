using Statistics: mean
using Oceananigans
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters, getname
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: compute_rpe_density, 
                                        calculate_KE,
                                        calculate_APE,
                                        calculate_RPE,
                                        calculate_Ω,
                                        calculate_N²,
                                        calculate_z★_diagnostics,
                                        calculate_deformation_radius,
                                        calculate_b_budget,
                                        compute_spectra,
                                        VolumeField

using JLD2

using Oceananigans.Advection: CrossAndSelfUpwinding, OnlySelfUpwinding, VelocityUpwinding
using Oceananigans.Advection: VelocityStencil, DefaultStencil

using BaroclinicAdjustment.Diagnostics: VerticalVorticityOperation, KineticEnergyOperation, DensityOperation

add_trailing_name(name) = name * "_snapshots.jld2"

function compute_energy_diagnostics(f::Dict)

    KE, BTKE, BCKE  = calculate_KE(f)
    aux = compute_rpe_density(f)

    APE = calculate_APE(aux)
    RPE = calculate_RPE(aux)

    return (; KE, BTKE, BCKE, APE, RPE)
end

function compute_zonal_mean(f::Dict)
    ū = FieldTimeSeries{Nothing, Center, Center}(f[:u].grid, f[:u].times)
    b̄ = FieldTimeSeries{Nothing, Center, Center}(f[:u].grid, f[:u].times)

    v̄ = FieldTimeSeries{Nothing, Face, Center}(f[:u].grid, f[:u].times)
    w̄ = FieldTimeSeries{Nothing, Center, Face}(f[:u].grid, f[:u].times)

    for time in 1:length(f[:u].times)
        set!(ū[time], mean(f[:u][time], dims = 1))
        set!(v̄[time], mean(f[:v][time], dims = 1))
        set!(w̄[time], mean(f[:w][time], dims = 1))
        set!(b̄[time], mean(f[:b][time], dims = 1))
    end

    return (; ū, v̄, w̄, b̄)
end

function save_contours(trailing_character = "_weaker", file_prefix = generate_names())

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    enstrophy  = Dict()
    buoyancy   = Dict()
    yzbuoyancy = Dict()
    vorticity  = Dict()
    energy     = Dict()
    referenpe  = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())

            ze = calculate_z★_diagnostics(fields[:b], 400)
            ρ  = DensityOperation(fields[:b][400])
            εe = compute!(Field(ze * ρ))
    
            GC.gc()
            ζ  = compute!(Field(VerticalVorticityOperation(fields, 400)))
            ζ² = compute!(Field(ζ^2))

            B = compute!(Field(Average(fields[:b][400], dims = 1)))

            yzbuoyancy[String(prefix)] = (interior(B, 1, :, :), )
            KE = compute!(Field(KineticEnergyOperation(fields, 400)))

            enstrophy[Symbol(prefix)] = (interior(ζ², :, :, 25), interior(ζ², :, :, 50))
            vorticity[Symbol(prefix)] = (interior(ζ, :, :, 25), interior(ζ, :, :, 50))
            buoyancy[Symbol(prefix)]  = (interior(fields[:b][400], :, :, 25), interior(fields[:b][400], :, :, 50))
            energy[Symbol(prefix)]    = (interior(KE, :, :, 25), interior(KE, :, :, 50))
            referenpe[Symbol(prefix)] = (interior(εe, :, :, 25), interior(εe, :, :, 50))
        end
    end

    write_file!("contours" * trailing_character * ".jld2", "enstrophy",   enstrophy )
    write_file!("contours" * trailing_character * ".jld2", "vorticity",   vorticity )
    write_file!("contours" * trailing_character * ".jld2", "buoyancy",    buoyancy  )    
    write_file!("contours" * trailing_character * ".jld2", "yzbuoyancy2",  yzbuoyancy)
    write_file!("contours" * trailing_character * ".jld2", "energy",      energy    )
    write_file!("contours" * trailing_character * ".jld2", "referencepe", referenpe )
end

function calculate_fluxes(var)

    fluxes = []        
    vol = VolumeField(var[:u].grid)
    mean_vol = mean(interior(vol))

    for t in 1:length(var[:u].times)
        @info "doing time $t"
        b  = var[:b][t]
        w  = var[:w][t]        
        v  = var[:v][t]
        WB = mean(w * b, dims = 1)
        VB = mean(v * b, dims = 1)

        push!(fluxes, (; wb = mean(WB * vol)[1, 1, 1] / mean_vol, vb = mean(VB * vol)[1, 1, 1] / mean_vol))
    end

    return fluxes
end

function calculate_diagnostics(trailing_character = "_weaker", file_prefix = generate_names())

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    # energies    = Dict()
    # vardiss     = Dict()
    # spectras    = Dict()
    # enstrophies = Dict()
    # stratif     = Dict()
    # budgetB     = Dict()
    fluxes      = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())
            
            # spectras[Symbol(prefix)] = []

            GC.gc()
            # energy    = compute_energy_diagnostics(fields)
            # budget    = calculate_b_budget(fields)
            # enstrophy = calculate_Ω(fields)
            # N²        = calculate_N²(fields)
            fluxe     = calculate_fluxes(fields)
            # for range in (45:55, 95:105, 145:155, 190:200)
            #     if length(fields[:u].times) >= range[end]
            #         spectra = compute_spectra(fields, range)
            #         push!(spectras[Symbol(prefix)], spectra)
            #     end
            # end
            # energies[Symbol(prefix)]    = energy
            # enstrophies[Symbol(prefix)] = enstrophy
            # stratif[Symbol(prefix)]     = N²
            # budgetB[Symbol(prefix)]     = budget
            fluxes[Symbol(prefix)]      = fluxe
        end
    end

    # write_file!("enstrophies" * trailing_character * ".jld2", enstrophies)
    # write_file!("stratif" *     trailing_character * ".jld2", stratif)
    # write_file!("energies" *    trailing_character * ".jld2", energies)
    # write_file!("vardiss" *     trailing_character * ".jld2", vardiss) 
    # write_file!("budgetB" *     trailing_character * ".jld2", budgetB)
    # write_file!("spectra" *     trailing_character * ".jld2", spectras)
    write_file!("fluxes" *     trailing_character * ".jld2", fluxes)

    return nothing
end
