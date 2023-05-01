using Statistics: mean
using Oceananigans
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: compute_rpe_density, 
                                        calculate_KE,
                                        calculate_APE,
                                        calculate_RPE,
                                        calculate_Ω,
                                        calculate_N²,
                                        calculate_deformation_radius,
                                        compute_spectra

using JLD2

add_trailing_name(name) = name * "_snapshots.jld2"

function compute_spurious_mixing(f::Dict)

    KE  = calculate_KE(f)
    aux = compute_rpe_density(f)

    APE = calculate_APE(aux)
    RPE = calculate_RPE(aux)

    return (; KE, APE, RPE)
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

function calculate_diagnostics(trailing_character = "_weaker")
    file_prefix = ["weno5vd", "leith", "lapleith", "bilap",
                   "smag", "weno5dd", "weno5vv", "weno9", "weno9dd",
                   "qgleith", "highres"]
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    energies  = Dict()
    spectras  = Dict()
    defradii  = Dict()

    enstrophies = Dict()
    stratif     = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())

            GC.gc()
            energy    = compute_spurious_mixing(fields)
            enstrophy = calculate_Ω(fields)
            N²        = calculate_N²(fields)
            defradius = calculate_deformation_radius(fields)
            # spectra   = compute_spectra(fields)

            energies[Symbol(prefix)]    = energy
            enstrophies[Symbol(prefix)] = enstrophy
            stratif[Symbol(prefix)]     = N²
            defradii[Symbol(prefix)]    = defradius
            # spectras[Symbol(prefix)] = spectra
        end
    end

    write_file!("enstrophies" * trailing_character * ".jld2", enstrophies)
    write_file!("stratif" *     trailing_character * ".jld2", stratif)
    write_file!("energies" *    trailing_character * ".jld2", energies) 
    write_file!("defradii" *    trailing_character * ".jld2", defradii)
    # write_file!("spectra" * trailing_character * ".jld2", spectras) 
    
    return nothing
end

function write_file!(name, var) 
    jldopen(name,"w") do f
        for (key, value) in var
            f[string(key)] = value
        end
    end
end