using Statistics: mean

using BaroclinicAdjustment.Diagnostics: compute_rpe_density, 
                                        calculate_KE,
                                        calculate_APE,
                                        calculate_RPE

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

function compute_spectra(f::Dict)

    # Espec = Array{Spectrum}(undef, )
    # Ωspec = 

    # for time in 1:length(f[:u].times)


end

function calculate_diagnostics()
    file_prefix = ["bilap", "weno5vd", "leith", "lapleith", 
                   "smag", "weno5dd", "weno9",
                   "qgleith"]
    filenames = add_trailing_characters.(file_prefix)
    filenames = add_trailing_name.(filenames)

    energies   = Dict()
    spectras   = Dict()
    zonalmeans = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        fields    = all_fieldtimeseries(filename)
        # fields    = add_kinetic_energy_and_vorticity_timeseries!(fields)
        energy    = compute_spurious_mixing(fields)
        zonalmean = compute_zonal_mean(fields)
        # spectra   = compute_spectra(fields)

        energies[Symbol(prefix)] = energy
        # spectras[Symbol(prefix)] = spectra
        zonalmeans[Symbol(prefix)] = zonalmean

        jldsave("energies_"  * filename,"w") 
        jldsave("zonalmean_" * filename,"w") 
    end
    
    return nothing
end
