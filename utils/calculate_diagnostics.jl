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
                                        calculate_deformation_radius,
                                        calculate_b_dissipation,
                                        compute_spectra

using JLD2

using Oceananigans.Advection: CrossUpwinding, SelfUpwinding, VelocityUpwinding
using Oceananigans.Advection: VelocityStencil, DefaultStencil

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

function generate_names()

    names = []

    # First five are the "Explicit" LES closures
    push!(names, "bilap", "leith", "lapleith", "smag", "qgleith")

    # Next twelve are the "Implicit" LES closures
    for order in [5, 9]
        for upwinding_treatment in (CrossUpwinding(), SelfUpwinding(), VelocityUpwinding())
            for vorticity_stencil in (VelocityStencil(), DefaultStencil())
                adv = VectorInvariant(; vorticity_scheme = WENO(; order), 
                                        vorticity_stencil,
                                        vertical_scheme = WENO(), 
                                        upwinding_treatment)
                push!(names, getname(adv))
            end
        end
    end

    push!(names, "weno5pAllD", "weno9pAllD")
    push!(names, "weno5Fl", "weno9Fl")
    push!(names, "weno5MD", "weno9MD")

    return names
end

function calculate_diagnostics(trailing_character = "_weaker")
    file_prefix = generate_names()

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    energies  = Dict()
    vardiss   = Dict()
    spectras  = Dict()

    enstrophies = Dict()
    stratif     = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())

            GC.gc()
            energy    = compute_energy_diagnostics(fields)
            variance  = calculate_b_dissipation(fields)
            enstrophy = calculate_Ω(fields)
            N²        = calculate_N²(fields)
            if length(fields[:u].times) > 390
                spectra = compute_spectra(fields, 380:400)
                spectras[Symbol(prefix)] = spectra
            end

            energies[Symbol(prefix)]    = energy
            vardiss[Symbol(prefix)]     = variance
            enstrophies[Symbol(prefix)] = enstrophy
            stratif[Symbol(prefix)]     = N²
        end
    end

    write_file!("enstrophies" * trailing_character * ".jld2", enstrophies)
    write_file!("stratif" *     trailing_character * ".jld2", stratif)
    write_file!("energies" *    trailing_character * ".jld2", energies)
    write_file!("vardiss" *     trailing_character * ".jld2", vardiss) 
    write_file!("spectra" *     trailing_character * ".jld2", spectras)
    
    return nothing
end

function write_file!(name, var) 
    if isfile(name)
        jldopen(name,"r+") do f
            for (key, value) in var
                f[string(key)] = value
            end
        end
    else
        jldopen(name,"w") do f
            for (key, value) in var
                f[string(key)] = value
            end
        end
    end
    return nothing
end
