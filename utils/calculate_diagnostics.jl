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
                                        compute_spectra

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

function generate_names()

    names = []

    # First five are the "Explicit" LES closures
    push!(names, "bilap", "leith", "lapleith", "smag", "qgleith")

    # Next twelve are the "Implicit" LES closures
    for order in [5, 9]
        for upwinding_treatment in (CrossAndSelfUpwinding(), OnlySelfUpwinding(), VelocityUpwinding())
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

function save_contours(trailing_character = "_weaker")
    file_prefix = generate_names()

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

function calculate_diagnostics(trailing_character = "_weaker")
    file_prefix = generate_names()

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    # energies  = Dict()
    # vardiss   = Dict()
    spectras  = Dict()

    # enstrophies = Dict()
    # stratif     = Dict()
    # budgetB     = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())

            GC.gc()
            # energy    = compute_energy_diagnostics(fields)
            # budget    = calculate_b_budget(fields)
            # enstrophy = calculate_Ω(fields)
            # N²        = calculate_N²(fields)
            if length(fields[:u].times) > 90
                spectra = compute_spectra(fields, 80:100)
                spectras[(Symbol(prefix), :90)]  = spectra
            end
            if length(fields[:u].times) > 190
                spectra = compute_spectra(fields, 180:200)
                spectras[(Symbol(prefix), :190)] = spectra
            end
            if length(fields[:u].times) > 290
                spectra = compute_spectra(fields, 280:300)
                spectras[(Symbol(prefix), :290)] = spectra
            end
            if length(fields[:u].times) > 390
                spectra = compute_spectra(fields, 380:400)
                spectras[(Symbol(prefix), :390)] = spectra
            end
            # energies[Symbol(prefix)]    = energy
            # enstrophies[Symbol(prefix)] = enstrophy
            # stratif[Symbol(prefix)]     = N²
            # budgetB[Symbol(prefix)]     = budget
        end
    end

    # write_file!("enstrophies" * trailing_character * ".jld2", enstrophies)
    # write_file!("stratif" *     trailing_character * ".jld2", stratif)
    # write_file!("energies" *    trailing_character * ".jld2", energies)
    # write_file!("vardiss" *     trailing_character * ".jld2", vardiss) 
    # write_file!("budgetB" *     trailing_character * ".jld2", budgetB)
    write_file!("spectra" *     trailing_character * ".jld2", spectras)

    return nothing
end

function write_file!(name, prefix, var) 
    if isfile(name)
        jldopen(name,"r+") do f
            for (key, value) in var
                f[prefix * "/" * string(key)] = value
            end
        end
    else
        jldopen(name,"w") do f
            for (key, value) in var
                f[prefix * "/" * string(key)] = value
            end
        end
    end
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
