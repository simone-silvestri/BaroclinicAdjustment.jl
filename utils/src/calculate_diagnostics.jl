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
                                        compute_spectra,
                                        VolumeField,
                                        KineticEnergyOperation,
                                        StratificationOperation,
                                        propagate,
                                        time_average

using JLD2
using CUDA

using Oceananigans.Advection: CrossAndSelfUpwinding, OnlySelfUpwinding, VelocityUpwinding
using Oceananigans.Advection: VelocityStencil, DefaultStencil

using BaroclinicAdjustment.Diagnostics: VerticalVorticityOperation, KineticEnergyOperation, DensityOperation

add_trailing_name(name) = name * "_snapshots.jld2"

function compute_energy_diagnostics(f::Dict, iterations)

    KE = calculate_KE(f)

    Etimeseries = compute_energy_timeseries(f)

    E = FieldTimeSeries{Face, Center, Center}(f[:u].grid, f[:u].times[iterations])

    for (i, t) in enumerate(iterations)
        set!(E[i], KineticEnergyOperation(f, t))
    end

    KEavg  = Diagnostics.time_average(E)
    EKE    = Diagnostics.propagate(E, KEavg, func = (e, ē) -> e - ē)
    EKEavg = Diagnostics.time_average(EKE)

    return (; KE, EKEavg, Etimeseries)
end

function compute_energy_timeseries(f)
    ū = propagate(f[:u]; func = x -> mean(x, dims = 1))
    v̄ = propagate(f[:v]; func = x -> mean(x, dims = 1))
    b̄ = propagate(f[:b]; func = x -> mean(x, dims = 1))

    u′ = propagate(f[:u], ū; func = (x, X) -> x - X)
    v′ = propagate(f[:v], v̄; func = (x, X) -> x - X)
    b′ = propagate(f[:b], b̄; func = (x, X) -> x - X)

    MEKE = propagate(ū , v̄ ; func = (u, v) -> mean(0.5 * (u^2 + v^2)))
    EKE  = propagate(u′, v′; func = (u, v) -> mean(0.5 * (u^2 + v^2)))
    MAPE = propagate(b̄     ; func =  B     -> mean(0.5 * B^2 / StratificationOperation(B)))
    EAPE = propagate(b′, b̄ ; func = (b, B) -> mean(0.5 * b^2 / StratificationOperation(B)))

    return (; MEKE, EKE, MAPE, EAPE)
end

function compute_zonal_mean(f::Dict, iterations)
    ū = time_average(f[:u], iterations)
    v̄ = time_average(f[:v], iterations)
    w̄ = time_average(f[:w], iterations)
    b̄ = time_average(f[:b], iterations)

    U = mean(ū, dims = 1) 
    V = mean(v̄, dims = 1)
    W = mean(w̄, dims = 1)
    B = mean(b̄, dims = 1)

    B★ = mean(B, dims = 2)
    
    return (; U, V, W, B, B★, ū, v̄, w̄, b̄)
end

function compute_variances(f::Dict, fm, iterations)
    u′ = propagate(f[:u], fm.U; func = (x, X) -> x - X)
    v′ = propagate(f[:v], fm.V; func = (x, X) -> x - X)
    w′ = propagate(f[:w], fm.W; func = (x, X) -> x - X)
    b′ = propagate(f[:b], fm.B; func = (x, X) -> x - X)
    
    u′² = time_average(propagate(u′; func = x -> x^2), iterations)
    v′² = time_average(propagate(v′; func = x -> x^2), iterations)
    w′² = time_average(propagate(w′; func = x -> x^2), iterations)
    b′² = time_average(propagate(b′; func = x -> x^2), iterations)

    u′v′ = time_average(propagate(u′, v′; func = (x, y) -> x * y), iterations)
    u′w′ = time_average(propagate(u′, w′; func = (x, y) -> x * y), iterations)
    v′w′ = time_average(propagate(v′, w′; func = (x, y) -> x * y), iterations)

    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y), iterations)
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y), iterations)
    w′b′ = time_average(propagate(w′, b′; func = (x, y) -> x * y), iterations)
    
    return (; u′², v′², w′², b′², u′v′, u′w′, v′w′, u′b′, v′b′, w′b′)
end

function compute_instability(fields, fm, iterations)

    b′ = propagate(fields[:b], fm.b̄; func = (x, X) -> x - X)
    u′ = propagate(fields[:u], fm.ū; func = (x, X) -> x - X)
    v′ = propagate(fields[:v], fm.v̄; func = (x, X) -> x - X)

    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y), iterations)
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y), iterations)

    return (; u′b′, v′b′)
end

function calculate_diagnostics(trailing_character = "_weaker", file_prefix = generate_names())

    iter = [5, 14, 19, 2]
    file_prefix = file_prefix[iter]
    
    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    postprocess = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename) && !(prefix == "qgleith" && trailing_character == "_eight_new")

            @info "doing file " filename
            fields = all_fieldtimeseries(filename; arch = CPU())

            GC.gc()
            energy    = compute_energy_diagnostics(fields, 50:200)
            GC.gc()
            enstrophy = calculate_Ω(fields)
            GC.gc()
            N²        = calculate_N²(fields)
            GC.gc()
            spectra   = compute_spectra(fields, 50:200)
            GC.gc()
            averages  = compute_zonal_mean(fields, 50:200)
            GC.gc()
            variance  = compute_variances(fields, averages, 50:200)
            GC.gc()
            instab    = compute_instability(fields, averages, 50:200)
            GC.gc()

            postprocess[:energies]  = energy
            postprocess[:enstrophy] = enstrophy
            postprocess[:stratif]   = N²
            postprocess[:spectra]   = spectra
            postprocess[:mean]      = averages
            postprocess[:variance]  = variance
            postprocess[:instab]    = instab

            write_file!(prefix * trailing_character * "_postprocess.jld2", postprocess)
        end
    end

    return nothing
end

# In case we want to postprocess on the GPU we need to save on CPU
using Oceananigans.Fields: location, AbstractField
using Oceananigans.Grids:  on_architecture

move_on_cpu(array::CuArray) = Array(array)
move_on_cpu(array) = array

move_on_cpu(fields::NamedTuple) = 
    NamedTuple{propertynames(fields)}(map(move_on_cpu, fields))

move_on_cpu(fields::Tuple) = map(move_on_cpu, fields)
move_on_cpu(field::AbstractField) = move_on_cpu(field, architecuture(field))

move_on_cpu(field, ::CPU) = field

function move_on_cpu(field, ::GPU)
    grid = on_architecture(CPU(), field.grid)
    cpu_field = Field{location(field)...}(grid)
    set!(cpu_field, field)
    return cpu_field
end
