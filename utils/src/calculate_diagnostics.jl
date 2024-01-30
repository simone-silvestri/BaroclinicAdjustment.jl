using Revise
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

    KEavg = Diagnostics.time_average(E)
    CUDA.@allowscalar begin
    	EKE = Diagnostics.propagate(E, KEavg; func = (e, ē) -> e - ē, path = "auxiliaries/energies.jld2", name = "old_EKE")
    end
    EKEavg = Diagnostics.time_average(EKE)

    return (; KE, EKEavg, Etimeseries)
end

function compute_energy_timeseries(f)
    ū = propagate(f[:u]; func = x -> mean(x, dims = 1), path = "auxiliaries/mean.jld2", name = "ū")
    v̄ = propagate(f[:v]; func = x -> mean(x, dims = 1), path = "auxiliaries/mean.jld2", name = "v̄")
    b̄ = propagate(f[:b]; func = x -> mean(x, dims = 1), path = "auxiliaries/mean.jld2", name = "b̄")

    u′ = propagate(f[:u], ū; func = (x, X) -> x - X, path = "auxiliaries/variance.jld2", name = "u′")
    v′ = propagate(f[:v], v̄; func = (x, X) -> x - X, path = "auxiliaries/variance.jld2", name = "v′")
    b′ = propagate(f[:b], b̄; func = (x, X) -> x - X, path = "auxiliaries/variance.jld2", name = "b′")

    B = propagate(f[:b]; func = x -> mean(x, dims = 2))

    B̄  = propagate(b̄, B; func = (b̄, B) -> b̄ - B, path = "auxiliaries/energies.jld2", name = "B̄")
    N² = propagate(B; func = B -> StratificationOperation(B), path = "auxiliaries/energies.jld2", name = "N²")

    MEKE = propagate(ū , v̄ ; func = (u, v)  -> mean(0.5 * (u^2 + v^2)), path = "auxiliaries/energies.jld2", name = "MEKE")
    EKE  = propagate(u′, v′; func = (u, v)  -> mean(0.5 * (u^2 + v^2)), path = "auxiliaries/energies.jld2", name = "EKE")
    MAPE = propagate(B̄ , N²; func = (B, N²) -> mean(0.5 * B^2 / N²),    path = "auxiliaries/energies.jld2", name = "MAPE")
    EAPE = propagate(b′, N²; func = (b, N²) -> mean(0.5 * b^2 / N²),    path = "auxiliaries/energies.jld2", name = "EAPE")

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
    CUDA.@allowscalar begin
        u′ = propagate(f[:u], fm.U; func = (x, X) -> x - X, path = "auxiliaries/new_variance.jld2", name = "u′")
        v′ = propagate(f[:v], fm.V; func = (x, X) -> x - X, path = "auxiliaries/new_variance.jld2", name = "v′")
        w′ = propagate(f[:w], fm.W; func = (x, X) -> x - X, path = "auxiliaries/new_variance.jld2", name = "w′")
        b′ = propagate(f[:b], fm.B; func = (x, X) -> x - X, path = "auxiliaries/new_variance.jld2", name = "b′")
    end
    
    u′² = time_average(propagate(u′; func = x -> x^2, path = "auxiliaries/new_variance.jld2", name = "u′²"), iterations)
    v′² = time_average(propagate(v′; func = x -> x^2, path = "auxiliaries/new_variance.jld2", name = "v′²"), iterations)
    w′² = time_average(propagate(w′; func = x -> x^2, path = "auxiliaries/new_variance.jld2", name = "w′²"), iterations)
    b′² = time_average(propagate(b′; func = x -> x^2, path = "auxiliaries/new_variance.jld2", name = "b′²"), iterations)

    u′v′ = time_average(propagate(u′, v′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "u′v′"), iterations)
    u′w′ = time_average(propagate(u′, w′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "u′w′"), iterations)
    v′w′ = time_average(propagate(v′, w′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "v′w′"), iterations)

    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "u′b′"), iterations)
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "v′b′"), iterations)
    w′b′ = time_average(propagate(w′, b′; func = (x, y) -> x * y, path = "auxiliaries/new_variance.jld2", name = "w′b′"), iterations)
    
    return (; u′², v′², w′², b′², u′v′, u′w′, v′w′, u′b′, v′b′, w′b′)
end

function compute_instability(fields, fm, iterations)

    CUDA.@allowscalar begin
        b′ = propagate(fields[:b], fm.b̄; func = (x, X) -> x - X, path = "auxiliaries/new_new_variance.jld2", name = "b′")
        u′ = propagate(fields[:u], fm.ū; func = (x, X) -> x - X, path = "auxiliaries/new_new_variance.jld2", name = "u′")
        v′ = propagate(fields[:v], fm.v̄; func = (x, X) -> x - X, path = "auxiliaries/new_new_variance.jld2", name = "v′")
    end
    
    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y, path = "auxiliaries/new_new_variance.jld2", name = "u′b′"), iterations)
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y, path = "auxiliaries/new_new_variance.jld2", name = "v′b′"), iterations)

    return (; u′b′, v′b′)
end

function write_down_fields(fields::Dict)
    new_fields = Dict()
    grid  = fields[:u].grid
    times = fields[:u].times
    path = fields[:u].data.path

    u = FieldTimeSeries{Face, Center, Center}(grid, times; backend = OnDisk(), path, name = "u")
    v = FieldTimeSeries{Center, Face, Center}(grid, times; backend = OnDisk(), path, name = "v")
    w = FieldTimeSeries{Center, Center, Face}(grid, times; backend = OnDisk(), path, name = "w")
    b = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path, name = "b")

    new_fields[:u] = u
    new_fields[:v] = v
    new_fields[:w] = w
    new_fields[:b] = b

    return new_fields
end

function calculate_diagnostics(trailing_character = "_weaker", file_prefix = generate_names())

    iter = [5, 14, 19, 2]
    file_prefix = file_prefix[iter]
    
    try 
        readdir("./auxiliaries/")
    catch
        cmd = `mkdir auxiliaries`
        run(cmd)
    end

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    postprocess = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename) 
	    
            arch = GPU()

            try run(`rm ./auxiliaries/fields.jld2`); catch; end

            @info "doing file " filename arch
            fields_previous = all_fieldtimeseries(filename; arch)

            fields = write_down_fields(fields_previous)


            lim = min(200, length(fields[:u].times))

            GC.gc(true)
            energy    = compute_energy_diagnostics(fields, 50:lim)
            GC.gc(true)
            enstrophy = calculate_Ω(fields)
            GC.gc(true)
            N²        = calculate_N²(fields)
            GC.gc(true)
            spectra   = compute_spectra(fields, 50:lim)
            GC.gc(true)
            averages  = compute_zonal_mean(fields, 50:lim)
            GC.gc(true)
            variance  = compute_variances(fields, averages, 50:lim)
            GC.gc(true)
            instab    = compute_instability(fields, averages, 50:lim)
            GC.gc(true)

	        postprocess[:energies]  = move_on_cpu(energy)
	        postprocess[:enstrophy] = move_on_cpu(enstrophy)
	        postprocess[:stratif]   = move_on_cpu(N²)
	        postprocess[:spectra]   = move_on_cpu(spectra)
	        postprocess[:mean]      = move_on_cpu(averages)
	        postprocess[:variance]  = move_on_cpu(variance)
	        postprocess[:instab]    = move_on_cpu(instab)

            write_file!(prefix * trailing_character * "_postprocess.jld2", postprocess)

            try run(`mv ./auxiliaries/mean.jld2 ./auxiliaries/$(file_prefix)_mean.jld2`); catch; end
            try run(`mv ./auxiliaries/variance.jld2 ./auxiliaries/$(file_prefix)_variance.jld2`); catch; end
            try run(`mv ./auxiliaries/energies.jld2 ./auxiliaries/$(file_prefix)_energies.jld2`); catch; end
            try run(`mv ./auxiliaries/new_variance.jld2 ./auxiliaries/$(file_prefix)_new_variance.jld2`); catch; end
            try run(`mv ./auxiliaries/new_new_variance.jld2 ./auxiliaries/$(file_prefix)_new_new_variance.jld2`); catch; end
        end
    end

    return nothing
end

# In case we want to postprocess on the GPU we need to save on CPU
using Oceananigans.Fields: location, AbstractField
using Oceananigans.Grids: on_architecture, architecture

move_on_cpu(array::CuArray) = Array(array)
move_on_cpu(array) = array

move_on_cpu(fields::NamedTuple) = 
    NamedTuple{propertynames(fields)}(map(move_on_cpu, fields))

move_on_cpu(fields::Tuple) = map(move_on_cpu, fields)
move_on_cpu(field::AbstractField)   = move_on_cpu(field, architecture(field))
move_on_cpu(field::FieldTimeSeries) = move_on_cpu(field, architecture(field))
move_on_cpu(field, ::CPU) = field

function move_on_cpu(fields::FieldTimeSeries) 
    grid = on_architecture(CPU(), fields.grid)
    cpu_fields = FieldTimeSeries{location(fields)...}(grid, fields.times)
    for i in 1:length(fields.times)
       set!(cpu_fields[i], move_on_cpu(fields[i]))
    end
    return cpu_fields
end

function move_on_cpu(field, ::GPU)
    grid = on_architecture(CPU(), field.grid)
    cpu_field = Field{location(field)...}(grid)
    set!(cpu_field, field)
    return cpu_field
end
