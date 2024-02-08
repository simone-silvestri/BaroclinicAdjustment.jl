using Statistics: mean
using Oceananigans
using Oceananigans.BoundaryConditions
using BaroclinicAdjustment
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

using Oceananigans.Grids: 
        cpu_face_constructor_x, 
        cpu_face_constructor_y, 
        cpu_face_constructor_z,
        pop_flat_elements,
        topology,
        architecture,
        metrics_precomputed

function with_precision(precision::DataType, old_grid::LatitudeLongitudeGrid)

    size = (old_grid.Nx, old_grid.Ny, old_grid.Nz)
    halo = (old_grid.Hx, old_grid.Hy, old_grid.Hz)

    topo = topology(old_grid)

    x = cpu_face_constructor_x(old_grid)
    y = cpu_face_constructor_y(old_grid)
    z = cpu_face_constructor_z(old_grid)

    # Remove elements of size and new_halo in Flat directions as expected by grid
    # constructor
    size     = pop_flat_elements(size, topo)
    new_halo = pop_flat_elements(halo, topo)

    new_grid = LatitudeLongitudeGrid(architecture(old_grid), precision;
                                     size = size, halo = new_halo,
                                     longitude = x, latitude = y, z = z, topology = topo,
                                     precompute_metrics = metrics_precomputed(old_grid),
                                     radius = old_grid.radius)

    return new_grid
end


add_trailing_name(name) = name * "_snapshots.jld2"

function compute_energy_diagnostics(f::Dict, iterations)      

    KE = calculate_KE(f)

    Etimeseries = compute_energy_timeseries(f)

    E = FieldTimeSeries{Center, Center, Center}(f[:u].grid, f[:u].times[iterations])

    for (i, t) in enumerate(iterations)
        set!(E[i], KineticEnergyOperation(f, t))
    end

    KEavg = Diagnostics.time_average(E)
    CUDA.@allowscalar begin
    	EKE = Diagnostics.propagate(E, KEavg; func = (e, ē) -> e - ē)
    end
    EKEavg = Diagnostics.time_average(EKE)

    return (; KE, EKEavg, Etimeseries)
end

function mean_eke(u, v, V)
    ū = mean(u, dims = 1)
    v̄ = mean(v, dims = 1)

    MEKE = 0.5 * (ū^2 + v̄^2) * V

    return sum(MEKE)
end

function eddy_eke(u, v, V)
    ū = mean(u, dims = 1)
    v̄ = mean(v, dims = 1)

    u′ = u - ū
    v′ = v - v̄

    EKE = 0.5 * (u′^2 + v′^2) * V

    return sum(EKE)
end

function mean_ape(b, V)
    
    b̄ = mean(b, dims = 1)
    B = mean(b̄, dims = 2)
    
    B̄  = b̄ - B
    N² = StratificationOperation(B)

    MAPE = 0.5 * B̄^2 / N² * V

    return sum(MAPE)
end

function eddy_ape(b, V)
    
    b̄ = mean(b, dims = 1)
    B = mean(b̄, dims = 2)
    
    B̄  = b̄ - B
    N² = StratificationOperation(B)

    b′ = b - B̄
    EAPE = 0.5 * b′^2 / N² * V

    return sum(EAPE)
end

function compute_energy_timeseries(f)

    grid = f[:u].grid
    Vᶜᶜᶜ = VolumeField(grid)
    V̄ᶜᶜᶜ = sum(Vᶜᶜᶜ, dims = 1)
    Vᵗ   = sum(interior(V̄ᶜᶜᶜ))

    MEKE = propagate(u, v, V̄ᶜᶜᶜ; func = (u, v, V) -> mean_eke(u, v, V) / Vᵗ)
    EKE  = propagate(u, v, Vᶜᶜᶜ; func = (u, v, V) -> eddy_eke(u, v, V) / Vᵗ)
    MAPE = propagate(b, V̄ᶜᶜᶜ;    func = (b, V)    -> mean_ape(b, V) / Vᵗ)
    EAPE = propagate(b, Vᶜᶜᶜ;    func = (b, V)    -> eddy_ape(b, V) / Vᵗ)

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

function compute_variances(f::Dict, fm, iterations; path = nothing)

    flucpath = path isa Nothing ? nothing : path * "fluc2_val.jld2"

    CUDA.@allowscalar begin
        u′ = propagate(f[:u], fm.U; func = (x, X) -> x - X, path = flucpath, name = "u") 
        v′ = propagate(f[:v], fm.V; func = (x, X) -> x - X, path = flucpath, name = "v")
        w′ = propagate(f[:w], fm.W; func = (x, X) -> x - X, path = flucpath, name = "w")
        b′ = propagate(f[:b], fm.B; func = (x, X) -> x - X, path = flucpath, name = "b")
    end
    
    u′² = time_average(propagate(u′; func = x -> x^2), iterations) 
    GC.gc(true)
    v′² = time_average(propagate(v′; func = x -> x^2), iterations)
    GC.gc(true)
    w′² = time_average(propagate(w′; func = x -> x^2), iterations)
    GC.gc(true)
    b′² = time_average(propagate(b′; func = x -> x^2), iterations)
    GC.gc(true)

    GC.gc(true)

    u′v′ = time_average(propagate(u′, v′; func = (x, y) -> x * y), iterations) 
    GC.gc(true)
    u′w′ = time_average(propagate(u′, w′; func = (x, y) -> x * y), iterations)
    GC.gc(true)
    v′w′ = time_average(propagate(v′, w′; func = (x, y) -> x * y), iterations)
    GC.gc(true)

    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y), iterations)
    GC.gc(true)
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y), iterations)
    GC.gc(true)
    w′b′ = time_average(propagate(w′, b′; func = (x, y) -> x * y), iterations)
    GC.gc(true)
    
    return (; u′², v′², w′², b′², u′v′, u′w′, v′w′, u′b′, v′b′, w′b′)
end

function compute_instability(fields, fm; path = nothing)

    flucpath = path isa Nothing ? nothing : path * "fluc3_val.jld2" 

    CUDA.@allowscalar begin
        b′ = propagate(fields[:b], fm.b̄; func = (x, X) -> x - X, path = flucpath, name = "u")
        u′ = propagate(fields[:u], fm.ū; func = (x, X) -> x - X, path = flucpath, name = "v")
        v′ = propagate(fields[:v], fm.v̄; func = (x, X) -> x - X, path = flucpath, name = "b")
    end
    
    u′b′ = time_average(propagate(u′, b′; func = (x, y) -> x * y))
    v′b′ = time_average(propagate(v′, b′; func = (x, y) -> x * y))

    return (; u′b′, v′b′)
end

function write_down_fields(fields::Dict; arch = nothing, path = nothing, chunk_size = nothing)
    new_fields = Dict()
    grid  = fields[:u].grid
    times = fields[:u].times

    grid = arch isa Nothing ? grid : on_architecture(arch, grid)
    path = path isa Nothing ? nothing : path * "all_values.jld2" 
    
    grid = with_precision(Float32, grid)

    if path isa Nothing
        u = FieldTimeSeries{Face, Center, Center}(grid,   times)
        v = FieldTimeSeries{Center, Face, Center}(grid,   times)
        w = FieldTimeSeries{Center, Center, Face}(grid,   times)
        b = FieldTimeSeries{Center, Center, Center}(grid, times) 
    else
        u = FieldTimeSeries{Face, Center, Center}(grid,   times; backend = OnDisk(), path, name = "u")
        v = FieldTimeSeries{Center, Face, Center}(grid,   times; backend = OnDisk(), path, name = "v")
        w = FieldTimeSeries{Center, Center, Face}(grid,   times; backend = OnDisk(), path, name = "w")
        b = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path, name = "b") 
    end

    utmp =  XFaceField(on_architecture(CPU(), grid))
    vtmp =  YFaceField(on_architecture(CPU(), grid))
    wtmp =  ZFaceField(on_architecture(CPU(), grid))
    btmp = CenterField(on_architecture(CPU(), grid))

    for t in eachindex(times)
      set!(utmp, Array(interior(fields[:u][t])))
      set!(vtmp, Array(interior(fields[:v][t])))
      set!(wtmp, Array(interior(fields[:w][t])))
      set!(btmp, Array(interior(fields[:b][t])))

      set!(u, utmp, t) 
      set!(v, vtmp, t)
      set!(w, wtmp, t)
      set!(b, btmp, t)
    end

    if !(chunk_size isa Nothing)
        u = FieldTimeSeries(path, "u"; backend = InMemory(; chunk_size))
        v = FieldTimeSeries(path, "v"; backend = InMemory(; chunk_size))
        w = FieldTimeSeries(path, "w"; backend = InMemory(; chunk_size))
        b = FieldTimeSeries(path, "b"; backend = InMemory(; chunk_size)) 
    end

    new_fields[:u] = u
    new_fields[:v] = v
    new_fields[:w] = w
    new_fields[:b] = b

    return new_fields
end

add_trailing_characters(string, trailing) = string * trailing

# fallback
calculate_diagnostics(test::TestCase, trailing_character; kwargs...) = 
    calculate_diagnostics([test.n], trailing_character; kwargs... )

function calculate_diagnostics(file_prefix::Vector = [], 
                               trailing_character = "_eigth";
                               arch = CPU(),
                               auxiliary_path = nothing,
                               src_path = nothing,
                               chunk_size = nothing)
    
    if !(auxiliary_path isa Nothing)
        try 
            files = readdir(auxiliary_path)
            for file in files
	           cmd = `rm $(auxiliary_path)/$(file)`
               run(cmd)
            end
        catch
            cmd = `mkdir $(auxiliary_path)`
            run(cmd)
        end
    end

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    postprocess = Dict()

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename) 
            @info "doing file " filename arch

            aux_arch = src_path isa Nothing ? arch : CPU()
            new_arch = src_path isa Nothing ? nothing : arch

            fields_previous = all_fieldtimeseries(filename; arch = aux_arch)
	        fields = write_down_fields(fields_previous; arch = new_arch, path = src_path, chunk_size)            

            lim = min(200, length(fields[:u].times))

            GC.gc(true)
            energy = compute_energy_diagnostics(fields, 50:lim)      
            GC.gc(true)
            spectra  = compute_spectra(fields, 50:lim)
            GC.gc(true)
            averages = compute_zonal_mean(fields, 50:lim)
            GC.gc(true)

            # From here I do not need!
            #= 
            variance = compute_variances(fields, averages, 50:lim; path = auxiliary_path)      
            GC.gc(true)
            instab   = compute_instability(fields, averages; path = auxiliary_path)      
            GC.gc(true)
            CUDA.@allowscalar begin
               N² = calculate_N²(fields)
            end
            GC.gc(true)
            enstrophy = calculate_Ω(fields)
            GC.gc(true)
	        postprocess[:enstrophy] = move_on_cpu(enstrophy)
	        postprocess[:stratif]   = move_on_cpu(N²)
	        postprocess[:variance]  = move_on_cpu(variance)
	        postprocess[:instab]    = move_on_cpu(instab)
            =#

	        postprocess[:energies]  = move_on_cpu(energy)
	        postprocess[:spectra]   = move_on_cpu(spectra)
	        postprocess[:mean]      = move_on_cpu(averages)

            write_file!(prefix * trailing_character * "_postprocess.jld2", postprocess)
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
move_on_cpu(field::AbstractField) = move_on_cpu(field, architecture(field))
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
