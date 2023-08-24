using Oceananigans.Utils
using Oceananigans.Fields: mean
using Oceananigans.Grids: halo_size
using Oceananigans.OutputReaders: OnDisk
using JLD2
using Oceananigans.Fields: default_indices

function checkpoint_fields(file)
    file = jldopen(file)
    grid = file["grid"]

    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    b = CenterField(grid)
    
    Hx, Hy, Hz = halo_size(grid)
    for (var, name) in zip((u, v, w, b), ("u", "v", "w", "b"))
        set!(var, file[name * "/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz])
        fill_halo_regions!(var)
    end

    return (; u, v, w, b)
end

assumed_location(var) = var == "u" ? (Face, Center, Center) : 
                        var == "v" ? (Center, Face, Center) : 
                        var == "w" ? (Center, Center, Face) : 
                        (Center, Center, Center)

function all_fieldtimeseries(filename, dir = nothing; arch = CPU(), variables = ("u", "v", "w", "b"), checkpointer = false)

    fields = Dict()

    if !(checkpointer)
        for var in variables
            fields[Symbol(var)] = FieldTimeSeries(filename, var; backend=OnDisk(), architecture=arch)
        end
    else
        files   = readdir(dir)
        myfiles = filter((x) -> x[1:length(filename)] == filename, files)
        numbers = parse.(Int, filter.(isdigit, myfiles))
        times   = numbers
        
        grid = jldopen(dir * myfiles[1])["grid"] 
        perm = sortperm(numbers)
        myfiles = myfiles[perm]
        for var in variables
            field = FieldTimeSeries{assumed_location(var)...}(grid, times)
            for (idx, file) in enumerate(myfiles)
                concrete_var = jldopen(dir * file)[var * "/data"]
                set!(field[idx], concrete_var)
            end

            fields[Symbol(var)] = field
        end
    end

    return fields
end

function limit_timeseries!(fields::Dict, times)
    new_fields = Dict()

    for (key, field) in fields
        new_fields[key] = limit_timeseries!(field, times)
    end

    return new_fields
end

function limit_timeseries!(field::FieldTimeSeries, times)

    loc = location(field)
    new_field = FieldTimeSeries{loc...}(field.grid, times)

    for (idx, time) in enumerate(field.times)
        id2 = findfirst(isequal(time), times)
        if !isnothing(id2)
            set!(new_field[id2], field[idx])
        end
    end

    return new_field
end

function reduce_output_size!(old_file_name, new_file_name; limit_to = 20, variables = ("u", "v", "w", "b"))

    var_dict = all_fieldtimeseries(old_file_name; variables)
    times = var_dict[Symbol(variables[1])].times[end - limit_to:end]
    var_dict = limit_timeseries!(var_dict, times)

    jldsave(new_file_name, vars = var_dict)
end

function add_kinetic_energy_from_timeseries!(fields::Dict, iterations = 1:length(fields[:u].times))

    E = FieldTimeSeries{Face, Center, Center}(fields[:u].grid, fields[:u].times[iterations])

    for (i, t) in enumerate(iterations)
        set!(E[i], KineticEnergyOperation(fields, i))
    end

    fields[:E] = E

    return nothing
end

function add_kinetic_energy_and_vorticity_timeseries!(fields::Dict)

    E = FieldTimeSeries{Center, Center, Center}(fields[:u].grid, fields[:u].times)
    ζ = FieldTimeSeries{Face, Face, Center}(fields[:u].grid, fields[:u].times)
    u = Field{Face, Center, Center}(fields[:u].grid)
    v = Field{Center, Face, Center}(fields[:u].grid)

    for t in 1:length(E.times)
        set!(u, fields[:u][t])
        set!(v, fields[:v][t])

        fill_halo_regions!((u, v))

        set!(ζ[t], VerticalVorticityOperation((; u, v)))
        set!(E[t], KineticEnergyOperation((; u, v)))
    end

    fields[:E] = E
    fields[:ζ] = ζ

    return nothing
end

function kinetic_energy(u::FieldTimeSeries, v::FieldTimeSeries)

    energy = Float64[]
    vol = VolumeField(u.grid)

    for (i, time) in enumerate(u.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(u.times[end]))"
        ke = u[i]^2 + v[i]^2
        push!(energy, sum(ke * vol))
    end

    return energy
end

function ACC_transport(u::FieldTimeSeries; stride = 1, start_time = 1, end_time = length(u.times))

    transport = Float64[]
    vol = VolumeField(u.grid)

    for i in start_time:stride:end_time
        @info "integrating time $(prettytime(u.times[i])) of $(prettytime(u.times[end]))"
        push!(transport, sum(u[i] * vol, dims = (2, 3))[1, 1, 1])
    end

    return transport
end

function heat_content(b::FieldTimeSeries; stride = 1, start_time = 1, end_time = length(b.times))

    heat = Float64[]
    vol = VolumeField(b.grid)

    for i in start_time:stride:end_time
        @info "integrating time $(prettytime(b.times[i])) of $(prettytime(b.times[end]))"
        push!(heat, sum(b[i] * vol))
    end

    return heat
end

using Oceananigans.Operators: Δzᶜᶠᶜ

function calculate_eulerian_MOC(v::Field)
        
    v̄ = mean(v, dims = 1)

    ψ = Field((Nothing, Face, Face), v.grid)

    for k in 2:v.grid.Nz
        dz =  Δzᶜᶠᶜ(1, 1, k-1, v.grid)
        for j in 1:size(v.grid, 2)
            ψ[1, j, k] = ψ[1, j, k - 1] + v̄[1, j, k - 1] * dz
        end
    end

    return ψ
end

function calculate_eulerian_MOC(v::FieldTimeSeries)

    v̄ = time_average(v)
    ψ = calculate_MOC(v̄)
    
    return ψ
end

function time_average(field::FieldTimeSeries, iterations = 1:length(field.times))
    avg = similar(field[1])
    fill!(avg, 0)

    for t in iterations
        avg .+= field[t] ./ length(iterations)
    end

    return avg
end

function calculate_fluctuations!(fields::Dict, variables)

    for var in variables

        field_avg  = time_average(fields[var])
        field_fluc = FieldTimeSeries{location(fields[var])...}(fields[var].grid, fields[var].times)

        for t in 1:length(fields[var].times)
            fluc = fields[var][t] - field_avg
            set!(field_fluc[t], fluc)
        end

        fields[Symbol(var, :fluc)] = field_fluc
    end

    return nothing
end
