module Diagnostics

export all_fieldtimeseries, limit_timeseries!, propagate
export VolumeField, AreaField, MetricField, KineticEnergyField, time_average

using Oceananigans
using Oceananigans.AbstractOperations: AbstractOperation
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Fields: boundary_conditions, indices, location

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.OutputReaders: InMemoryFieldTimeSeries, OnDiskFieldTimeSeries

function FieldTimeSeries{LX, LY, LZ}(grid, times, FT=eltype(grid);
                                     indices = (:, :, :), 
                                     backend = InMemory(),
                                     path = nothing,
                                     name = nothing,
                                     boundary_conditions = nothing) where {LX, LY, LZ}

    Nt   = length(times)
    loc  = map(instantiate, (LX, LY, LZ))
    data = new_fieldtimeseries_data(FT, grid, loc, indices, Nt, path, name, backend)
    K = typeof(backend)
    return FieldTimeSeries{LX, LY, LZ, K}(data, grid, boundary_conditions, times, indices)
end

new_fieldtimeseries_data(FT, grid, loc, indices, Nt, path, name, ::OnDisk) = Oceananigans.OutputReaders.OnDiskData(path, name)

function new_fieldtimeseries_data(FT, grid, loc, indices, Nt, path, name, ::InMemory)
    space_size = total_size(grid, loc, indices)
    underlying_data = zeros(FT, arch, space_size..., Nt)
    return offset_data(underlying_data, grid, loc, indices)
end

set!(time_series::InMemoryFieldTimeSeries, f, index::Int) = set!(time_series[index], f)
set!(time_series::OnDiskFieldTimeSeries, f::AbstractOperation, index::Int) = set!(time_series, compute!(Field(f)), index)

# When we set! a OnDiskFieldTimeSeries we automatically write down the memory path
function set!(time_series::OnDiskFieldTimeSeries, f::AbstractArray, index::Int)
    path = time_series.data.path
    name = time_series.data.name
    jldopen(path, "a+") do file
        initialize_file!(file, name, time_series)
        maybe_write_property!(file, "timeseries/t/$index", time_series.times[index])
        maybe_write_property!(file, "timeseries/$(name)/$(index)", Array(parent(f)))
    end
end

function initialize_file!(file, name, time_series)
    maybe_write_property!(file, "serialized/grid", time_series.grid)
    maybe_write_property!(file, "timeseries/$(name)/serialized/location", location(time_series))
    maybe_write_property!(file, "timeseries/$(name)/serialized/indices", indices(time_series))
    maybe_write_property!(file, "timeseries/$(name)/serialized/boundary_conditions", boundary_conditions(time_series))
    return nothing
end

# Write property only if it does not already exist
function maybe_write_property!(file, property, data)
    try
        test = file[property]
    catch 
        file[property] = data
    end
end

import Oceananigans.AbstractOperations: restrict_index_for_interpolation

restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Face}) = from_index
restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Center}) = from_index

function propagate(fields...; func, path = nothing, name = nothing)

    fields_op = Tuple(field[1] for field in fields)
    operation = func(fields_op...)

    Nx, Ny, Nz = size(fields[1]) # take size of the first field which is the one that dominates

    if !(path isa Nothing) && !(name isa Nothing)
        field_output = FieldTimeSeries{location(operation)...}(fields[1].grid, fields[1].times; path, name, 
							       backend = OnDisk(), indices = (1:Nx, 1:Ny, 1:Nz))
    else    
        field_output = FieldTimeSeries{location(operation)...}(fields[1].grid, fields[1].times)
    end

    set!(field_output, operation, 1)

    for i in 2:length(field_output.times)
        @info "propagating on index $i"
        fields_op = retrieve_operand.(fields, i)
        operation = func(fields_op...)
        set!(field_output, operation, i)
    end

    return field_output
end

retrieve_operand(f::Number, i)          = f
retrieve_operand(f::Field, i)           = f
retrieve_operand(f::FieldTimeSeries, i) = f[i]

include("spurious_mixing.jl")
include("diagnostic_fields.jl")
include("integrated_diagnostics.jl")
include("spectra.jl")
include("compute_integrated_diagnostics.jl")

end
