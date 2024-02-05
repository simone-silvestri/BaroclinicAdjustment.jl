module Diagnostics

export all_fieldtimeseries, limit_timeseries!, propagate
export VolumeField, AreaField, MetricField, KineticEnergyField, time_average

using Oceananigans
using Oceananigans.AbstractOperations: AbstractOperation
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Grids: architecture
using Oceananigans.Fields: boundary_conditions, indices, location, total_size, offset_data
using Oceananigans.OutputReaders: InMemoryFieldTimeSeries, OnDiskFieldTimeSeries

import Oceananigans.Fields: set!, AbstractField
import Oceananigans.AbstractOperations: restrict_index_for_interpolation

restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Face})   = from_index
restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Center}) = from_index

restrict_index_for_interpolation(from_index, ::Type{Face},   ::Type{Nothing}) = from_index
restrict_index_for_interpolation(from_index, ::Type{Center}, ::Type{Nothing}) = from_index

eestrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Nothing}) = from_index

function propagate(fields...; func, path = nothing, name = nothing)

    fields_op = Tuple(retrieve_operand(field, 1) for field in fields)
    operation = func(fields_op...)

    if operation isa Number
       field_output = []
       for i in 2:length(fields[1].times)
          fields_op = retrieve_operand.(fields, i)
          operation = func(fields_op...)
          push!(field_output, operation)
       end
       return field_output
    end

    field1 = Field(operation)

    if !(path isa Nothing)
        field_output = FieldTimeSeries{location(field1)...}(fields[1].grid, fields[1].times; path, name, backend = OnDisk(), indices = field1.indices)
    else    
        field_output = FieldTimeSeries{location(field1)...}(fields[1].grid, fields[1].times, indices = field1.indices)
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

import Oceananigans.Fields: set!
using Oceananigans.OutputReaders: initialize_file!, maybe_write_property!
using Oceananigans.OutputReaders: InMemoryFieldTimeSeries

set!(time_series::InMemoryFieldTimeSeries, f, index::Int) = set!(time_series[index], f)
 
# When we set! a OnDiskFieldTimeSeries we automatically write down the memory path
function set!(time_series::OnDiskFieldTimeSeries, f::AbstractOperation, index::Int)
    path = time_series.path
    name = time_series.name
    
    b = compute!(Field(f))
    jldopen(path, "a+") do file
        initialize_file!(file, name, time_series)
        maybe_write_property!(file, "timeseries/t/$index", time_series.times[index])
        maybe_write_property!(file, "timeseries/$(name)/$(index)", Array(parent(b)))
    end
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
