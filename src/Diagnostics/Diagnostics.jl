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

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: FieldTimeSeries

import Oceananigans.AbstractOperations: restrict_index_for_interpolation

restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Face})   = from_index
restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Center}) = from_index

restrict_index_for_interpolation(from_index, ::Type{Face},   ::Type{Nothing}) = from_index
restrict_index_for_interpolation(from_index, ::Type{Center}, ::Type{Nothing}) = from_index

restrict_index_for_interpolation(from_index, ::Type{Nothing}, ::Type{Nothing}) = from_index

function propagate(fields...; func, path = nothing, name = nothing)

    fields_op = Tuple(field[1] for field in fields)
    operation = func(fields_op...)

    field1 = Field(operation)

    if !(path isa Nothing) && !(name isa Nothing)
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

retrieve_operand(f::Number, i)          = f
retrieve_operand(f::Field, i)           = f
retrieve_operand(f::FieldTimeSeries, i) = f[i]

include("spurious_mixing.jl")
include("diagnostic_fields.jl")
include("integrated_diagnostics.jl")
include("spectra.jl")
include("compute_integrated_diagnostics.jl")

end
