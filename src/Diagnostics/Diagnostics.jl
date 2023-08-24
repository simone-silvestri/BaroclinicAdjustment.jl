module Diagnostics

export all_fieldtimeseries, limit_timeseries!, propagate
export VolumeField, AreaField, MetricField, KineticEnergyField, time_average

using Oceananigans
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll

function propagate(fields...; func)

    fields_op = Tuple(field[1] for field in fields)
    operation = func(fields_op...)

    field_output = FieldTimeSeries{location(operation)...}(fields[1].grid, fields[1].times)
    
    set!(field_output[1], operation)

    for i in 2:length(field_output.times)
        @info "propagating on index $i"
        fields_op = retrieve_operand.(fields, i)
        operation = func(fields_op...)
        set!(field_output[i], operation)
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
