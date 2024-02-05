module Postprocess

export calculate_diagnostics 

using BaroclinicAdjustment
using Oceananigans

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

include("calculate_diagnostics.jl")

end # module Utils
