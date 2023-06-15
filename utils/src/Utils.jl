module Utils

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

    # Last one are the "Strays" (Incorrect Stencils, Flux Form, Multi-Dimensional)
    push!(names, "weno5pAllD", "weno9pAllD")
    push!(names, "weno5Fl", "weno9Fl")
    push!(names, "weno5MD", "weno9MD")

    return names
end

include("calculate_diagnostics.jl")
include("surface_videos.jl")

end # module Utils