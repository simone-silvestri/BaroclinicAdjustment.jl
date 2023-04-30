using GLMakie
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: VerticalVorticity, KineticEnergy, DeformationRadius

using Oceananigans
using JLD2

add_trailing_name(name) = name * "_snapshots.jld2"

function record_video!(name, fig, iter) 
    GLMakie.record(fig, name * ".mp4", 1:200, framerate = 11) do i
        @info "step $i"; 
        iter[] = i; 
    end
end

function produce_videos(trailing_character = "_weaker")
    file_prefix = ["weno5vd", "leith", "lapleith", "bilap",
                   "smag", "weno5dd", "weno5vv", "weno9", "weno9dd",
                   "qgleith", "highres"]
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    iter = Observable(1)

    for (prefix, filename) in zip(file_prefix, filenames)
        f = all_fieldtimeseries(filename)

        @info "vorticity video $filename"    
        ζ = @lift(interior(VerticalVorticity(f, $iter), :, :, 50))

        fig = Figure()
        ax  = Axis(fig[1, 1])
        heatmap!(ax, ζ, colorrange = (-1e-5, 1e-5))
        record_video!("vort_" * prefix * trailing_character * ".jld2", fig, iter)

        @info "kinetic energy video $filename"    
        E = @lift(interior(KineticEnergy(f, $iter), :, :, 50))

        fig = Figure()
        ax  = Axis(fig[1, 1])
        heatmap!(ax, E, colorrange = (0, 1.0))
        record_video!("energy_" * prefix * trailing_character * ".jld2", fig, iter)

        @info "deformation radius video $filename"    
        R = @lift(interior(DeformationRadius(f, $iter), :, :, 1))

        fig = Figure()
        ax  = Axis(fig[1, 1])
        heatmap!(ax, R, colorrange = (0, 1e4))
        record_video!("Lr_" * prefix * trailing_character * ".jld2", fig, iter)
    end
end