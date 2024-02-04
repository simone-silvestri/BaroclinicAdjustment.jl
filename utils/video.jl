using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes, znodes, nodes

using GLMakie
using JLD2

using Statistics

var = all_fieldtimeseries("weno9pV_sixteen_new_postprocess.jld2", "./")

iter = Observable(1)
Nt = length(var[:u].times)
# @show maximum(wn16[:u][end])
# @show minimum(wn16[:u][200])

KE = @lift(interior(Diagnostics.Potential(var, $iter), :, :, 45))

kwargs = (; colormap = :berlin, colorrange = (-6e-5, 6e-5))

fig = Figure()
ax  = Axis(fig[1, 1])
heatmap!(ax, KE; kwargs...)

GLMakie.record(fig, "video.mp4", 1:Nt) do i
    @info "doing iteration $i"
    iter[] = i
end
