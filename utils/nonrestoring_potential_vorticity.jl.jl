using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes, znodes, nodes

using GLMakie
using JLD2

using Statistics


PVw = Float64[]
PVo = Float64[]
varw = all_fieldtimeseries("weno9pV_sixteen_norest_snapshots.jld2", "./"; with_halos = true, newpath = "test1.jld2")
varo = all_fieldtimeseries("omp25_sixteen_norest_snapshots.jld2", "./";   with_halos = true, newpath = "test2.jld2")

vol = Diagnostics.VolumeField(varw[:u].grid)

vtot = sum(interior(vol))
for time in 1:101
    @info "doing time $time"
    push!(PVw, sum(interior(Diagnostics.PotentialVorticity(varw, time)) .* interior(vol)) ./ vtot)
    push!(PVo, sum(interior(Diagnostics.PotentialVorticity(varo, time)) .* interior(vol)) ./ vtot)
end


