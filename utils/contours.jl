using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes, znodes, nodes

using GLMakie
using JLD2
using GeoMakie

using Statistics

using CairoMakie
CairoMakie.activate!()

function enhance_vars(var::Matrix, rep = 1)
    var = repeat(var,  inner = (rep, 1))
    return repeat(var, inner = (1, rep))
end

enhance_vars(var::Vector, rep = 1) = repeat(var,  inner = rep)

qg08 = all_fieldtimeseries("qgleith_eight_new_snapshots.jld2", "./")
wn08 = all_fieldtimeseries("weno9pV_eight_new_snapshots.jld2", "./")
qg12 = all_fieldtimeseries("qgleith_sixteen_new_snapshots.jld2", "./")
wn12 = all_fieldtimeseries("weno9pV_sixteen_new_snapshots.jld2", "./")
# qg16 = all_fieldtimeseries("test_qgleith.jld2", "./")
# wn16 = all_fieldtimeseries("test_weno9pV.jld2", "./")
# qg16 = all_fieldtimeseries("qgleith_thirtytwo_new_snapshots.jld2", "./")
# wn16 = all_fieldtimeseries("weno9pV_thirtytwo_new_snapshots.jld2", "./")

ζqg08 = Diagnostics.VerticalVorticity(qg08, 40)
ζwn08 = Diagnostics.VerticalVorticity(wn08, 40)
ζqg12 = Diagnostics.VerticalVorticity(qg12, 40)
ζwn12 = Diagnostics.VerticalVorticity(wn12, 40)
ζqg16 = Diagnostics.VerticalVorticity(qg16, 40)
ζwn16 = Diagnostics.VerticalVorticity(wn16, 40)

λ08, φ08, _ = enhance_vars.(Array.(nodes(ζqg08)))
λ12, φ12, _ = enhance_vars.(Array.(nodes(ζqg12)))
λ16, φ16, _ = enhance_vars.(Array.(nodes(ζqg16)))

λ08 = λ08[1:end-1]
λ12 = λ12[1:end-1]
λ16 = λ16[1:end-1]

φ08 = φ08[1:end-1]
φ12 = φ12[1:end-1]
φ16 = φ16[1:end-1]

ζqg08 = enhance_vars(Array(interior(ζqg08, :, :, 45)))
ζwn08 = enhance_vars(Array(interior(ζwn08, :, :, 45)))
ζqg12 = enhance_vars(Array(interior(ζqg12, :, :, 45)))
ζwn12 = enhance_vars(Array(interior(ζwn12, :, :, 45)))
ζqg16 = enhance_vars(Array(interior(ζqg16, :, :, 45)))
ζwn16 = enhance_vars(Array(interior(ζwn16, :, :, 45)))

g08 = LatitudeLongitudeGrid(size = (160, 160, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))
g12 = LatitudeLongitudeGrid(size = (240, 240, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))
g16 = LatitudeLongitudeGrid(size = (320, 320, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))

lonlims08 = extrema(λ08)
latlims08 = extrema(φ08)
lonlims12 = extrema(λ12)
latlims12 = extrema(φ12)
lonlims16 = extrema(λ16)
latlims16 = extrema(φ16)

kwargs08 = (; )
kwargs12 = (; )
kwargs16 = (; )

ticks(a::Vector) = (a, [L"\text{%$b}" for b in a])

# possible cmaps:
# diverging_bkr_55_10_c35_n256
# diverging_linear_bjr_30_55_c53_n256

surfargs = (colorrange = (-6e-5, 6e-5), colormap = :bwr)

fig = Figure(fontsize=28, resolution = (930, 800) .* 1.5)

ga = fig[1, 1:3] = GridLayout()
gb = fig[2, 1:3] = GridLayout()

ax = Axis(ga[1, 1]; title = L"\text{QG2, 14 kilometers}", ylabel = L"\text{Latitude}", yticks = ([-60, -50, -40], [L"60^\text{o} \text{S}", L"50^\text{o} \text{S}", L"40^\text{o} \text{S}"]), xticks = ([-10, -5, 0, 5, 10], ["", "", "", "", ""]), kwargs08...)            
heatmap!(ax, λ08, φ08, ζqg08; surfargs...)
ax = Axis(ga[1, 2]; title = L"\text{QG2, 7 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], ["", "", ""]), kwargs12...)            
heatmap!(ax, λ12, φ12, ζqg12; surfargs...)
ax = Axis(ga[1, 3]; title = L"\text{QG2, 3.5 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], ["", "", ""]), kwargs16...)            # 
heatmap!(ax, λ16, φ16, ζqg16; surfargs...)

ax = Axis(gb[1, 1]; title = L"\text{W9V, 14 kilometers}", ylabel = L"\text{Latitude}", yticks = ([-60, -50, -40], [L"60^\text{o} \text{S}", L"50^\text{o} \text{S}", L"40^\text{o} \text{S}"]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs08...)            # 
heatmap!(ax, λ08, φ08, ζwn08; surfargs...)
ax = Axis(gb[1, 2]; title = L"\text{W9V, 7 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs12...)            # 
heatmap!(ax, λ12, φ12, ζwn12; surfargs...)
ax = Axis(gb[1, 3]; title = L"\text{W9V, 3.5 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs16...)            # 
hm = heatmap!(ax, λ16, φ16, ζwn16; surfargs...)

cb = Colorbar(fig[1:2, 4], hm, height = Relative(3/4), ticks = ([-4e-5, -2e-5, 0, 2e-5, 4e-5], [L"-4", L"-2", L"0", L"2", L"4"]), label = L"\text{Vertical vorticity } \left[s^{-1} \cdot 10^{-5} \right]",)

rowgap!(ga, 0)
rowgap!(gb, 0)

# CairoMakie.save("vorticity_cont.png", fig, px_per_unit = 8)
CairoMakie.save("vorticity_cont.eps", fig, px_per_unit = 8)