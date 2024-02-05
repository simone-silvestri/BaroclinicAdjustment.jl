using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes, znodes, nodes

using GLMakie
using JLD2

using Statistics

# using CairoMakie
# CairoMakie.activate!()

function enhance_vars(var::Matrix, rep = 1)
    var = repeat(var,  inner = (rep, 1))
    return repeat(var, inner = (1, rep))
end

enhance_vars(var::Vector, rep = 1) = repeat(var,  inner = rep)

qg08 = all_fieldtimeseries("ebs_eight_snapshots.jld2", "./")
wn08 = all_fieldtimeseries("weno9pV_eight_new_snapshots.jld2", "./")
qg16 = all_fieldtimeseries("test.jld2", "./")
wn16 = all_fieldtimeseries("qgleith_eight_new_snapshots.jld2", "./")
qg32 = all_fieldtimeseries("weno9pAllD_eight_new_snapshots.jld2", "./")
wn32 = all_fieldtimeseries("upwind_eight_snapshots.jld2", "./")

ζqg08 = Diagnostics.DeformationRadius(qg08, 119)
ζwn08 = Diagnostics.DeformationRadius(wn08, 119)
ζqg16 = Diagnostics.DeformationRadius(qg16, 119)
ζwn16 = Diagnostics.DeformationRadius(wn16, 119)
ζqg32 = Diagnostics.DeformationRadius(qg32, 119)
ζwn32 = Diagnostics.DeformationRadius(wn32, 119)
surfargs = (colorrange = (0, 8e-6), colormap = :magma)

level = 1

λ08, φ08, _ = enhance_vars.(Array.(nodes(ζqg08)))
λ16, φ16, _ = enhance_vars.(Array.(nodes(ζqg16)))
λ32, φ32, _ = enhance_vars.(Array.(nodes(ζqg32)))

λ08 = λ08[1:end-1]
λ16 = λ16[1:end-1]
λ32 = λ32[1:end-1]

φ08 = φ08[1:end-1]
φ16 = φ16[1:end-1]
φ32 = φ32[1:end-1]

ζqg08 = enhance_vars(Array(interior(ζqg08, :, :, level)))
ζwn08 = enhance_vars(Array(interior(ζwn08, :, :, level)))
ζqg16 = enhance_vars(Array(interior(ζqg16, :, :, level)))
ζwn16 = enhance_vars(Array(interior(ζwn16, :, :, level)))
ζqg32 = enhance_vars(Array(interior(ζqg32, :, :, level)))
ζwn32 = enhance_vars(Array(interior(ζwn32, :, :, level)))

g08 = LatitudeLongitudeGrid(size = (160, 160, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))
g16 = LatitudeLongitudeGrid(size = (320, 320, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))
g32 = LatitudeLongitudeGrid(size = (320, 320, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))

lonlims08 = extrema(λ08)
latlims08 = extrema(φ08)
lonlims16 = extrema(λ16)
latlims16 = extrema(φ16)
lonlims32 = extrema(λ32)
latlims32 = extrema(φ32)

kwargs08 = (; )
kwargs16 = (; )
kwargs32 = (; )

ticks(a::Vector) = (a, [L"\text{%$b}" for b in a])

# possible cmaps:
# diverging_bkr_55_10_c35_n256
# diverging_linear_bjr_30_55_c53_n256

fig = Figure(fontsize=28, resolution = (930, 800) .* 1.5)

ga = fig[1, 1:3] = GridLayout()
gb = fig[2, 1:3] = GridLayout()

ax = Axis(ga[1, 1]; title = L"\text{QG2, 14 kilometers}", ylabel = L"\text{Latitude}", yticks = ([-60, -50, -40], [L"60^\text{o} \text{S}", L"50^\text{o} \text{S}", L"40^\text{o} \text{S}"]), xticks = ([-10, -5, 0, 5, 10], ["", "", "", "", ""]), kwargs08...)            
heatmap!(ax, λ08, φ08, ζqg08; surfargs...)
ax = Axis(ga[1, 2]; title = L"\text{QG2, 7 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], ["", "", ""]), kwargs16...)            
heatmap!(ax, λ16, φ16, ζqg16; surfargs...)
ax = Axis(ga[1, 3]; title = L"\text{QG2, 3.5 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], ["", "", ""]), kwargs32...)            # 
heatmap!(ax, λ32, φ32, ζqg32; surfargs...)

ax = Axis(gb[1, 1]; title = L"\text{W9V, 14 kilometers}", ylabel = L"\text{Latitude}", yticks = ([-60, -50, -40], [L"60^\text{o} \text{S}", L"50^\text{o} \text{S}", L"40^\text{o} \text{S}"]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs08...)            # 
heatmap!(ax, λ08, φ08, ζwn08; surfargs...)
ax = Axis(gb[1, 2]; title = L"\text{W9V, 7 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs16...)            # 
hm = heatmap!(ax, λ16, φ16, ζwn16; surfargs...)
ax = Axis(gb[1, 3]; title = L"\text{W9V, 3.5 kilometers}", yticks = ([-60, -50, -40], ["", "", ""]), xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0", L"5^\text{o} \text{E}"]), xlabel = L"\text{Longitude}", kwargs32...)            # 
heatmap!(ax, λ32, φ32, ζwn32; surfargs...)

cb = Colorbar(fig[1:2, 4], hm, height = Relative(3/4), ticks = ([-4e-5, -2e-5, 0, 2e-5, 4e-5], [L"-4", L"-2", L"0", L"2", L"4"]), label = L"\text{Vertical vorticity } \left[s^{-1} \cdot 10^{-5} \right]",)

rowgap!(ga, 0)
rowgap!(gb, 0)

# CairoMakie.save("vorticity_cont.png", fig, px_per_unit = 8)
# CairoMakie.save("vorticity_cont.eps", fig, px_per_unit = 8)