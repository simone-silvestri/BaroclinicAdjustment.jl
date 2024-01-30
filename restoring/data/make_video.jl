using GLMakie
using JLD2

using Statistics
using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics

using Oceananigans
using Oceananigans.Grids: nodes

# using CairoMakie
# CairoMakie.activate!()

function enhance_vars(var::Matrix, rep = 1)
    var = repeat(var,  inner = (rep, 1))
    return repeat(var, inner = (1, rep))
end

enhance_vars(var::Vector, rep = 1) = repeat(var,  inner = rep)

# qg16 = all_fieldtimeseries("test.jld2", "./"; variables = ["u", "v", "w", "b", "zeta"])
qg16 = all_fieldtimeseries("qgleith_sixteen_video_level_30.jld2", "./"; variables = ["u", "v", "w", "b", "zeta"])
wn16 = all_fieldtimeseries("weno9pV_sixteen_video_level_30.jld2", "./"; variables = ["u", "v", "w", "b", "zeta"])

ζqg16 = qg16[:zeta]
ζwn16 = wn16[:zeta]

λ16, φ16, _ = enhance_vars.(Array.(nodes(ζqg16)))

λ16 = λ16[1:end-1]
φ16 = φ16[1:end-1]

iter = Observable(1)

ζqgl = @lift(Array(interior(Diagnostics.VerticalVorticity(qg16, $iter), :, :, 1)))
ζwnl = @lift(Array(interior(Diagnostics.VerticalVorticity(wn16, $iter), :, :, 1)))

g16 = LatitudeLongitudeGrid(size = (320, 320, 50), latitude = (-60, -40), longitude = (-10, 10), z = (-1000, 0))

lonlims16 = extrema(λ16)
latlims16 = extrema(φ16)

kwargs16 = (; )

ticks(a::Vector) = (a, [L"\text{%$b}" for b in a])

# possible cmaps:
# diverging_bkr_55_10_c35_n256
# diverging_linear_bjr_30_55_c53_n256

surfargs = (colorrange = (-4e-5, 4e-5), colormap = :bwr)

fig = Figure(fontsize=28, resolution = (930, 700) .* 1.5)

ga = fig[1, 1:3] = GridLayout()

ax = Axis(ga[1, 1]; title = L"\text{QG2, 7 kilometers}", 
                    ylabel = L"\text{Latitude}", 
                    yticks = ([-55, -50, -45], [L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}"]), 
                    xlabel = L"\text{Longitude}",
                    xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0^\text{o}", L"5^\text{o} \text{E}"]), kwargs16...)            
heatmap!(ax, λ16, φ16, ζqgl; surfargs...)
ax = Axis(ga[1, 2]; title = L"\text{W9V, 7 kilometers}", 
                    ylabel = L"\text{Latitude}", 
                    yticks = ([-55, -50, -45], [L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}"]), 
                    xlabel = L"\text{Longitude}",
                    xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0^\text{o}", L"5^\text{o} \text{E}"]), kwargs16...)            
hm = heatmap!(ax, λ16, φ16, ζwnl; surfargs...)

cb = Colorbar(fig[1, 4], hm, height = Relative(3/4), ticks = ([-4e-5, -2e-5, 0, 2e-5, 4e-5], [L"-4", L"-2", L"0", L"2", L"4"]), label = L"\text{Vertical vorticity } \left[s^{-1} \cdot 10^{-5} \right]",)

rowgap!(ga, 0)
# rowgap!(gb, 0)

Nt = length(qg16[:b].times)

record(fig, "vorticity.mp4", 1:2:Nt; framerate = 12) do i
    iter[] = i
    @info "doing iter $i"
end

surfargs = (colorrange = (-0.002, 0.0025), colormap = :thermal)

ζqgl = @lift(Array(interior(qg16[:b][$iter], :, :, 1)))
ζwnl = @lift(Array(interior(wn16[:b][$iter], :, :, 1)))

fig = Figure(fontsize=28, resolution = (930, 700) .* 1.5)

ga = fig[1, 1:3] = GridLayout()

ax = Axis(ga[1, 1]; title = L"\text{QG2, 7 kilometers}", 
                    ylabel = L"\text{Latitude}", 
                    yticks = ([-55, -50, -45], [L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}"]), 
                    xlabel = L"\text{Longitude}",
                    xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0^\text{o}", L"5^\text{o} \text{E}"]), kwargs16...)            
heatmap!(ax, λ16, φ16, ζqgl; surfargs...)
ax = Axis(ga[1, 2]; title = L"\text{W9V, 7 kilometers}", 
                    ylabel = L"\text{Latitude}", 
                    yticks = ([-55, -50, -45], [L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}"]), 
                    xlabel = L"\text{Longitude}",
                    xticks = ([-5, 0, 5], [L"5^\text{o} \text{W}", L"0^\text{o}", L"5^\text{o} \text{E}"]), kwargs16...)            
hm = heatmap!(ax, λ16, φ16, ζwnl; surfargs...)

cb = Colorbar(fig[1, 4], hm, height = Relative(3/4), ticks = ([-2e-3, 0, 2e-3], [L"0", L"1", L"2"]), label = L"\text{Temperature } [ ^\circ C ]",)

rowgap!(ga, 0)
# rowgap!(gb, 0)

record(fig, "buoyancy.mp4", 1:2:Nt; framerate = 12) do i
    iter[] = i
    @info "doing iter $i"
end

# CairoMakie.save("vorticity_cont.png", fig, px_per_unit = 8)
# CairoMakie.save("vorticity_cont.eps", fig, px_per_unit = 8)