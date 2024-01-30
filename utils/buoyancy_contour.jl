using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using GLMakie
using JLD2

using Statistics
using ColorSchemes

colorqg = :deepskyblue2
colorwn = :firebrick2
colorwd = :orange1

qg08 = jldopen("qgleith_eight_new_postprocess.jld2")
wn08 = jldopen("upwind_eight_postprocess.jld2")
qg16 = jldopen("weno9pAllD_eight_new_postprocess.jld2")
wn16 = jldopen("weno9pV_eight_new_postprocess.jld2")
qg32 = jldopen("qgleith_thirtytwo_new_postprocess.jld2")
wn32 = jldopen("weno9pV_thirtytwo_new_postprocess.jld2")

wd08 = jldopen("weno9pAllD_sixteen_new_postprocess.jld2")
wd16 = jldopen("bilap_sixteen_postprocess.jld2")

λ08 = range(-10,  10, length = 160)
φ08 = range(-60, -40, length = 160)
φ12 = range(-60, -40, length = 240)
λ16 = range(-10,  10, length = 160)
φ16 = range(-60, -40, length = 160)
λ32 = range(-10,  10, length = 640)
φ32 = range(-60, -40, length = 640)
z  = range(-990, -10, length = 50)

Bd = wd08["mean"].B
BD = wd16["mean"].B
Bw = wn08["mean"].B
Bp = qg08["mean"].B
Bq = qg16["mean"].B
B6 = wn16["mean"].B
BQ = qg32["mean"].B
BW = wn32["mean"].B
z  = range(-990, -10, length = 50)

col = deepcopy(ColorSchemes.PRGn_5.colors)

new_col = []
for i in 1:2
    push!(new_col, RGBf(col[i].r, col[i].g, col[i].b * 0.65))
end
for i in 3:5
    push!(new_col, col[i])
end

mycmap = ColorScheme(typeof(new_col[1]).(new_col)) 

fig = Figure(resolution = (1050, 550), fontsize = 20)
ga = fig[1:4, 5]
ax = Axis(fig[2:4, 1:4], xlabel = L"\text{Latitude}", ylabel = L"\text{Depth km}",
          xticks = ([-60, -55, -50, -45, -40], [L"60^\text{o} \text{S}", L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}", L"40^\text{o} \text{S}"]),
          yticks = ([-250, -500, -750], [L"0.25", L"0.5", L"0.75"]))
hm = contourf!(ax, φ32, z, interior(BQ, 1, :, :), levels = range(-0.004, 0.006, length = 9), colormap = mycmap)
contour!(ax, φ08, z, interior(Bp, 1, :, :), levels = range(-0.004, 0.006, length = 9), linewidth = 3.5, color = colorqg)
contour!(ax, φ08, z, interior(Bw, 1, :, :), levels = range(-0.004, 0.006, length = 9), linewidth = 3.5, color = colorwn)

# lines!(ax, [1], label = L"\text{QG2 14 km}",  linewidth = 3.5, color = colorqg)
# lines!(ax, [1], label = L"\text{W9V 14 km}",  linewidth = 3.5, color = colorwn)

cplot = contour!(ax, φ16, z, interior(Bq, 1, :, :), levels = range(-0.004, 0.006, length = 9), linewidth = 2, color = colorqg)
beginnings = Point2f[]
# First plot in contour is the line plot, first arguments are the points of the contour
segments = cplot.plots[2][1][]
for (i, p) in enumerate(segments[1:15:end])
    # the segments are separated by NaN, which signals that a new contour starts
    if !isnan(p)
        push!(beginnings, p)
    end
end
scatter!(ax, beginnings, marker = :star5, markersize=12, color=colorqg)
# scq = lines!(ax, beginnings[1:1], marker = :star5, linewidth = 1, label = L"\text{Explicit 7km}", markersize=12, color=colorqg)

cplot = contour!(ax, φ16, z, interior(B6, 1, :, :), levels = range(-0.004, 0.006, length = 9), linewidth = 2, color = colorwn)
beginnings = Point2f[]
# First plot in contour is the line plot, first arguments are the points of the contour
segments = cplot.plots[2][1][]
for (i, p) in enumerate(segments[1:15:end])
    # the segments are separated by NaN, which signals that a new contour starts
    if !isnan(p)
        push!(beginnings, p)
    end
end
scatter!(ax, beginnings, marker = :star5, markersize=12, color=colorwn)
# scw = lines!(ax, beginnings[1:1], marker = :star5, linewidth = 1, label = L"\text{Implicit 7km}",  markersize=12, color=colorwn)

# cb = Colorbar(fig[1, 1:4], hm, vertical = false, label = L"\text{Filled contour: reference at 3.5-kilometer resolution}",
#               ticks = ([-0.0025, 0, 0.0025, 0.005], [L"-0.0025", L"0", L"0.0025", L"0.005"]))

# xlims!(ax, extrema(φ16))
# ylims!(ax, extrema(z))

# ax2 = Axis(ga[1:3, 1], xaxisposition = :top,
#            ylabel = L"\text{Depth km}", 
#            yticks = ([0, -500, -1000], [L"0", L"0.5", L"1"]),
#            xticks = ([0, 0.0005, 0.001], [L"0", L"0.0005", L"0.001"]),
#            xlabel = L"\overline{b} - \overline{b}_{\text{Ref}}\text{ at 42.5}^\text{o}")

# # lines!(ax2, interior(Bd, 1, 140, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorwd, linewidth = 2)
# # lines!(ax2, interior(Bp, 1, 140, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorqg, linewidth = 2)
# # lines!(ax2, interior(Bw, 1, 210, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorwn, linewidth = 2)
# # lines!(ax2, interior(BD, 1, 280, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorwd, linewidth = 1)
# lines!(ax2, interior(Bq, 1, 280, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorqg, linewidth = 2)
# lines!(ax2, interior(B6, 1, 280, 1:50) .- interior(BQ, 1, 560, 1:50), z[1:50], color = colorwn, linewidth = 2)
# scatter!(ax2, interior(Bq, 1, 280, 1:4:50) .- interior(BQ, 1, 560, 1:4:50), z[1:4:50], marker = :star5, markersize=10, color=colorqg)
# scatter!(ax2, interior(B6, 1, 280, 1:4:50) .- interior(BQ, 1, 560, 1:4:50), z[1:4:50], marker = :star5, markersize=10, color=colorwn)

# leg = Legend(ga[4, 1], ax)
