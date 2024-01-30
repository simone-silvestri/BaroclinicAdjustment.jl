using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using GLMakie
using JLD2

using Statistics

colorqg = :deepskyblue2
colorwn = :firebrick2
colorwd = :orange1
# 
qg08 = jldopen("qgleith_eight_new_postprocess.jld2")
wn08 = jldopen("weno9pV_eight_new_postprocess.jld2")
qg16 = jldopen("qgleith_sixteen_new_postprocess.jld2")
wn16 = jldopen("weno9pV_twelve_new_new_postprocess.jld2")
qg32 = jldopen("qgleith_thirtytwo_new_postprocess.jld2")
wn32 = jldopen("weno9pV_thirtytwo_new_postprocess.jld2")

wd08 = jldopen("weno9pAllD_eight_new_postprocess.jld2")
wd16 = jldopen("weno9pAllD_sixteen_new_postprocess.jld2")

λ08 = range(-10,  10, length = 160)
φ08 = range(-60, -40, length = 160)
λ16 = range(-10,  10, length = 320)
φ16 = range(-60, -40, length = 320)
λ32 = range(-10,  10, length = 640)
φ32 = range(-60, -40, length = 640)
z  = range(-990, -10, length = 50)

vbw = mean(wn08["variance"].v′b′, dims = (1, 2))
vbp = mean(qg08["variance"].v′b′, dims = (1, 2))
vbq = mean(qg16["variance"].v′b′, dims = (1, 2))
vb6 = mean(wn16["variance"].v′b′, dims = (1, 2))
vbQ = mean(qg32["variance"].v′b′, dims = (1, 2))
vbW = mean(wn32["variance"].v′b′, dims = (1, 2))

vbd = mean(wd08["variance"].v′b′, dims = (1, 2))
vbD = mean(wd16["variance"].v′b′, dims = (1, 2))

GC.gc(true)

wbw = mean(wn08["variance"].w′b′, dims = (1, 2))
wbp = mean(qg08["variance"].w′b′, dims = (1, 2))
wbq = mean(qg16["variance"].w′b′, dims = (1, 2))
wb6 = mean(wn16["variance"].w′b′, dims = (1, 2))
wbQ = mean(qg32["variance"].w′b′, dims = (1, 2))
wbW = mean(wn32["variance"].w′b′, dims = (1, 2))

wbd = mean(wd08["variance"].w′b′, dims = (1, 2))
wbD = mean(wd16["variance"].w′b′, dims = (1, 2))

GC.gc(true)

uww = mean(wn08["variance"].u′w′, dims = (1, 2))
uwp = mean(qg08["variance"].u′w′, dims = (1, 2))
uwq = mean(qg16["variance"].u′w′, dims = (1, 2))
uw6 = mean(wn16["variance"].u′w′, dims = (1, 2))
uwQ = mean(qg32["variance"].u′w′, dims = (1, 2))
uwW = mean(wn32["variance"].u′w′, dims = (1, 2))

uwd = mean(wd08["variance"].u′w′, dims = (1, 2))
uwD = mean(wd16["variance"].u′w′, dims = (1, 2))

GC.gc(true)

vww = mean(wn08["variance"].v′w′, dims = (1, 2))
vwp = mean(qg08["variance"].v′w′, dims = (1, 2))
vwq = mean(qg16["variance"].v′w′, dims = (1, 2))
vw6 = mean(wn16["variance"].v′w′, dims = (1, 2))
vwQ = mean(qg32["variance"].v′w′, dims = (1, 2))
vwW = mean(wn32["variance"].v′w′, dims = (1, 2))

vwd = mean(wd08["variance"].v′w′, dims = (1, 2))
vwD = mean(wd16["variance"].v′w′, dims = (1, 2))

GC.gc(true)

figb = Figure(resolution = (1200, 500) .* 0.8, fontsize = 20)

ga = figb[1, 1:6] = GridLayout()

supertitle = Label(ga[0, 1:2], L"\text{Eddy fluxes}", fontsize = 20)
supertitle = Label(ga[0, 3:4], L"\text{Eddy stresses}", fontsize = 20)

ax1 = Axis(ga[1, 1], 
          xlabel = L"\overline{v^\prime b^\prime}(z) \cdot 10^5",
          ylabel = L"\text{Depth km}",
          yticks = ([-1000, -750, -500, -250, 0], [L"1", L"0.75", L"0.5", L"0.25", L"0"]),
          xticks = ([-8e-5, -6e-5, -4e-5, -2e-5], [L"-8", L"-6", L"-4", L"-2"]))
lines!(ax1, interior(vbp, 1, 1, 1:50), z[1:50], linewidth = 2, color = colorqg)
lines!(ax1, interior(vbw, 1, 1, 1:50), z[1:50], linewidth = 2, color = colorwn)
lines!(ax1, interior(vbq, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax1, interior(vb6, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
lines!(ax1, interior(vbQ, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax1, interior(vbW, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
# lines!(ax1, interior(vbd, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
# lines!(ax1, interior(vbD, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
scatter!(ax1, interior(vbq, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorqg)
scatter!(ax1, interior(vb6, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorwn)
scatter!(ax1, interior(vbQ, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorqg)
scatter!(ax1, interior(vbW, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorwn)

ax2 = Axis(ga[1, 2], 
          xlabel = L"\overline{w^\prime b^\prime}(z) \cdot 10^8",
          yticks = ([-1000, -750, -500, -250, 0], ["", "", "", "", ""]),
          xticks = ([0, 2, 4], [L"0", L"2", L"4"]))
lines!(ax2, 1e8 .* interior(wbp, 1, 1, 1:50), z[1:50], linewidth = 2, color = colorqg, label = L"\text{QG2 14km}")
lines!(ax2, 1e8 .* interior(wbw, 1, 1, 1:50), z[1:50], linewidth = 2, color = colorwn, label = L"\text{W9V 14km}")
lines!(ax2, 1e8 .* interior(wbq, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax2, 1e8 .* interior(wb6, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
lines!(ax2, 1e8 .* interior(wbQ, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax2, 1e8 .* interior(wbW, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
# lines!(ax2, 1e8 .* interior(wbd, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
# lines!(ax2, 1e8 .* interior(wbD, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
scatter!(ax2, 1e8 .* interior(wbq, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorqg)
scatter!(ax2, 1e8 .* interior(wb6, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorwn)
scatter!(ax2, 1e8 .* interior(wbQ, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorqg)
scatter!(ax2, 1e8 .* interior(wbW, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorwn)
scatterlines!(ax2, 1e8 .* interior(wbq, 1, 1, 1:1), z[1:1], linewidth = 1, marker = :star5,  markersize = 10, label = L"\text{QG2 7km}",    color = colorqg)
scatterlines!(ax2, 1e8 .* interior(wb6, 1, 1, 1:1), z[1:1], linewidth = 1, marker = :star5,  markersize = 10, label = L"\text{W9V 7km}",   color = colorwn)
scatterlines!(ax2, 1e8 .* interior(wbQ, 1, 1, 1:1), z[1:1], linewidth = 1, marker = :circle, markersize = 8,  label = L"\text{QG2 3.5km}",  color = colorqg)
scatterlines!(ax2, 1e8 .* interior(wbW, 1, 1, 1:1), z[1:1], linewidth = 1, marker = :circle, markersize = 8,  label = L"\text{W9V 3.5km}", color = colorwn)

ax = Axis(ga[1, 3], 
          xlabel = L"\overline{u^\prime w^\prime}(z) \cdot 10^6",
          yticks = ([-1000, -750, -500, -250, 0], ["", "", "", "", ""]),
          xticks = ([-1e-6, -5e-7, 0], [L"-1", L"-0.5", L"0"]))

lines!(ax, interior(uwp, 1, 1, :), z[1:50], linewidth = 2, color = colorqg)
lines!(ax, interior(uww, 1, 1, :), z[1:50], linewidth = 2, color = colorwn)
lines!(ax, interior(uwq, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax, interior(uw6, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
lines!(ax, interior(uwQ, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax, interior(uwW, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
# lines!(ax, interior(uwd, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
# lines!(ax, interior(uwD, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
scatter!(ax, interior(uwq, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorqg)
scatter!(ax, interior(uw6, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorwn)
scatter!(ax, interior(uwQ, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorqg)
scatter!(ax, interior(uwW, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorwn)

ax = Axis(ga[1, 4],
          xlabel = L"\overline{v^\prime w^\prime}(z) \cdot 10^6",
          yticks = ([-1000, -750, -500, -250, 0], ["", "", "", "", ""]),
          xticks = ([-4e-6, -2e-6, 0], [L"-4", L"-2", L"0"]))

lines!(ax, interior(vwp, 1, 1, :), z[1:50], linewidth = 2, color = colorqg)
lines!(ax, interior(vww, 1, 1, :), z[1:50], linewidth = 2, color = colorwn)
lines!(ax, interior(vwq, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax, interior(vw6, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
lines!(ax, interior(vwQ, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorqg)
lines!(ax, interior(vwW, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwn)
# lines!(ax, interior(vwd, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
# lines!(ax, interior(vwD, 1, 1, 1:50), z[1:50], linewidth = 1, color = colorwd)
scatter!(ax, interior(vwq, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorqg)
scatter!(ax, interior(vw6, 1, 1, 1:3:50), z[1:3:50], marker = :star5,  markersize = 10, color = colorwn)
scatter!(ax, interior(vwQ, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorqg)
scatter!(ax, interior(vwW, 1, 1, 1:3:50), z[1:3:50], marker = :circle, markersize = 8,  color = colorwn)

colgap!(ga, 5)
rowgap!(ga, 5)
leg = Legend(figb[1, 7], ax2, labelsize = 16)


# grid = Bw.grid
# B₀   = Field{Nothing, Center, Center}(grid)
# # Parameters
# param = (; vb = 4e-6, Δb = 0.005, Lz = 1000, Lφ = 20, φinit = 60)
# set!(B₀, (y, z) -> BaroclinicAdjustment.bᵢ(1, y, z, param))



