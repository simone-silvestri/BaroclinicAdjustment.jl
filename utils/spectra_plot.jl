using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using GLMakie
using JLD2

using Statistics

const colorQG = :deepskyblue
const colorWN = :firebrick2

jrange = 40:50

qg08 = jldopen("qgleith_eight_new_postprocess.jld2")
wn08 = jldopen("weno9pV_eight_new_postprocess.jld2")
qg16 = jldopen("qgleith_sixteen_new_postprocess.jld2")
wn16 = jldopen("weno9pV_sixteen_new_postprocess.jld2")
qg32 = jldopen("qgleith_thirtytwo_new_postprocess.jld2")
wn32 = jldopen("weno9pV_thirtytwo_new_postprocess.jld2")

λ8 = range(-10,   10, length = 160)
φ8 = range(-60,  -40, length = 160)
λ6 = range(-10,   10, length = 320)
φ6 = range(-60,  -40, length = 320)
z  = range(-990, -10, length = 50)

Ew = abs.(real.(wn08["spectra"].V .+ wn08["spectra"].U))
Eq = abs.(real.(qg16["spectra"].V .+ qg16["spectra"].U))
Ep = abs.(real.(qg08["spectra"].V .+ qg08["spectra"].U))
E6 = abs.(real.(wn16["spectra"].V .+ wn16["spectra"].U))
EQ = abs.(real.(qg32["spectra"].V .+ qg32["spectra"].U))
EW = abs.(real.(wn32["spectra"].V .+ wn32["spectra"].U))

WBw = abs.(real.(wn08["spectra"].WB))
WBq = abs.(real.(qg16["spectra"].WB))
WBp = abs.(real.(qg08["spectra"].WB))
WB6 = abs.(real.(wn16["spectra"].WB))
WBQ = abs.(real.(qg32["spectra"].WB))
WBW = abs.(real.(wn32["spectra"].WB))

Bw = abs.(real.(wn08["spectra"].Ω))
Bq = abs.(real.(qg16["spectra"].Ω))
Bp = abs.(real.(qg08["spectra"].Ω))
B6 = abs.(real.(wn16["spectra"].Ω))
BQ = abs.(real.(qg32["spectra"].Ω))
BW = abs.(real.(wn32["spectra"].Ω))

figs = Figure(resolution = (700, 300) .* 1.5)
ax = Axis(figs[1:2, 1], ylabel = L"E(k)", 
          xscale = log10, yscale = log10, 
          xticks = ([1e-5, 1e-4, 1e-3], ["", "", ""]),
          yticks = ([1e-4, 1e-2, 1], [L"10^{-4}", L"10^{-2}", L"10^0"]))
freqw = Ew[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 
freqq = Eq[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 
freqp = Ep[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 
freq6 = E6[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 
freqQ = EQ[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 
freqW = EW[1, 1].freq[2:end] .* 1/8 ./ 14kilometers 

lines!(ax, freqw, mean(Ew[:, jrange]).spec[2:end], linewidth = 2, color = colorWN)
lines!(ax, freqp, mean(Ep[:, jrange]).spec[2:end], linewidth = 2, color = colorQG)
lines!(ax, freqq, mean(Eq[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorQG)
lines!(ax, freq6, mean(E6[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorWN)
lines!(ax, freqQ, mean(EQ[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dashdot, color = colorQG)
lines!(ax, freqQ, mean(EW[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dashdot, color = colorWN)

ax = Axis(figs[3, 1], ylabel = L"\text{W9V/QG}", xlabel = L"\text{Wavenumber m}^{-1}", xscale = log10, yticks = ([1, 1.5, 2, 2.5], [L"1.0", L"1.5", L"2.0", L"2.5"]), xticks = ([1e-5, 1e-4, 1e-3], [L"10^{-5}", L"10^{-4}", L"10^{-3}"]))
lines!(ax, freqw, mean(Ew[:, jrange]).spec[2:end] ./ mean(Ep[:, jrange]).spec[2:end], linewidth = 2, color = :gray)
lines!(ax, freq6, mean(E6[:, jrange]).spec[2:end] ./ mean(Eq[:, jrange]).spec[2:end], linewidth = 2, color = :gray, linestyle = :dash)
lines!(ax, freqQ, mean(EW[:, jrange]).spec[2:end] ./ mean(EQ[:, jrange]).spec[2:end], linewidth = 2, color = :gray, linestyle = :dashdot)

vlines!(ax, [2π / 8000], xmin = 0.3, xmax = 2.7)

ax = Axis(figs[1:2, 2], ylabel = L"\mathcal{G} (k)", 
          xscale = log10, yscale = log10, 
          xticks = ([1e-5, 1e-4, 1e-3], ["", "", ""]),
          yticks = ([1e-11, 1e-10, 1e-9], [L"10^{-11}", L"10^{-10}", L"10^{-9}"]))

lines!(ax, freqp, mean(Bp[:, jrange]).spec[2:end], label = L"\text{QG, 8 km}", linewidth = 2, color = colorQG)
lines!(ax, freqw, mean(Bw[:, jrange]).spec[2:end], label = L"\text{W9V, 8 km}",  linewidth = 2, color = colorWN)
lines!(ax, freqq, mean(Bq[:, jrange]).spec[2:end], label = L"\text{QG, 4 km}",  linewidth = 2, linestyle = :dash, color = colorQG)
lines!(ax, freq6, mean(B6[:, jrange]).spec[2:end], label = L"\text{W9V, 4 km}", linewidth = 2, linestyle = :dash, color = colorWN)
lines!(ax, freqQ, mean(BQ[:, jrange]).spec[2:end], label = L"\text{W9V, 4 km}", linewidth = 2, linestyle = :dashdot, color = colorQG)
lines!(ax, freqW, mean(BW[:, jrange]).spec[2:end], label = L"\text{W9V, 4 km}", linewidth = 2, linestyle = :dashdot, color = colorWN)

axislegend(ax, position = :lb)

ax = Axis(figs[3, 2], ylabel = L"\text{W9V/QG}", xlabel = L"\text{Wavenumber m}^{-1}", xscale = log10, yticks = ([1, 1.5, 2, 2.5], [L"1.0", L"1.5", L"2.0", L"2.5"]), xticks = ([1e-5, 1e-4, 1e-3], [L"10^{-5}", L"10^{-4}", L"10^{-3}"]))
lines!(ax, freqw, mean(Bw[:, jrange]).spec[2:end] ./ mean(Bp[:, jrange]).spec[2:end], label = L"\text{8 km}", linewidth = 2, color = :gray)
lines!(ax, freq6, mean(B6[:, jrange]).spec[2:end] ./ mean(Bq[:, jrange]).spec[2:end], label = L"\text{4 km}", linewidth = 2, color = :gray, linestyle = :dash)
lines!(ax, freqQ, mean(BW[:, jrange]).spec[2:end] ./ mean(BQ[:, jrange]).spec[2:end], label = L"\text{4 km}", linewidth = 2, color = :gray, linestyle = :dashdot)

vlines!(ax, [2π / 8000], xmin = 0.3, xmax = 2.7)

axislegend(ax, position = :ct)

ax = Axis(figs[1:2, 3], ylabel = L"\overline{w^\prime b^\prime} (k)", 
          xscale = log10, yscale = log10, 
          xticks = ([1e-5, 1e-4, 1e-3], ["", "", ""]),
          yticks = ([1e-10, 1e-9, 1e-8, 1e-7], [L"10^{-10}", L"10^{-9}", L"10^{-8}", L"10^{-7}"]))

lines!(ax, freqw, mean(WBw[:, jrange]).spec[2:end], linewidth = 2, color = colorWN)
lines!(ax, freqp, mean(WBp[:, jrange]).spec[2:end], linewidth = 2, color = colorQG)
lines!(ax, freqq, mean(WBq[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorQG)
lines!(ax, freq6, mean(WB6[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorWN)
lines!(ax, freqQ, mean(WBQ[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorWN)
lines!(ax, freqW, mean(WBW[:, jrange]).spec[2:end], linewidth = 2, linestyle = :dash, color = colorWN)

ax = Axis(figs[3, 3], ylabel = L"\text{W9V/QG}", xlabel = L"\text{Wavenumber m}^{-1}", xscale = log10, yticks = ([0.5, 1, 1.5], [L"0.5", L"1.0", L"1.5"]), xticks = ([1e-5, 1e-4, 1e-3], [L"10^{-5}", L"10^{-4}", L"10^{-3}"]))
lines!(ax, freqw, mean(WBw[:, jrange]).spec[2:end] ./ mean(WBp[:, jrange]).spec[2:end], linewidth = 2, color = :gray)
lines!(ax, freq6, mean(WB6[:, jrange]).spec[2:end] ./ mean(WBq[:, jrange]).spec[2:end], linewidth = 2, color = :gray, linestyle = :dash)
lines!(ax, freqQ, mean(WBW[:, jrange]).spec[2:end] ./ mean(WBQ[:, jrange]).spec[2:end], linewidth = 2, color = :gray, linestyle = :dashdot)

vlines!(ax, [2π / 8000], xmin = 0.3, xmax = 2.7)

rowgap!(figs.layout, 0)
colgap!(figs.layout, 10)

# CairoMakie.save("spectra_new.eps", figs)