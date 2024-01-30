using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using GLMakie
using JLD2

const colorQG = :deepskyblue
const colorWN = :firebrick2
const colorWD = :orange1

using Statistics

wn08 = jldopen("weno9pV_eight_new_postprocess.jld2")
qg08 = jldopen("qgleith_eight_new_postprocess.jld2")
wd08 = jldopen("weno9pAllD_eight_new_postprocess.jld2")
qg16 = jldopen("upwind_sixteen_postprocess.jld2")
wd16 = jldopen("qgleith_sixteen_new_postprocess.jld2")
wn16 = jldopen("weno9pV_sixteen_new_postprocess.jld2")
qg32 = jldopen("qgleith_thirtytwo_new_postprocess.jld2")
wn32 = jldopen("weno9pV_thirtytwo_new_postprocess.jld2")

λ8 = range(-10,   10, length = 160)
φ8 = range(-60,  -40, length = 160)
λ6 = range(-10,   10, length = 320)
φ6 = range(-60,  -40, length = 320)
z  = range(-990, -10, length = 50)

function convert_to_array(f::FieldTimeSeries, averaging = 10)
    arr = zeros(length(f.times))
    for t in 1:averaging ÷ 2 + 1
        arr[t] = f[t][1, 1, 1]
    end
    for t in averaging ÷ 2 + 1:length(f.times) - averaging ÷ 2 - 1
        arr[t] = sum(Tuple(f[i][1, 1, 1] for i in t - averaging ÷ 2:t + averaging ÷ 2)) / (averaging + 1)
    end
    for t in length(f.times) - averaging ÷ 2:length(f.times)
        arr[t] = f[t][1, 1, 1]
    end
    return arr
end

function convert_to_array(f::Vector, averaging = 40)
    arr = zeros(length(f))
    f = Float64.(f)
    for t in 1:averaging ÷ 2 + 1
        arr[t] = sum(Tuple(f[i] for i in 1:2t-1)) / (2t-1)
    end
    for t in averaging ÷ 2 + 1:length(f) - averaging ÷ 2 - 1
        arr[t] = sum(Tuple(f[i] for i in t - averaging ÷ 2:t + averaging ÷ 2)) / (averaging + 1)
    end
    for t in length(f) - averaging ÷ 2:length(f)
        arr[t] = sum(Tuple(f[i] for i in 2t-length(f):length(f))) / (2 * (length(f) - t) + 1)
    end
    return arr
end

MAPEw = convert_to_array(wn08["energies"].Etimeseries.MAPE)
EAPEw = convert_to_array(wn08["energies"].Etimeseries.EAPE)
 MKEw = convert_to_array(wn08["energies"].Etimeseries.MEKE)
 EKEw = convert_to_array(wn08["energies"].Etimeseries.EKE)

MAPEq = convert_to_array(qg16["energies"].Etimeseries.MAPE)
EAPEq = convert_to_array(qg16["energies"].Etimeseries.EAPE)
 MKEq = convert_to_array(qg16["energies"].Etimeseries.MEKE)
 EKEq = convert_to_array(qg16["energies"].Etimeseries.EKE)

MAPEp = convert_to_array(qg08["energies"].Etimeseries.MAPE)
EAPEp = convert_to_array(qg08["energies"].Etimeseries.EAPE)
 MKEp = convert_to_array(qg08["energies"].Etimeseries.MEKE)
 EKEp = convert_to_array(qg08["energies"].Etimeseries.EKE)

MAPE6 = convert_to_array(wn16["energies"].Etimeseries.MAPE)
EAPE6 = convert_to_array(wn16["energies"].Etimeseries.EAPE)
 MKE6 = convert_to_array(wn16["energies"].Etimeseries.MEKE)
 EKE6 = convert_to_array(wn16["energies"].Etimeseries.EKE)

MAPEQ = convert_to_array(qg32["energies"].Etimeseries.MAPE)
EAPEQ = convert_to_array(qg32["energies"].Etimeseries.EAPE)
 MKEQ = convert_to_array(qg32["energies"].Etimeseries.MEKE)
 EKEQ = convert_to_array(qg32["energies"].Etimeseries.EKE)

 MAPEW = convert_to_array(wn32["energies"].Etimeseries.MAPE)
 EAPEW = convert_to_array(wn32["energies"].Etimeseries.EAPE)
  MKEW = convert_to_array(wn32["energies"].Etimeseries.MEKE)
  EKEW = convert_to_array(wn32["energies"].Etimeseries.EKE)
 
MAPEd = convert_to_array(wd08["energies"].Etimeseries.MAPE)
EAPEd = convert_to_array(wd08["energies"].Etimeseries.EAPE)
 MKEd = convert_to_array(wd08["energies"].Etimeseries.MEKE)
 EKEd = convert_to_array(wd08["energies"].Etimeseries.EKE)

MAPED = convert_to_array(wd16["energies"].Etimeseries.MAPE)
EAPED = convert_to_array(wd16["energies"].Etimeseries.EAPE)
 MKED = convert_to_array(wd16["energies"].Etimeseries.MEKE)
 EKED = convert_to_array(wd16["energies"].Etimeseries.EKE)

totq = MAPEq .+ EAPEq .+ MKEq .+ EKEq
totw = MAPEw .+ EAPEw .+ MKEw .+ EKEw
totp = MAPEp .+ EAPEp .+ MKEp .+ EKEp
# tot6 = MAPE6 .+ EAPE6 .+ MKE6 .+ EKE6

fig = Figure(resolution = (900, 400), fontsize = 10)
ga = fig[1:2, 1:4] = GridLayout()


ax = Axis(ga[2, 3:4], ylabel = L"\text{Total APE}", 
          xlabel = L"\text{days}",
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0.3, 0.5, 0.7], [L"0.3", L"0.5", L"0.7"]))

lines!(ax, MAPEp .+ EAPEp, linewidth = 2, color = colorQG, label = L"\text{QG2 14km}")
lines!(ax, MAPEw .+ EAPEw, linewidth = 2, color = colorWN,  label = L"\text{W9V 14km}")
lines!(ax, MAPEd .+ EAPEd, linewidth = 1, color = colorWD)
lines!(ax, MAPED .+ EAPED, linewidth = 1, color = colorWD)
lines!(ax, MAPEq .+ EAPEq, linewidth = 1, color = colorQG)
lines!(ax, MAPEQ .+ EAPEQ, linewidth = 1, color = colorQG)
lines!(ax, MAPE6 .+ EAPE6, linewidth = 1, color = colorWN)
lines!(ax, MAPEW .+ EAPEW, linewidth = 1, color = colorWN)
scatter!(ax, (1:10:200), MAPEq[1:10:end] .+ EAPEq[1:10:end], marker = :star5,  markersize = 10, color = colorQG)
scatter!(ax, (1:10:200), MAPE6[1:10:end] .+ EAPE6[1:10:end], marker = :star5,  markersize = 10, color = colorWN)
scatter!(ax, (1:10:200), MAPEQ[1:10:end] .+ EAPEQ[1:10:end], marker = :circle, markersize = 8,  color = colorQG)
scatter!(ax, (1:10:200), MAPEW[1:10:end] .+ EAPEW[1:10:end], marker = :circle, markersize = 8,  color = colorWN)
scatterlines!(ax, MAPEq[1:1] .+ EAPEq[1:1], linewidth = 1, marker = :star5,  markersize = 10, label = L"\text{Explicit 7km}",    color = colorQG)
scatterlines!(ax, MAPE6[1:1] .+ EAPE6[1:1], linewidth = 1, marker = :star5,  markersize = 10, label = L"\text{Implicit 7km}",   color = colorWN)
scatterlines!(ax, MAPEQ[1:1] .+ EAPEQ[1:1], linewidth = 1, marker = :circle, markersize = 8,  label = L"\text{QG2 3.5km}",  color = colorQG)
scatterlines!(ax, MAPEW[1:1] .+ EAPEW[1:1], linewidth = 1, marker = :circle, markersize = 8,  label = L"\text{W9V 3.5km}", color = colorWN)

leg = Legend(fig[1:2, 5], ax)

ax = Axis(ga[2, 1:2], ylabel = L"\text{Total KE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

lines!(ax, MKEp .+ EKEp, linewidth = 2, color = colorQG)
lines!(ax, MKEw .+ EKEw, linewidth = 2, color = colorWN)
lines!(ax, MKEd .+ EKEd, linewidth = 1, color = colorWD)
lines!(ax, MKED .+ EKED, linewidth = 1, color = colorWD)
lines!(ax, MKEq .+ EKEq, linewidth = 1, color = colorQG)
lines!(ax, MKEQ .+ EKEQ, linewidth = 1, color = colorQG)
lines!(ax, MKE6 .+ EKE6, linewidth = 1, color = colorWN)
lines!(ax, MKEW .+ EKEW, linewidth = 1, color = colorWN)
scatter!(ax, (1:7:200), MKEq[1:7:end] .+ EKEq[1:7:end], marker = :star5,  markersize = 8, color = colorQG)
scatter!(ax, (1:7:200), MKE6[1:7:end] .+ EKE6[1:7:end], marker = :star5,  markersize = 8, color = colorWN)
scatter!(ax, (1:10:200), MKEQ[1:10:end] .+ EKEQ[1:10:end], marker = :circle, markersize = 8,  color = colorQG)
scatter!(ax, (1:10:200), MKEW[1:10:end] .+ EKEW[1:10:end], marker = :circle, markersize = 8,  color = colorWN)

ax = Axis(ga[1, 3:4], ylabel = L"\text{Eddy APE}", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.03, 0.06], [L"0.0", L"0.03", L"0.06"]))

lines!(ax, EAPEp, linewidth = 2, color = colorQG)
lines!(ax, EAPEw, linewidth = 2, color = colorWN)
lines!(ax, EAPEd, linewidth = 1, color = colorWD)
lines!(ax, EAPED, linewidth = 1, color = colorWD)
lines!(ax, EAPEq, linewidth = 1, color = colorQG)
lines!(ax, EAPE6, linewidth = 1, color = colorWN)
lines!(ax, EAPEQ, linewidth = 1, color = colorQG)
lines!(ax, EAPEW, linewidth = 1, color = colorWN)
scatter!(ax, (1:10:200), EAPEq[1:10:end], marker = :star5,  markersize = 10, color = colorQG)
scatter!(ax, (1:10:200), EAPE6[1:10:end], marker = :star5,  markersize = 10, color = colorWN)
scatter!(ax, (1:10:200), EAPEQ[1:10:end], marker = :circle, markersize = 8,  color = colorQG)
scatter!(ax, (1:10:200), EAPEW[1:10:end], marker = :circle, markersize = 8,  color = colorWN)

ax = Axis(ga[1, 1:2], ylabel = L"\text{Eddy KE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

lines!(ax, EKEp, linewidth = 2, color = colorQG, label = L"\text{QG2 14 km}")
lines!(ax, EKEw, linewidth = 2, color = colorWN, label = L"\text{WN2 14 km}")
lines!(ax, EKEd, linewidth = 1, color = colorWD)
lines!(ax, EKED, linewidth = 1, color = colorWD)
lines!(ax, EKEq, linewidth = 1, color = colorQG)
lines!(ax, EKE6, linewidth = 1, color = colorWN)
lines!(ax, EKEQ, linewidth = 1, color = colorQG)
lines!(ax, EKEW, linewidth = 1, color = colorWN)
scatter!(ax, (1:10:200), EKEq[1:10:end], marker = :star5,  markersize = 10, color = colorQG)
scatter!(ax, (1:10:200), EKE6[1:10:end], marker = :star5,  markersize = 10, color = colorWN)
scatter!(ax, (1:10:200), EKEQ[1:10:end], marker = :circle, markersize = 8,  color = colorQG)
scatter!(ax, (1:10:200), EKEW[1:10:end], marker = :circle, markersize = 8,  color = colorWN)

N²w = wn08["stratif"]
N²p = qg08["stratif"]
N²q = qg16["stratif"]
N²6 = wn16["stratif"]
N²Q = qg32["stratif"]
N²W = wn32["stratif"]

# ax = Axis(ga[2, 5:6], ylabel = L"\text{Average }N^2", 
#           xlabel = L"\text{days}", 
#           xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
#           yticks = ([0., 0.1, 0.2, 0.3], [L"0.0", L"0.1", L"0.2", L"0.3"]))

# lines!(ax, N²p, linewidth = 2, color = colorQG, label = L"\text{QG 14 km}")
# lines!(ax, N²w, linewidth = 2, color = colorWN, label = L"\text{WN 14 km}")
# lines!(ax, N²q, linewidth = 2, color = colorQG, label = L"\text{QG 7 km}", linestyle = :dash)
# lines!(ax, N²6, linewidth = 2, color = colorWN, label = L"\text{WN 7 km}", linestyle = :dash)
# lines!(ax, N²Q, linewidth = 2, color = colorQG, label = L"\text{QG 3.5 km}",linestyle = :dashdot)
# lines!(ax, N²W, linewidth = 2, color = colorWN, label = L"\text{WN 3.5 km}",linestyle = :dashdot)

colgap!(ga, 5)
rowgap!(ga, 5)

# using CairoMakie
# CairoMakie.activate!()
 