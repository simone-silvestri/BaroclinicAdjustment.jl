using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using CairoMakie
CairoMakie.activate!()
using JLD2

using Statistics

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

function convert_to_array(f::Vector, averaging = 10)
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

function get_energy(filename, averaging)

    var = jldopen(filename)

    MAPE = convert_to_array(var["energies"].Etimeseries.MAPE, averaging)
    EAPE = convert_to_array(var["energies"].Etimeseries.EAPE, averaging)
     MKE = convert_to_array(var["energies"].Etimeseries.MEKE, averaging)
     EKE = convert_to_array(var["energies"].Etimeseries.EKE, averaging)
    
     
     TKE =  MKE .+  EKE
    TAPE = MAPE .+ EAPE

    return (; MAPE, EAPE, MKE, EKE, TKE, TAPE)
end

avg = 10
wn08 = get_energy("weno9pV_eight_new_postprocess.jld2", avg)
qg08 = get_energy("qgleith2_eight_postprocess.jld2", avg)
up08 = get_energy("upwind_eight_postprocess.jld2", avg)
om08 = get_energy("omp25_eight_postprocess.jld2", avg)
wd08 = get_energy("weno9pAllD_correct_eight_postprocess.jld2", avg)

qg16 = get_energy("qgleith2_sixteen_postprocess.jld2", avg)
up16 = get_energy("upwind_sixteen_postprocess.jld2", avg)
om16 = get_energy("omp25_sixteen_postprocess.jld2", avg)
wd16 = get_energy("weno9pAllD_correct_sixteen_postprocess.jld2", avg)
wn16 = get_energy("weno9pV_sixteen_new_postprocess.jld2", avg)

qg32 = get_energy("qgleith2_thirtytwo_postprocess.jld2", avg)
up32 = get_energy("upwind_thirtytwo_postprocess.jld2", avg)
om32 = get_energy("omp25_thirtytwo_postprocess.jld2", avg)
wd32 = get_energy("upwind_thirtytwo_postprocess.jld2", avg)
wn32 = get_energy("weno9pV_thirtytwo_new_postprocess.jld2", avg)

others_expl_08 = (qg08, om08)
others_expl_16 = (qg16, om16)
others_expl_32 = (qg32, om32)

others_impl_08 = (wd08, up08)
others_impl_16 = (wd16, up16)
others_impl_32 = (wd32, up32)

TKEglob08_exp = [sum(var.TKE) for var in others_expl_08]
TKEglob16_exp = [sum(var.TKE) for var in others_expl_16]
TKEglob32_exp = [sum(var.TKE) for var in others_expl_32]

TKEglob08_imp = [sum(var.TKE) for var in others_impl_08]
TKEglob16_imp = [sum(var.TKE) for var in others_impl_16]
TKEglob32_imp = [sum(var.TKE) for var in others_impl_32]

mn08_e = others_expl_08[argmin(TKEglob08_exp)]
mx08_e = others_expl_08[argmax(TKEglob08_exp)]
mn16_e = others_expl_16[argmin(TKEglob16_exp)]
mn16_e = others_expl_16[argmin(TKEglob16_exp)]
mx32_e = others_expl_32[argmax(TKEglob32_exp)]
mx32_e = others_expl_32[argmax(TKEglob32_exp)]

mn08_i = others_impl_08[argmin(TKEglob08_imp)]
mx08_i = others_impl_08[argmax(TKEglob08_imp)]
mn16_i = others_impl_16[argmin(TKEglob16_imp)]
mx16_i = others_impl_16[argmax(TKEglob16_imp)]
mn32_i = others_impl_32[argmin(TKEglob32_imp)]
mx32_i = others_impl_32[argmax(TKEglob32_imp)]

maxTKE08_e = Float64[maximum(var.TKE[i] for var in others_expl_08) for i in 1:200]
minTKE08_e = Float64[minimum(var.TKE[i] for var in others_expl_08) for i in 1:200]
maxTKE16_e = Float64[maximum(var.TKE[i] for var in others_expl_16) for i in 1:200]
minTKE16_e = Float64[minimum(var.TKE[i] for var in others_expl_16) for i in 1:200]
maxTKE32_e = Float64[maximum(var.TKE[i] for var in others_expl_32) for i in 1:200]
minTKE32_e = Float64[minimum(var.TKE[i] for var in others_expl_32) for i in 1:200]

maxEKE08_e = Float64[maximum(var.EKE[i] for var in others_expl_08) for i in 1:200]
minEKE08_e = Float64[minimum(var.EKE[i] for var in others_expl_08) for i in 1:200]
maxEKE16_e = Float64[maximum(var.EKE[i] for var in others_expl_16) for i in 1:200]
minEKE16_e = Float64[minimum(var.EKE[i] for var in others_expl_16) for i in 1:200]
maxEKE32_e = Float64[maximum(var.EKE[i] for var in others_expl_32) for i in 1:200]
minEKE32_e = Float64[minimum(var.EKE[i] for var in others_expl_32) for i in 1:200]


maxEAPE08_e = Float64[maximum(var.EAPE[i] for var in others_expl_08) for i in 1:200]
minEAPE08_e = Float64[minimum(var.EAPE[i] for var in others_expl_08) for i in 1:200]
maxEAPE16_e = Float64[maximum(var.EAPE[i] for var in others_expl_16) for i in 1:200]
minEAPE16_e = Float64[minimum(var.EAPE[i] for var in others_expl_16) for i in 1:200]
maxEAPE32_e = Float64[maximum(var.EAPE[i] for var in others_expl_32) for i in 1:200]
minEAPE32_e = Float64[minimum(var.EAPE[i] for var in others_expl_32) for i in 1:200]

maxTKE08_i = Float64[maximum(var.TKE[i] for var in others_impl_08) for i in 1:200]
minTKE08_i = Float64[minimum(var.TKE[i] for var in others_impl_08) for i in 1:200]
maxTKE16_i = Float64[maximum(var.TKE[i] for var in others_impl_16) for i in 1:200]
minTKE16_i = Float64[minimum(var.TKE[i] for var in others_impl_16) for i in 1:200]
maxTKE32_i = Float64[maximum(var.TKE[i] for var in others_impl_32) for i in 1:200]
minTKE32_i = Float64[minimum(var.TKE[i] for var in others_impl_32) for i in 1:200]

maxEKE08_i = Float64[maximum(var.EKE[i] for var in others_impl_08) for i in 1:200]
minEKE08_i = Float64[minimum(var.EKE[i] for var in others_impl_08) for i in 1:200]
maxEKE16_i = Float64[maximum(var.EKE[i] for var in others_impl_16) for i in 1:200]
maxEKE16_i = Float64[maximum(var.EKE[i] for var in others_impl_16) for i in 1:200]
minEKE32_i = Float64[minimum(var.EKE[i] for var in others_impl_32) for i in 1:200]
minEKE32_i = Float64[minimum(var.EKE[i] for var in others_impl_32) for i in 1:200]

maxEAPE08_i = Float64[maximum(var.EAPE[i] for var in others_impl_08) for i in 1:200]
minEAPE08_i = Float64[minimum(var.EAPE[i] for var in others_impl_08) for i in 1:200]
maxEAPE16_i = Float64[maximum(var.EAPE[i] for var in others_impl_16) for i in 1:200]
minEAPE16_i = Float64[minimum(var.EAPE[i] for var in others_impl_16) for i in 1:200]
maxEAPE32_i = Float64[maximum(var.EAPE[i] for var in others_impl_32) for i in 1:200]
minEAPE32_i = Float64[minimum(var.EAPE[i] for var in others_impl_32) for i in 1:200]

color08 = Makie.color("blueviolet")
color16 = Makie.color("steelblue")   
color32 = Makie.color("lightskyblue") 

fig = Figure(resolution = (900, 400), fontsize = 15)
ga = fig[1:2, 1:6] = GridLayout()

ax = Axis(ga[1, 1:2], title = L"\text{Total KE}", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, minTKE08_e, maxTKE08_e, color = (color08, 0.4))
band!(ax, 1:200, minTKE16_e, maxTKE16_e, color = (color16, 0.4))
band!(ax, 1:200, minTKE32_e, maxTKE32_e, color = (color32, 0.4))
lines!(ax, mn08_e.TKE, linewidth = 1,    color = color08)
lines!(ax, mx08_e.TKE, linewidth = 1,    color = color08)
lines!(ax, mn16_e.TKE, linewidth = 1,    color = color16)
lines!(ax, mx16_e.TKE, linewidth = 1,    color = color16)
lines!(ax, mn32_e.TKE, linewidth = 1,    color = color32)
lines!(ax, mx32_e.TKE, linewidth = 1,    color = color32)

lines!(ax, wn08.TKE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.TKE, linewidth = 2,    color = color16)
lines!(ax, wn32.TKE, linewidth = 2,    color = color32)

ax = Axis(ga[1, 3:4], title = L"\text{Eddy KE}", 
          xlabel = "", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, minEKE08_e, maxEKE08_e, color = (color08, 0.4))
band!(ax, 1:200, minEKE16_e, maxEKE16_e, color = (color16, 0.4))
band!(ax, 1:200, minEKE32_e, maxEKE32_e, color = (color32, 0.4))
lines!(ax, mn08_e.EKE, linewidth = 1,    color = color08)
lines!(ax, mx08_e.EKE, linewidth = 1,    color = color08)
lines!(ax, mn16_e.EKE, linewidth = 1,    color = color16)
lines!(ax, mx16_e.EKE, linewidth = 1,    color = color16)
lines!(ax, mn32_e.EKE, linewidth = 1,    color = color16)
lines!(ax, mx32_e.EKE, linewidth = 1,    color = color16)

# lines!(ax, eb08.EKE, linewidth = 1,    color = color08, linestyle = :dash)
lines!(ax, wn08.EKE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.EKE, linewidth = 2,    color = color16)
lines!(ax, wn32.EKE, linewidth = 2,    color = color32)
lines!(ax, wn208.EKE, linewidth = 0.2, color = color08)
lines!(ax, wn216.EKE, linewidth = 0.2, color = color16)
lines!(ax, wn232.EKE, linewidth = 0.2, color = color32)

colgap!(ga, 5)
rowgap!(ga, 5)

ax = Axis(ga[1, 5:6], title = L"\text{Eddy APE}", 
          xlabel = "", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.03, 0.06], [L"0.0", L"0.03", L"0.06"]))

band!(ax, 1:200, minEAPE08_e, maxEAPE08_e, color = (color08, 0.4))
band!(ax, 1:200, minEAPE16_e, maxEAPE16_e, color = (color16, 0.4))
band!(ax, 1:200, minEAPE32_e, maxEAPE32_e, color = (color16, 0.4))
lines!(ax, mn08_e.EAPE, linewidth = 1,     color = color08)
lines!(ax, mx08_e.EAPE, linewidth = 1,     color = color08)
lines!(ax, mn16_e.EAPE, linewidth = 1,     color = color16)
lines!(ax, mx16_e.EAPE, linewidth = 1,     color = color16)
lines!(ax, mn32_e.EAPE, linewidth = 1,     color = color32)
lines!(ax, mx32_e.EAPE, linewidth = 1,     color = color32)

lines!(ax, wn08.EAPE, linewidth = 2, color = color08, label = "let's go")
lines!(ax, wn16.EAPE, linewidth = 2, color = color16)
lines!(ax, wn32.EAPE, linewidth = 2, color = color32)

colgap!(ga, 5)
rowgap!(ga, 5)

ax = Axis(ga[2, 1:2], # title = L"\text{Total KE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, minTKE08_i, maxTKE08_i, color = (color08, 0.4))
band!(ax, 1:200, minTKE16_i, maxTKE16_i, color = (color16, 0.4))
band!(ax, 1:200, minTKE32_i, maxTKE32_i, color = (color16, 0.4))
lines!(ax, mn08_i.TKE, linewidth = 1,    color = color08)
lines!(ax, mx08_i.TKE, linewidth = 1,    color = color08)
lines!(ax, mn16_i.TKE, linewidth = 1,    color = color16)
lines!(ax, mx16_i.TKE, linewidth = 1,    color = color16)
lines!(ax, mn32_i.TKE, linewidth = 1,    color = color32)
lines!(ax, mx32_i.TKE, linewidth = 1,    color = color32)

lines!(ax, wn08.TKE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.TKE, linewidth = 2,    color = color16)
lines!(ax, wn32.TKE, linewidth = 2,    color = color32)

ax = Axis(ga[2, 3:4], # title = L"\text{Eddy KE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, minEKE08_i, maxEKE08_i, color = (color08, 0.4))
band!(ax, 1:200, minEKE16_i, maxEKE16_i, color = (color16, 0.4))
band!(ax, 1:200, minEKE32_i, maxEKE32_i, color = (color32, 0.4))
lines!(ax, mn08_i.EKE, linewidth = 1,    color = color08)
lines!(ax, mx08_i.EKE, linewidth = 1,    color = color08)
lines!(ax, mn16_i.EKE, linewidth = 1,    color = color16)
lines!(ax, mx16_i.EKE, linewidth = 1,    color = color16)
lines!(ax, mn32_i.EKE, linewidth = 1,    color = color32)
lines!(ax, mx32_i.EKE, linewidth = 1,    color = color32)

lines!(ax, wn08.EKE, linewidth = 2,      color = color08, label = "let's go")
lines!(ax, wn16.EKE, linewidth = 2,      color = color16)
lines!(ax, wn32.EKE, linewidth = 2,      color = color32)

ax = Axis(ga[2, 5:6], # title = L"\text{Eddy APE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.03, 0.06], [L"0.0", L"0.03", L"0.06"]))

band!(ax, 1:200, minEAPE08_i, maxEAPE08_i, color = (color08, 0.4))
band!(ax, 1:200, minEAPE16_i, maxEAPE16_i, color = (color16, 0.4))
band!(ax, 1:200, minEAPE32_i, maxEAPE32_i, color = (color32, 0.4))
lines!(ax, mn08_i.EAPE, linewidth = 1,     color = color08)
lines!(ax, mx08_i.EAPE, linewidth = 1,     color = color08)
lines!(ax, mn16_i.EAPE, linewidth = 1,     color = color16)
lines!(ax, mx16_i.EAPE, linewidth = 1,     color = color16)
lines!(ax, mn32_i.EAPE, linewidth = 1,     color = color32)
lines!(ax, mx32_i.EAPE, linewidth = 1,     color = color32)

lines!(ax, wn08.EAPE, linewidth = 2,       color = color08, label = "let's go")
lines!(ax, wn16.EAPE, linewidth = 2,       color = color16)
lines!(ax, wn32.EAPE, linewidth = 2,       color = color32)

colgap!(ga, 5)
rowgap!(ga, 5)

CairoMakie.save("new_energy_plots.eps", fig)
 