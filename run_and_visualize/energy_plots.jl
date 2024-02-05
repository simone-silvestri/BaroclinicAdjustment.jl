using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using Oceananigans
using Oceananigans.Units

using GLMakie
using JLD2
using ColorTypes
using FixedPointNumbers

import Base

Base.:+(c::ColorTypes.RGB{FixedPointNumbers.N0f8}, f::Float64) = ColorTypes.RGB{FixedPointNumbers.N0f8}(min(1,c.b + f), min(1,c.r + f), min(1,c.g + f))
Base.:+(c::ColorTypes.RGB{FixedPointNumbers.N0f8}, f::Tuple) = ColorTypes.RGB{FixedPointNumbers.N0f8}(min(1,c.b + f[1]), min(1,c.r + f[2]), min(1,c.g + f[3]))

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

avg = 50
wn08 = get_energy("weno9pV_eight_new_postprocess.jld2", avg)
qg08 = get_energy("qgleith2_eight_postprocess.jld2", avg)
up08 = get_energy("upwind_eight_postprocess.jld2", avg)
eb08 = get_energy("ebs_eight_postprocess.jld2", avg)
om08 = get_energy("omp25_eight_postprocess.jld2", avg)
wd08 = get_energy("weno9pAllD_eight_new_postprocess.jld2", avg)
qg16 = get_energy("qgleith2_sixteen_postprocess.jld2", avg)
up16 = get_energy("upwind_sixteen_postprocess.jld2", avg)
om16 = get_energy("omp25_sixteen_postprocess.jld2", avg)
wd16 = get_energy("weno9pAllD_sixteen_new_postprocess.jld2", avg)
wn16 = get_energy("weno9pV_sixteen_new_postprocess.jld2", avg)
qg32 = get_energy("qgleith_thirtytwo_new_postprocess.jld2", avg)
wn32 = get_energy("weno9pV_thirtytwo_new_postprocess.jld2", avg)

avg2 = 2
wn208 = get_energy("weno9pV_eight_new_postprocess.jld2", avg2)
qg208 = get_energy("qgleith2_eight_postprocess.jld2", avg2)
up208 = get_energy("upwind_eight_postprocess.jld2", avg2)
eb208 = get_energy("ebs_eight_postprocess.jld2", avg2)
om208 = get_energy("omp25_eight_postprocess.jld2", avg2)
wd208 = get_energy("weno9pAllD_eight_new_postprocess.jld2", avg2)
qg216 = get_energy("qgleith2_sixteen_postprocess.jld2", avg2)
up216 = get_energy("upwind_sixteen_postprocess.jld2", avg)
om216 = get_energy("omp25_sixteen_postprocess.jld2", avg2)
wd216 = get_energy("weno9pAllD_sixteen_new_postprocess.jld2", avg2)
wn216 = get_energy("weno9pV_sixteen_new_postprocess.jld2", avg2)
qg232 = get_energy("qgleith_thirtytwo_new_postprocess.jld2", avg2)
wn232 = get_energy("weno9pV_thirtytwo_new_postprocess.jld2", avg2)

others08 = (om08, wd08, up08)
others16 = (om16, wd16, up16)

TKEglob08 = [sum(var.TKE) for var in others08]
TKEglob16 = [sum(var.TKE) for var in others16]

@show argmin(TKEglob08)
@show argmax(TKEglob08)
@show argmin(TKEglob16)
@show argmax(TKEglob16)

mn08 = others08[argmin(TKEglob08)]
mx08 = others08[argmax(TKEglob08)]
mn16 = others16[argmin(TKEglob16)]
mx16 = others16[argmax(TKEglob16)]

maxTKE08 = Float64[maximum(var.TKE[i] for var in others08) for i in 1:200]
minTKE08 = Float64[minimum(var.TKE[i] for var in others08) for i in 1:200]
maxTKE16 = Float64[maximum(var.TKE[i] for var in others16) for i in 1:200]
minTKE16 = Float64[minimum(var.TKE[i] for var in others16) for i in 1:200]

color08 = Makie.color("chocolate2") 
color16 = Makie.color("blueviolet")
color32 = Makie.color("steelblue")   

fig = Figure(resolution = (900, 400) .* 2, fontsize = 15 * 2)
ga = fig[1:2, 1:6] = GridLayout()
gb = fig[1:2, 7:8] = GridLayout()

ax = Axis(ga[2, 3:4], ylabel = L"\text{Total APE}", 
          xlabel = L"\text{days}",
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0.3, 0.5, 0.7], [L"0.3", L"0.5", L"0.7"]))

band!(ax, 1:200, mn08.TAPE ./ mn08.TAPE[1] .* 0.67, mx08.TAPE ./ mx08.TAPE[1] .* 0.67,  color = (color08, 0.4))
band!(ax, 1:200, mn16.TAPE ./ mn16.TAPE[1] .* 0.67, mx16.TAPE ./ mx16.TAPE[1] .* 0.67,  color = (color16, 0.4))
lines!(ax, mn08.TAPE ./ mn08.TAPE[1] .* 0.67,  linewidth = 1,   color = color08)
lines!(ax, mx08.TAPE ./ mx08.TAPE[1] .* 0.67,  linewidth = 1,   color = color08)
# lines!(ax, eb08.TAPE ./ eb08.TAPE[1] .* 0.67,  linewidth = 2,   color = color08)
lines!(ax, mn16.TAPE ./ mn16.TAPE[1] .* 0.67,  linewidth = 1,   color = color16)
lines!(ax, mx16.TAPE ./ mx16.TAPE[1] .* 0.67,  linewidth = 1,   color = color16)
lines!(ax, wn08.TAPE ./ wn08.TAPE[1] .* 0.67,  linewidth = 2,   color = color08, label = "let's go")
lines!(ax, wn16.TAPE ./ wn16.TAPE[1] .* 0.67,  linewidth = 2,   color = color16)
lines!(ax, wn32.TAPE ./ wn32.TAPE[1] .* 0.67,  linewidth = 2,   color = color32)
lines!(ax, wn208.TAPE ./ wn208.TAPE[1] .* 0.67, linewidth = 0.2, color = color08)
lines!(ax, wn216.TAPE ./ wn216.TAPE[1] .* 0.67, linewidth = 0.2, color = color16)
lines!(ax, wn232.TAPE ./ wn232.TAPE[1] .* 0.67, linewidth = 0.2, color = color32)

# leg = Legend(fig[1:2, 5], ax)

ax = Axis(ga[2, 1:2], ylabel = L"\text{Total KE}", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, minTKE08, maxTKE08,   color = (color08, 0.4))
band!(ax, 1:200, minTKE16, maxTKE16,   color = (color16, 0.4))
lines!(ax, mn08.TKE, linewidth = 1,    color = color08)
lines!(ax, mx08.TKE, linewidth = 1,    color = color08)
# lines!(ax, eb08.TKE, linewidth = 2,    color = color08)
lines!(ax, mn16.TKE, linewidth = 1,    color = color16)
lines!(ax, mx16.TKE, linewidth = 1,    color = color16)
lines!(ax, wn08.TKE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.TKE, linewidth = 2,    color = color16)
lines!(ax, wn32.TKE, linewidth = 2,    color = color32)
lines!(ax, wn208.TKE, linewidth = 0.2, color = color08)
lines!(ax, wn216.TKE, linewidth = 0.2, color = color16)
lines!(ax, wn232.TKE, linewidth = 0.2, color = color32)

ax = Axis(ga[1, 3:4], ylabel = L"\text{Eddy APE}", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.03, 0.06], [L"0.0", L"0.03", L"0.06"]))

band!(ax, 1:200, mn08.EAPE, mx08.EAPE,  color = (color08, 0.4))
band!(ax, 1:200, mn16.EAPE, mx16.EAPE,  color = (color16, 0.4))
lines!(ax, mn08.EAPE, linewidth = 1,    color = color08)
lines!(ax, mx08.EAPE, linewidth = 1,    color = color08)
lines!(ax, mn16.EAPE, linewidth = 1,    color = color16)
lines!(ax, mx16.EAPE, linewidth = 1,    color = color16)
lines!(ax, wn08.EAPE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.EAPE, linewidth = 2,    color = color16)
lines!(ax, wn32.EAPE, linewidth = 2,    color = color32)
lines!(ax, wn208.EAPE, linewidth = 0.2, color = color08)
lines!(ax, wn216.EAPE, linewidth = 0.2, color = color16)
lines!(ax, wn232.EAPE, linewidth = 0.2, color = color32)

ax = Axis(ga[1, 1:2], ylabel = L"\text{Eddy KE}", 
          xlabel = "", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

band!(ax, 1:200, mn08.EKE, mx08.EKE,   color = (color08, 0.4))
band!(ax, 1:200, mn16.EKE, mx16.EKE,   color = (color16, 0.4))
lines!(ax, mn08.EKE, linewidth = 1,    color = color08)
lines!(ax, mx08.EKE, linewidth = 1,    color = color08)
lines!(ax, mn16.EKE, linewidth = 1,    color = color16)
lines!(ax, mx16.EKE, linewidth = 1,    color = color16)
lines!(ax, wn08.EKE, linewidth = 2,    color = color08, label = "let's go")
lines!(ax, wn16.EKE, linewidth = 2,    color = color16)
lines!(ax, wn32.EKE, linewidth = 2,    color = color32)
lines!(ax, wn208.EKE, linewidth = 0.2, color = color08)
lines!(ax, wn216.EKE, linewidth = 0.2, color = color16)
lines!(ax, wn232.EKE, linewidth = 0.2, color = color32)

colgap!(ga, 5)
rowgap!(ga, 5)

ax = Axis(gb[1, 1], ylabel = "", 
          xlabel = "", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.1, 0.2], [L"0.0", L"0.1", L"0.2"]))

lines!(ax, qg208.TKE, linewidth = 1.5,    color = color08)
# lines!(ax, eb208.TKE, linewidth = 2.0,    color = color08)
lines!(ax, up208.TKE, linewidth = 2.5,  color = color08, label = "let's go")
lines!(ax, om208.TKE, linewidth = 0.5, color = color08)
lines!(ax, wd208.TKE, linewidth = 1.0, color = color08)

ax = Axis(gb[2, 1], ylabel = "", 
          xticks = ([0, 50, 100, 150, 200], ["", "", "", "", ""]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

lines!(ax, qg216.TKE, linewidth = 1.5,    color = color16)
# lines!(ax, qg16.TKE, linewidth = 1.5,    color = color16)
# lines!(ax, eb16.TKE, linewidth = 1,    color = color08)
lines!(ax, up216.TKE, linewidth = 2.5, color = color16, label = "let's go")
lines!(ax, om216.TKE, linewidth = 0.2, color = color16)
lines!(ax, wd216.TKE, linewidth = 1.0, color = color16)

ax = Axis(gb[3, 1], ylabel = "", 
          xlabel = L"\text{days}", 
          xticks = ([0, 50, 100, 150, 200], [L"0", L"250", L"500", L"750", L"1000"]),
          yticks = ([0., 0.2, 0.4], [L"0.0", L"0.2", L"0.4"]))

lines!(ax, qg232.TKE, linewidth = 1,    color = color32)
# lines!(ax, eb32.TKE, linewidth = 1,    color = color08)
# lines!(ax, up32.TKE, linewidth = 2.5,  color = color08, label = "let's go")
# lines!(ax, om32.TKE, linewidth = 0.2, color = color08)
lines!(ax, wn232.TKE, linewidth = 0.2, color = color32)

colgap!(ga, 5)
rowgap!(ga, 5)






# using CairoMakie
# CairoMakie.activate!()
 