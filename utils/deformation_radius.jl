
H = 1000

Δλ = 0.1
Δφ = 0.1

λ = 0:Δλ:20
φ = -60:Δφ:-40
z = -H:5:0

Δb = 5e-3
N² = 4e-6

γ(φ) = 2π * (φ + 60) / 20 - π/2

g(φ) = γ(φ) < 0 ? 0 :
       γ(φ) > π ? 1 :
       (γ(φ) - sin(γ(φ)) * cos(γ(φ))) / π

b(φ, z) = N² * z + Δb * g(φ)

bᵢ = [b(φ[i], z[j]) for i in 1:length(φ), j in 1:length(z)]

day = 24 * 3600
Ω = 2π / day
R = 6400e3

f(φ) = 2Ω * sind(φ)

uᵢ = 0*bᵢ
for i in 2:length(φ)-1, j in 1:length(z)
    ∂g∂φᵢ = (g(φ[i+1]) - g(φ[i-1])) / 2Δφ
    uᵢ[i, j] = - 180 / (π * R * f(φ[i])) * Δb * ∂g∂φᵢ * (z[j] + H)
end

using ColorSchemes

using GLMakie
using CairoMakie
GLMakie.activate!()
CairoMakie.activate!()

col = deepcopy(ColorSchemes.PRGn_5.colors)

new_col = []
for i in 1:2
    push!(new_col, RGBf(col[i].r, col[i].g, col[i].b * 0.65))
end
for i in 3:5
    push!(new_col, col[i])
end

mycmap = ColorScheme(typeof(new_col[1]).(new_col))

using JLD2
Rd = jldopen("restoring/deformation_radius.jld2")["deformation_radius"]

using Makie

function my_label_formatter(level::Real)::String 
    makie_string = Makie.contour_label_formatter(level)
    return LaTeXStrings.latexstring(makie_string)
end

fig = Figure(resolution = (1050, 500), fontsize = 20)
ax = Axis(fig[2:4, 1:4],
          xlabel = L"\textrm{Latitude}",
          ylabel = L"\textrm{Depth km}",
          xticks = ([-60, -55, -50, -45, -40], [L"60^\text{o} \text{S}", L"55^\text{o} \text{S}", L"50^\text{o} \text{S}", L"45^\text{o} \text{S}", L"40^\text{o} \text{S}"]),
          yticks = (collect(-1:0.5:0), [L"-1.0", L"-0.5", L"0"]))
hm = contourf!(ax, φ, z / 1e3, bᵢ, levels = range(-0.004, 0.006, length = 9), colormap = mycmap)
with_theme(theme_latexfonts()) do
    contour!(ax, φ, z / 1e3, uᵢ, linewidth = 3, color = :black, liestyle = :dashed, levels=0.01:0.02:0.08, labels = true, labelformatter = Makie.contour_label_formatter, labelsize = 20)
end
xlims!(ax, (-60, -40))
ylims!(ax, (-1, 0))

cb = Colorbar(fig[1, 1:4], hm, vertical = false, label = L"\text{Shading: initial buoyancy m/s}^2\text{. Contour: initial velocity m/s}",
              ticks = ([-0.0025, 0, 0.0025, 0.005], [L"-0.0025", L"0", L"0.0025", L"0.005"]))

ax = Axis(fig[1:4, 5:6],
          xlabel = L"\textrm{Time days}",
          ylabel = L"\textrm{Deformation radius km}",
          yticks = ([5.75, 6, 6.25, 6.5, 6.75], [L"5.75", L"6.00", L"6.25", L"6.50", L"6.75"]),
          xticks = ([0, 250, 500, 750, 1000], [L"0", L"250", L"500", L"750", L"1000"]))
lines!(ax, range(1, 1000, length = 201)[2:end], Rd[2:end] ./ 1000, linewidth = 3, color = :black)

save("initial_conditions.eps", fig)

# using Oceananigans
# using BaroclinicAdjustment
# using BaroclinicAdjustment.Diagnostics: DeformationRadius, all_fieldtimeseries, VolumeField

# fts = all_fieldtimeseries("../restoring/data/weno9pV_sixteen_new_snapshots.jld2");

# grid   = fts[:u].grid
# volume = VolumeField(grid, (Center, Center, Nothing))

# Rd = Float64[]
# for i in eachindex(fts[:b].times)
#     push!(Rd, sum(volume * DeformationRadius(fts, i)) / sum(volume))
# end

