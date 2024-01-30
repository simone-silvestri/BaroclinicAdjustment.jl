using BaroclinicAdjustment
using JLD2
using GeoMakie
using Oceananigans
using Oceananigans.Units

labels_to_plot = ["bilap",
                  "qgleith", 
                  "weno5pV",
                  "weno9pV", 
                  "weno9pAllD",
                  "leith",
                  "smag"]

label_all_three = ["qgleith", "weno9pV", "leith", "smag"] #, "weno5pV"]
label_fifty = ["qgleith"]

heatmaps_to_plot = ["qgleith", 
                    "weno9pV", 
                    "weno9pAllD"]

plot_titles_quarter  = [L"\text{QG }  0.25-\text{degree}", 
                        L"\text{W9V } 0.25-\text{degree}", 
                        L"\text{W9D } 0.25-\text{degree}"]

plot_titles_eight  = [L"\text{QG }  0.125-\text{degree}", 
                      L"\text{W9V } 0.125-\text{degree}", 
                      L"\text{W9D } 0.125-\text{degree}"]

plot_titles_sixteen  = [L"\text{QG }  0.0625-\text{degree}", 
                        L"\text{W9V } 0.0625-\text{degree}", 
                        L"\text{W9D } 0.0625-\text{degree}"]

function all_figures()
    fig1 = Figure(resolution = (900, 800))
    fig1 = plot_energies(fig1)

    res = ceil(Int, 450 * 0.8)   
    fig2 = Figure(resolution = (800, res))
    fig2 = plot_energy_stratif(fig2)

    res = ceil(Int, 400 * 0.8)   
    fig3 = Figure(resolution = (800, res))
    fig3 = plot_enstrophies(fig3)

    return fig1, fig2, fig3
end
                
function all_contours()
    fig1 = Figure(resolution = (5000, 5000))
    fig1 = energy_contours(fig1)

    fig2 = Figure(resolution = (5000, 5000))
    fig2 = vorticity_contours(fig2)

    fig3 = Figure(resolution = (5000, 5000))
    fig3 = enstrophy_contours(fig3)

    fig4 = Figure(resolution = (5000, 5000))
    fig4 = buoyancy_contours(fig4)

    return fig1, fig2, fig3, fig4
end

function all_contours_white()
    fig1 = Figure(resolution = (800, 800))
    fig1 = white_contours(fig1)

    return fig1
end

create_axis(fig, row, column, title, xlabel, ylabel, xticks, yticks = ([], []); yscale = identity) = 
    Axis(fig[row, column]; 
    yscale,
    xticks,
    yticks,
    xgridvisible = false, ygridvisible = false,
    xlabel, ylabel,
    title)

# xminorticksvisible = true, yminorticksvisible = true, 
# xticksmirrored = true,
# yticksmirrored = true,
# xminortickalign = 0,
# yminortickalign = 0,
# xtickalign = Makie.automatic,
# ytickalign = Makie.automatic,
# xminorticks = IntervalsBetween(9),
# yminorticks = IntervalsBetween(5),
    
function plot_stuff_all_three!(ax, E4, E8, E16, E50)

    days  = range(1, 200, length = 101)
    days2 = range(1, 200, length = 11)

    D4  = []
    D8  = []
    D16 = []
    for i in eachindex(E4)
        e4  = E4[i]
        e8  = E8[i]
        e16 = E16[i]
        push!(D4,  @lift($e4[1:10:end]))
        push!(D8,  @lift($e8[1:10:end]))
        push!(D16, @lift($e16[1:10:end]))
    end

    days5 = [length(E50[i].val) < 101 ?  range(1, 100, length = 26) : days for i in eachindex(E50)]

    lines!(ax, days5[1], E50[1], linewidth = 7, color = RGBf(0.8, 0.8, 0.8))
    if length(E50) > 1
        lines!(ax, days5[2], E50[2], linewidth = 7, color = :red)
    end
    if length(E50) > 2
        lines!(ax, days5[3], E50[3], linewidth = 7, color = :deepskyblue)
    end
    if length(E50) > 3
        lines!(ax, days5[4], E50[4], linewidth = 7, color = :purple)
    end
    
    lines!(ax, days, E4[1], linewidth = 0.75, color = :green3)
    scatter!(ax, days2, D4[1],  markersize = 10, marker = :rect,      color = Symbol("green3"))    
    lines!(ax, days, E4[2], linewidth = 0.75, color = :purple)
    scatter!(ax, days2, D4[2],  markersize = 10, marker = :ltriangle, color = Symbol("purple"))

    lines!(ax, days, E8[1], linewidth = 1.5, color = :green3)
    scatter!(ax, days2, D8[1],  markersize = 12.5, marker = :rect,      color = Symbol("green3"))    
    lines!(ax, days, E8[2], linewidth = 1.5, color = :purple)
    scatter!(ax, days2, D8[2],  markersize = 12.5, marker = :ltriangle, color = Symbol("purple"))

    lines!(ax, days, E16[1], linewidth = 3.0, color = :green3)
    scatter!(ax, days2, D16[1], markersize = 15, marker = :rect,      color = Symbol("green3"))    
    lines!(ax, days, E16[2], linewidth = 3.0, color = :purple)
    scatter!(ax, days2, D16[2], markersize = 15, marker = :ltriangle, color = Symbol("purple"))

    if length(E4) > 2
        lines!(ax, days, E4[3],  linewidth = 0.75, color = :red)
        lines!(ax, days, E8[3],  linewidth = 1.5,  color = :red)
        lines!(ax, days, E16[3], linewidth = 3.0,  color = :red)
    end
    if length(E4) > 3
        lines!(ax, days, E4[4],  linewidth = 0.75, color = :deepskyblue)
        lines!(ax, days, E8[4],  linewidth = 1.5,  color = :deepskyblue)
        lines!(ax, days, E16[4], linewidth = 3.0,  color = :deepskyblue)
    end

    if length(D4) > 2
        scatter!(ax, days2, D4[3], markersize  = 10,   marker = :rtriangle, color = Symbol("red"))
        scatter!(ax, days2, D8[3], markersize  = 12.5, marker = :rtriangle, color = Symbol("red"))
        scatter!(ax, days2, D16[3], markersize = 15,   marker = :rtriangle, color = Symbol("red"))
    end
    if length(D4) > 3
        scatter!(ax, days2, D4[4], markersize  = 10,   marker = :diamond, color = Symbol("deepskyblue"))
        scatter!(ax, days2, D8[4], markersize  = 12.5, marker = :diamond, color = Symbol("deepskyblue"))
        scatter!(ax, days2, D16[4], markersize = 15,   marker = :diamond, color = Symbol("deepskyblue"))
    end

    return nothing
end

function plot_stuff!(ax, E)
    days = range(1, 100, length = 101)

    lines!(ax, days, E[1], linewidth = 1, color = Symbol("orange1"),     label = "Bilaplacian")
    lines!(ax, days, E[2], linewidth = 1, color = Symbol("green3"),      label = "Leith")
    lines!(ax, days, E[3], linewidth = 1, color = Symbol("deepskyblue"), label = "QGLeith")
    lines!(ax, days, E[4], linewidth = 1, color = Symbol("purple"),      label = "Weno9pV")
    lines!(ax, days, E[5], linewidth = 1, color = Symbol("blue"),        label = "Weno9fV")
    if length(E) > 5
        lines!(ax, days, E[6], linewidth = 1, color = Symbol("red"),     label = "Weno9sV")
    end
    if length(E) > 6
        lines!(ax, days, E[7], linewidth = 1, color = Symbol("yellow"), label = "Weno9pD")
    end
    if length(E) > 7
        lines!(ax, days, E[8], linewidth = 1, color = Symbol("gray10"), label = "Weno9pD")
    end

    D = []
    for i in eachindex(E)
        e = E[i]
        push!(D, @lift($e[1:10:end]))
    end

    days = range(1, 100, length = 11)
    scatter!(ax, days, D[1], markersize = 10, marker = :circle,    color = Symbol("orange1"))   
    scatter!(ax, days, D[2], markersize = 10, marker = :rect,      color = Symbol("green3"))    
    scatter!(ax, days, D[3], markersize = 10, marker = :diamond,   color = Symbol("deepskyblue"))    
    scatter!(ax, days, D[4], markersize = 10, marker = :ltriangle, color = Symbol("purple"))
    scatter!(ax, days, D[5], markersize = 10, marker = :utriangle, color = Symbol("blue"))       
    if length(D) > 5
        scatter!(ax, days, D[6], markersize = 10, marker = :rtriangle, color = Symbol("red"))      
    end
    if length(D) > 6
        scatter!(ax, days, D[7], markersize = 10, marker = :star4, color = Symbol("yellow"))
    end
    if length(D) > 7
        scatter!(ax, days, D[8], markersize = 10, marker = :hexagon, color = Symbol("gray10"))
    end

    return nothing
end

function plot_stratif_RPE(; row, column, fig = nothing, title = "", 
                          xlabel = L"days", ylabel = "",
                          xlims = nothing, ylims = nothing, xticks = ([], []), yticks = ticks([]))

    iter = Observable(201)

    keys = label_all_three

    en4  = jldopen("energies_quarter.jld2")
    en8  = jldopen("energies_eight.jld2")
    en16 = jldopen("energies_sixteen.jld2")
    en50 = jldopen("energies_fifty.jld2")

    st4  = jldopen("stratif_quarter.jld2")
    st8  = jldopen("stratif_eight.jld2")
    st16 = jldopen("stratif_sixteen.jld2")
    st50 = jldopen("stratif_fifty.jld2")

    KE4  = []
    KE8  = []
    KE16 = []
    KE50 = []
    for i in eachindex(keys)
        push!(KE4,  @lift( st4[keys[i]][1:2:$iter] ./ ( en4[keys[i]].RPE[1:2:$iter] .-  en4[keys[i]].RPE[1])./ 4e-6 .* 1e17))
        push!(KE8,  @lift( st8[keys[i]][1:2:$iter] ./ ( en8[keys[i]].RPE[1:2:$iter] .-  en8[keys[i]].RPE[1])./ 4e-6 .* 1e17))
        push!(KE16, @lift(st16[keys[i]][1:2:$iter] ./ (en16[keys[i]].RPE[1:2:$iter] .- en16[keys[i]].RPE[1])./ 4e-6 .* 1e17))
    end
    for key in label_fifty
        myiter = length(st50[key]) < 201 ? Observable(51) : Observable(201)
        push!(KE50, @lift(st50[key][1:2:$myiter] ./ (en50[key].RPE[1:2:$myiter] .- en50[key].RPE[1]) ./ 4e-6 .* 1e17))
    end

    if isnothing(fig)
        fig = Figure(resolution = (1000, 500))
    end

    ax1 = create_axis(fig, row, column, title, xlabel, ylabel, xticks, yticks)

    plot_stuff_all_three!(ax1, KE4, KE8, KE16, KE50)
    # if isnothing(xlims)
    #     xlims!(ax1, (0, 100))
    # else        
    #     xlims!(ax1, xlims)
    # end

    # if !isnothing(ylims)
    #     ylims!(ax1, ylims)
    # end

    return fig
end

function plot_all_three(prefix; row, column, fig = nothing, title = "", 
                                xlabel = L"days", ylabel = "", normalize = true, 
                                normalize_val = 1, energy = false, 
                                xlims = nothing, ylims = nothing, xticks = ([], []), yticks = ticks([]))
    
    iter = Observable(401)

    keys = label_all_three

    file4  = jldopen(prefix * "_quarter.jld2")
    file8  = jldopen(prefix * "_eight.jld2")
    file16 = jldopen(prefix * "_sixteen.jld2")
    file50 = jldopen(prefix * "_fifty.jld2")

    KE4  = []
    KE8  = []
    KE16 = []
    KE50 = []
    if energy
        for i in eachindex(keys)
            push!(KE4,  !(normalize) ? @lift( file4[keys[i]].KE[1:4:$iter]) : @lift( file4[keys[i]].KE[1:4:$iter] ./ normalize_val))
            push!(KE8,  !(normalize) ? @lift( file8[keys[i]].KE[1:4:$iter]) : @lift( file8[keys[i]].KE[1:4:$iter] ./ normalize_val))
            push!(KE16, !(normalize) ? @lift(file16[keys[i]].KE[1:4:$iter]) : @lift(file16[keys[i]].KE[1:4:$iter] ./ normalize_val))
        end
        for key in label_fifty
            myiter = length(file50[key].KE) < 201 ? Observable(101) : Observable(201)
            push!(KE50, !(normalize) ? @lift(file50[key].KE[1:4:$myiter]) : @lift(file50[key].KE[1:4:$myiter] ./ normalize_val))
        end
    else
        for i in eachindex(keys)
            push!(KE4,  !(normalize) ? @lift( file4[keys[i]][1:4:$iter]) : @lift( file4[keys[i]][1:4:$iter] ./ normalize_val))
            push!(KE8,  !(normalize) ? @lift( file8[keys[i]][1:4:$iter]) : @lift( file8[keys[i]][1:4:$iter] ./ normalize_val))
            push!(KE16, !(normalize) ? @lift(file16[keys[i]][1:4:$iter]) : @lift(file16[keys[i]][1:4:$iter] ./ normalize_val))
        end
        for key in label_fifty
            myiter = length(file50[key]) < 201 ? Observable(101) : Observable(201)
            push!(KE50, !(normalize) ? @lift(file50[key][1:4:$myiter]) : @lift(file50[key][1:4:$myiter] ./ normalize_val))
        end
    end

    if isnothing(fig)
        fig = Figure(resolution = (1000, 500))
    end

    ax1 = create_axis(fig, row, column, title, xlabel, ylabel, xticks, yticks)

    plot_stuff_all_three!(ax1, KE4, KE8, KE16, KE50)
    if isnothing(xlims)
        xlims!(ax1, (0, 100))
    else        
        xlims!(ax1, xlims)
    end

    if !isnothing(ylims)
        ylims!(ax1, ylims)
    end

    return fig
end

function plot_file(prefix, filename; 
                    row = 1, column = 1, 
                    fig = nothing, title = "", 
                    xlabel = "days", ylabel = "", normalize = true, 
                    normalize_val = 1, dissipation = false,
                    xlims = nothing, ylims = nothing,
                    xticks = ([], []), yticks = ([], []))

    iter = Observable(201)

    keys = labels_to_plot
    file = jldopen(prefix * filename * ".jld2")

    KE = []
    for i in eachindex(keys)
        if dissipation
            push!(KE, @lift(- file[keys[i]].time_derivative[1:1:($iter÷2)] ./ file[keys[i]].dissipation[1:1:($iter÷2)]))
        else
            push!(KE, !(normalize) ? @lift(file[keys[i]][1:2:$iter]) : @lift(file[keys[i]][1:2:$iter] ./ normalize_val))
        end
    end

    if isnothing(fig)
        fig = Figure(resolution = (1000, 500))
    end

    ax1 = create_axis(fig, row, column, title, xlabel, ylabel, xticks, yticks)

    plot_stuff!(ax1, KE)
    if isnothing(xlims)
        xlims!(ax1, (0, 100))
    else        
        xlims!(ax1, xlims)
    end

    if !isnothing(ylims)
        ylims!(ax1, ylims)
    end

    return fig
end

ticks(a::Vector) = (a, [L"%$b" for b in a])

function plot_energies(filename; row = 1, fig = nothing, xlabel = "", legend = true, normalize = true, xticks = ([], []), yticks = (ticks([]) for i in 1:3))
    iter = Observable(201)

    keys = []

    for l in labels_to_plot
        push!(keys, l)
    end

    file  = jldopen("energies" * filename * ".jld2")
    strat = jldopen("stratif" * filename * ".jld2")

    KE = []
    ST = []
    for i in eachindex(keys)
        push!(KE, !(normalize) ? @lift(file[keys[i]].KE[1:2:$iter] .* 1e3) : @lift(file[keys[i]].KE[1:2:$iter] .* 1e3 ./ 1e17))
    end
    for i in eachindex(keys)
        push!(ST, @lift(strat[keys[i]][1:2:$iter] ./ 4e-6))
    end

    if isnothing(fig)
        fig = Figure(resolution = (1000, 500))
    end

    title_ke  = L"\text{KE}" 
    title_rpe = L"\text{RPE} - \text{RPE}_0"
    title_frac = L"\text{KE}/(\text{RPE} - \text{RPE}_0)"
    ax1 = create_axis(fig, row, 1, (row == 1 ? title_ke : ""), xlabel, "", xticks, yticks[1])

    plot_stuff!(ax1, KE)
    xlims!(ax1, (0, 100))

    RPE = []
    for i in eachindex(keys)
        push!(RPE, !(normalize) ? @lift(file[keys[i]].RPE[1:2:$iter] .- file[keys[i]].RPE[1]) : @lift((file[keys[i]].RPE[1:2:$iter] .- file[keys[i]].RPE[1]) ./ 1e17))
    end

    ax2 = create_axis(fig, row, 2, (row == 1 ? title_rpe : ""), xlabel, "", xticks, yticks[2])

    plot_stuff!(ax2, RPE)
    xlims!(ax2, (0, 100))
    
    DE = []
    for i in eachindex(keys)
        ke  = KE[i]
        rpe = RPE[i]
        push!(DE, @lift($ke ./ $rpe))
    end
    
    ax3 = create_axis(fig, row, 3, (row == 1 ? title_frac : ""), xlabel, "", xticks, yticks[3])

    plot_stuff!(ax3, DE)
    xlims!(ax3, (10, 100))
    ylims!(ax3, (DE[5].val[end] * 0.5, DE[3].val[end - 50] * 1.5))
    legend && axislegend(ax3, position = :rt, framevisible = false)

    return fig
end  

Lz = 1kilometers

grid4 = LatitudeLongitudeGrid(topology = (Periodic, Bounded, Bounded),
                              size = (80, 80, 50), 
                              longitude = (-10, 10),
                              latitude = (-60, -40),
                              z = (-Lz, 0))

grid8 = LatitudeLongitudeGrid(topology = (Periodic, Bounded, Bounded),
                              size = (160, 160, 50), 
                              longitude = (-10, 10),
                              latitude = (-60, -40),
                              z = (-Lz, 0))

grid16 = LatitudeLongitudeGrid(topology = (Periodic, Bounded, Bounded),
                               size = (320, 320, 50), 
                               longitude = (-10, 10),
                               latitude = (-60, -40),
                               z = (-Lz, 0))

vol4  = BaroclinicAdjustment.Diagnostics.VolumeField(grid4,  (Face, Center, Center))                 
vol8  = BaroclinicAdjustment.Diagnostics.VolumeField(grid8,  (Face, Center, Center))                     
vol16 = BaroclinicAdjustment.Diagnostics.VolumeField(grid16, (Face, Center, Center))    

mean_vol4  = sum(interior(vol4,  2:80,  :, :), dims = 3)[:, :, 1]
mean_vol8  = sum(interior(vol8,  2:160, :, :), dims = 3)[:, :, 1]
mean_vol16 = sum(interior(vol16, 2:320, :, :), dims = 3)[:, :, 1]

function heatmap_yzplane(filename, variable; fig = nothing, 
    row = 1, maxrange = nothing, maxtol = 0.6,
    minrange = 0, colormap = :viridis,
    hide_spines = true,
    plot_titles = false, level = 1)

    file = jldopen(filename * ".jld2")

    @show(keys(file))

    titles = labels_to_plot
    map = Vector(undef, length(titles))

    for i in eachindex(titles)
        map[i] = file[variable * "/" * titles[i]][1]
    end

    @show Nx, Ny = size(map[1])
    x = range(-10,  10, length = Nx)
    y = range(-1000, 0, length = Ny)

    if isnothing(maxrange)
        maxrange = maximum(maximum.(map)) * maxtol
    end
    if isnothing(minrange)
        minrange = minimum(minimum.(map)) * 1 / 0.60
    end
    dlevel = (maxrange - minrange) / level
    levels = collect(minrange:dlevel:maxrange)

    if isnothing(fig)
        fig = Figure(resolution = (3000, 500))
    end
    
    for i in 1:length(map) - 1
        if plot_titles
            ax = Axis(fig[row, i]; title = titles[i])
        else
            ax = Axis(fig[row, i])
        end
        contourf!(ax, x, y, map[i]; levels = levels, colormap)
        contour!(ax, x, y, map[i]; levels = levels, color = :black)
        hidedecorations!(ax)
        xlims!(ax, (-10, 10))
        ylims!(ax, (-1000, 0))
        hide_spines && hidespines!(ax)
    end

    i = length(map)
    if plot_titles
        ax = Axis(fig[row, i]; title = titles[i])
    else
        ax = Axis(fig[row, i])
    end
    hm = contourf!(ax, x, y, map[i]; levels = levels, colormap)
    contour!(ax, x, y, map[i]; levels = levels, color = :black)
    hidedecorations!(ax)
    xlims!(ax, (-10, 10))
    ylims!(ax, (-1000, 0))
    hide_spines && hidespines!(ax)

    Colorbar(fig[row, length(map)+1], hm)

    return fig
end

function heatmap_energies(filename = "", variable = nothing; fig = nothing, 
                          row = 1, maxrange = nothing, maxtol = 0.6,
                          minrange = 0, colormap = :viridis, mean_vol = 1,
                          hide_spines = false,  normalize_val = 1.0,
                          plot_titles = nothing, level = 1)

    titles = heatmaps_to_plot
    map = Vector(undef, length(titles))

    if filename isa WhiteContour
        for i in eachindex(titles)
            map[i] = ones(10, 10)
        end
    else
        file = jldopen(filename * ".jld2")
        @show(keys(file))
        for i in eachindex(titles)
            if isnothing(variable)
                map[i] = file[titles[i]].BTKE[1][:, :, 1] ./ mean_vol 
            else
                map[i] = file[variable * "/" * titles[i]][level][:, :, 1] ./ normalize_val
            end
        end
    end

    @show Nx, Ny = size(map[1])
    x = range(-10,  10, length = Nx)
    y = range(-60, -40, length = Ny)

    x2 = range(-10,  10, length = Nx * 8)
    y2 = range(-60, -40, length = Ny * 8)

    map2 = Vector(undef, length(titles))

    for i in eachindex(titles)
        map2[i] = repeat(map[i],  inner = (8, 1))
        map2[i] = repeat(map2[i], inner = (1, 8))
    end

    if isnothing(maxrange)
        maxrange = maximum(maximum.(map)) * maxtol
    end

    if isnothing(minrange)
        minrange = minimum(minimum.(map)) * 1 / 0.60
    end
    colorrange = (minrange, maxrange)

    if isnothing(fig)
        fig = Figure(resolution = (3000, 500))
    end
    
    for i in 1:length(map) - 1
        if !isnothing(plot_titles)
            ax = GeoAxis(fig[row, i]; title = plot_titles[i], dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))
        else
            ax = GeoAxis(fig[row, i]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))
        end
        !(filename isa WhiteContour) && surface!(ax, x2, y2, map2[i]; colorrange, colormap, shading = false)
        hidedecorations!(ax)
        hide_spines && hidespines!(ax)
    end

    i = length(map)
    if !isnothing(plot_titles)
        ax = GeoAxis(fig[row, i]; title = plot_titles[i], dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))
    else
        ax = GeoAxis(fig[row, i]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))
    end
    !(filename isa WhiteContour) && surface!(ax, x, y, map[i]; colorrange, colormap, shading = false)
    hidedecorations!(ax)
    hide_spines && hidespines!(ax)

    figtmp = Figure()
    axtmp = GeoAxis(figtmp[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))
    hm    = surface!(axtmp, x, y, map[i]; colorrange, colormap, shading = false)
    
    Colorbar(fig[row, length(map)+1], hm, width = 10, height = Relative(2/3))

    return fig
end

function BTenergy_contours(fig)
    fig = heatmap_energies("energies_quarter"; fig, row = 1, mean_vol = mean_vol4,  colormap = :viridis,               plot_titles = plot_titles_quarter)
    fig = heatmap_energies("energies_eight";   fig, row = 2, mean_vol = mean_vol8,  colormap = :viridis,               plot_titles = plot_titles_eight  )
    fig = heatmap_energies("energies_sixteen"; fig, row = 3, mean_vol = mean_vol16, colormap = :viridis, maxtol = 0.2, plot_titles = plot_titles_sixteen)
    return fig
end

function energy_contours(fig)
    fig = heatmap_energies("contours_quarter", "energy"; fig, row = 1, normalize_val = mean_vol4,  colormap = :viridis,               plot_titles = plot_titles_quarter)
    fig = heatmap_energies("contours_eight", "energy";   fig, row = 2, normalize_val = mean_vol8,  colormap = :viridis,               plot_titles = plot_titles_eight  )
    fig = heatmap_energies("contours_sixteen", "energy"; fig, row = 3, normalize_val = mean_vol16, colormap = :viridis, maxtol = 0.2, plot_titles = plot_titles_sixteen)
    return fig
end

function vorticity_contours(fig, level = 1)
    fig = heatmap_energies("contours_quarter", "vorticity"; fig, row = 1, level, colormap = Reverse(:balance), normalize_val = 1e-6, maxrange = 7.5, minrange = - 7.5, plot_titles = plot_titles_quarter)
    fig = heatmap_energies("contours_eight", "vorticity";   fig, row = 2, level, colormap = Reverse(:balance), normalize_val = 1e-6, maxrange = 7.5, minrange = - 7.5, plot_titles = plot_titles_eight  )
    fig = heatmap_energies("contours_sixteen", "vorticity"; fig, row = 3, level, colormap = Reverse(:balance), normalize_val = 1e-6, maxrange = 7.5, minrange = - 7.5, plot_titles = plot_titles_sixteen)
    return fig
end

function enstrophy_contours(fig, level = 1)
    fig = heatmap_energies("contours_quarter", "enstrophy"; fig, row = 1, level, colormap = :solar, maxrange = 2e-1, normalize_val = 1e-10, plot_titles = plot_titles_quarter)
    fig = heatmap_energies("contours_eight", "enstrophy";   fig, row = 2, level, colormap = :solar, maxrange = 1,    normalize_val = 1e-10, plot_titles = plot_titles_eight  )
    fig = heatmap_energies("contours_sixteen", "enstrophy"; fig, row = 3, level, colormap = :solar, maxrange = 2,    normalize_val = 1e-10, plot_titles = plot_titles_sixteen)
    return fig
end

struct WhiteContour end

function white_contours(fig, level = 1)
    fig = heatmap_energies(WhiteContour(); fig, row = 1, level, colormap = :solar, maxrange = 1, normalize_val = 1e-10, plot_titles = plot_titles_quarter)
    fig = heatmap_energies(WhiteContour(); fig, row = 2, level, colormap = :solar, maxrange = 1, normalize_val = 1e-10, plot_titles = plot_titles_eight  )
    fig = heatmap_energies(WhiteContour(); fig, row = 3, level, colormap = :solar, maxrange = 1, normalize_val = 1e-10, plot_titles = plot_titles_sixteen)
    return fig
end

function buoyancy_contours(fig, level = 2)
    fig = heatmap_energies("contours_quarter", "buoyancy"; fig, row = 1, level, colormap = :thermometer, maxrange = 0.006, minrange = 0.002, plot_titles = plot_titles_quarter)
    fig = heatmap_energies("contours_eight", "buoyancy";   fig, row = 2, level, colormap = :thermometer, maxrange = 0.006, minrange = 0.002, plot_titles = plot_titles_eight  )
    fig = heatmap_energies("contours_sixteen", "buoyancy"; fig, row = 3, level, colormap = :thermometer, maxrange = 0.006, minrange = 0.002, plot_titles = plot_titles_sixteen)
    return fig
end

function yzbuoyancy_contours(level = 1)
    fig = heatmap_yzplane("contours_quarter", "yzbuoyancy2"; hide_spines = false, level, colormap = :thermometer, plot_titles = true, maxrange = 0.006, minrange = -0.004)
    fig = heatmap_yzplane("contours_eight", "yzbuoyancy2";   hide_spines = false, level, fig, row = 2, colormap = :thermometer, maxrange = 0.006,       minrange = -0.004)
    fig = heatmap_yzplane("contours_sixteen", "yzbuoyancy2"; hide_spines = false, level, fig, row = 3, colormap = :thermometer, maxrange = 0.006,       minrange = -0.004)
    return fig
end

function refpe_contours(level = 1)
    fig = heatmap_energies("contours_quarter", "referencepe"; level, colormap = :thermometer, plot_titles = true)
    fig = heatmap_energies("contours_eight", "referencepe"; level, fig, row = 2, colormap = :thermometer) 
    fig = heatmap_energies("contours_sixteen", "referencepe"; level, fig, row = 3, colormap = :thermometer)
    return fig
end

days_ticks      = ([0, 25, 50, 75, 100], [L"0", L"25", L"50", L"75", L"100"])
E_quarter_ticks = (ticks([0, 1.0, 2.0]), ticks([0, 0.1, 0.2, 0.3]), ticks([2.0, 4.0, 6.0]))
E_eight_ticks   = (ticks([0, 2.0, 4.0]), ticks([0, 0.1, 0.2, 0.3]), ticks([5, 10, 15, 20, 25]))
E_sixteen_ticks = (ticks([0, 2.0, 4.0, 6.0]), ticks([0, 0.1, 0.2, 0.3]), ticks([40, 60, 80]))

stratif_ticks = ticks([1.0, 1.2, 1.4, 1.6, 1.8])

Ω_quarter_ticks = ticks([0, 0.1, 0.2, 0.3])
Ω_eight_ticks   = ticks([0, 0.2, 0.4, 0.6])
Ω_sixteen_ticks = ticks([0, 0.5, 1.0])

include("deformation_radius.jl")

function plot_deformation(fig)
    r4 = 22.765
    r8 = r4 / 2
    r16 = r8 / 2
    r50 = r8 / 50 * 8

    ax = create_axis(fig, 1, 1, L"\text{average  }L_d", L"\text{days}", L"\text{kilometers}", 
                     days_ticks, ticks([1, 2, 4, 8, 16, 32]), yscale = log2)
    lines!(ax, def_time, def_radius, linewidth = 7, color = RGBf(0.8, 0.8, 0.8))
    lines!(ax, [0, 100], [r4, r4], linewidth = 3, color = RGBf(0.5, 0.5, 0.5), linestyle = :dash)
    text!(ax, 95, r4 * 0.85, text = L"\text{mean spacing at 0.25-degree resolution}", align = (:right, :center))
    lines!(ax, [0, 100], [r8, r8], linewidth = 3, color = RGBf(0.5, 0.5, 0.5), linestyle = :dash)
    text!(ax, 95, r8 * 0.85, text = L"\text{mean spacing at 0.125-degree resolution}", align = (:right, :center))
    lines!(ax, [0, 100], [r16, r16], linewidth = 3, color = RGBf(0.5, 0.5, 0.5), linestyle = :dash)
    text!(ax, 95, r16 * 0.85, text = L"\text{mean spacing at 0.0625-degree resolution}", align = (:right, :center))
    lines!(ax, [0, 1600], [r50, r50], linewidth = 3, color = RGBf(0.5, 0.5, 0.5), linestyle = :dash)
    text!(ax, 95, r50 * 0.85, text = L"\text{mean spacing at 0.02-degree resolution}", align = (:right, :center))
    xlims!(ax, (0, 100))
    ylims!(ax, (1, 32))
    return fig
end

function plot_energies(fig)
    fig = plot_energies("_quarter"; fig, row = 1, legend = false, xlabel = "",             xticks = ([0, 25, 50, 75, 100], ["", "", "", "", ""]), yticks = E_quarter_ticks)
    fig = plot_energies("_eight";   fig, row = 2, legend = false, xlabel = "",             xticks = ([0, 25, 50, 75, 100], ["", "", "", "", ""]), yticks = E_eight_ticks)
    fig = plot_energies("_sixteen"; fig, row = 3, legend = false, xlabel = L"\text{days}", xticks = days_ticks,                                   yticks = E_sixteen_ticks)
    return fig
end

function plot_energy_stratif(fig)
    fig = plot_all_three("energies"; fig, row = 1, column = 1, energy = true, title = L"\text{KE}", 
                                     normalize_val = 1e14, yticks = ticks([0.0, 2.5, 5.0, 7.5]), xticks = days_ticks, xlabel = L"\text{days}")
    fig = plot_all_three("stratif"; fig, normalize_val = 4e-6, row = 1, column = 2, energy = false, 
                                    title = L"N^2 / N^2_0", xticks = days_ticks, yticks = stratif_ticks, xlabel = L"\text{days}")
    # fig = plot_stratif_RPE(; fig, row = 1, column = 3, xticks = days_ticks, yticks = stratif_ticks, xlabel = L"\text{days}")
    return fig
end

function plot_enstrophies(fig)
    fig = plot_file("enstrophies", "_quarter"; fig, row = 1, column = 1, normalize = true, xlabel = L"\text{days}", title = L"\int_V \zeta^2 dV / 10^6", normalize_val = 1e6, xticks = days_ticks, yticks = Ω_quarter_ticks)
    fig = plot_file("enstrophies","_eight";    fig, row = 1, column = 2, normalize = true, xlabel = L"\text{days}", title = L"\int_V \zeta^2 dV / 10^6", normalize_val = 1e6, xticks = days_ticks, yticks = Ω_eight_ticks)
    fig = plot_file("enstrophies", "_sixteen"; fig, row = 1, column = 3, normalize = true, xlabel = L"\text{days}", title = L"\int_V \zeta^2 dV / 10^6", normalize_val = 1e6, xticks = days_ticks, yticks = Ω_sixteen_ticks)
    return fig
end

function plot_dissipation()
    fig = plot_file("budgetB", "_quarter"; dissipation = true, xlims = (50, 200), normalize = true, normalize_val = 1e6)
    fig = plot_file("budgetB","_eight";    dissipation = true, xlims = (50, 200), fig, row = 1, column = 2, normalize = true, normalize_val = 1e6)
    fig = plot_file("budgetB", "_sixteen"; dissipation = true, xlims = (50, 200), fig, row = 1, column = 3, normalize = true, normalize_val = 1e6)
    return fig
end

function plot_stratif()
    fig = plot_file("stratif", "_quarter"; normalize = true, normalize_val = 1e-6, title = L"N^2 \cdot 10^6")
    fig = plot_file("stratif","_eight";    fig, row = 1, column = 2, normalize = true, normalize_val = 1e-6, title = L"N^2 \cdot 10^6")
    fig = plot_file("stratif", "_sixteen"; fig, row = 1, column = 3, normalize = true, normalize_val = 1e-6, title = L"N^2 \cdot 10^6")
    return fig
end
