using GLMakie
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: VerticalVorticity, KineticEnergy, DeformationRadius

using Oceananigans
using Oceananigans.Units
using JLD2

add_trailing_name(name) = name * "_snapshots.jld2"

function record_video!(name, fig, iter, Nt) 
    GLMakie.record(fig, name * ".mp4", 1:Nt, framerate = 11) do i
        @info "step $i"; 
        iter[] = i; 
    end
end

function surface_videos(trailing_character = "_weaker")
    file_prefix = ["weno7vdMD"] #, "weno9vd3"] #"weno5vd", key8leith, "lapleith", key8bilap,
                   #"smag", key8weno, "weno9", "weno9dd", "highres"]

    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    iter = Observable(1)

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)
            f = all_fieldtimeseries(filename)

            Nt = length(f[:b].times)
            @info "vorticity video $filename"    
            ζ = @lift(interior(VerticalVorticity(f, $iter), :, :, 50))

            fig = Figure()
            ax  = Axis(fig[1, 1])
            heatmap!(ax, ζ, colorrange = (-1e-5, 1e-5))
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("vort_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
            
            @info "kinetic energy video $filename"    
            E = @lift(interior(KineticEnergy(f, $iter), :, :, 50))

            fig = Figure()
            ax  = Axis(fig[1, 1])
            heatmap!(ax, E, colorrange = (0, 1.0), colormap = :solar)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("energy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            @info "buoyancy video $filename"    
            b = @lift(interior(f[:b][$iter], :, :, 50))

            fig = Figure()
            ax  = Axis(fig[1, 1])
            heatmap!(ax, b, colorrange = (0, 0.006), colormap = :thermometer)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("buoyancy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            # @info "deformation radius video $filename"    
            # R = @lift(interior(DeformationRadius(f, $iter), :, :, 1))

            # fig = Figure()
            # ax  = Axis(fig[1, 1])
            # heatmap!(ax, R, colorrange = (1e1, 1e4), colormap = :thermometer)
            # hidedecorations!(ax)
            # hidespines!(ax)
            # record_video!("Lr_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
        end
    end
end

function plot_stuff!(ax, Ew4, El4, Eb4, Ew16, El16, Eb16, Pw4, 
                     Pl4, Pb4, Pw16, Pl16, Pb16, iter_arr,
                     color1, color2, color3)
    lines!(ax, Ew4, linewidth = 2, color = color3, label = "3D WenoR")
    lines!(ax, El4, linewidth = 2, color = color1, label = "Leith")
    # lines!(ax, Eb4, linewidth = 2, color = color2, label = "Biharmonic")
    lines!(ax, Ew16, linewidth = 2, linestyle = :dash, color = color3)
    lines!(ax, El16, linewidth = 2, linestyle = :dash, color = color1)
    # lines!(ax, Eb16, linewidth = 2, linestyle = :dash, color = color2)

    scatter!(ax, iter_arr, Pw4, color = color3, markersize = 8)
    scatter!(ax, iter_arr, Pl4, color = color1, markersize = 8)
    # scatter!(ax, iter_arr, Pb4, color = color2, markersize = 12)
    scatter!(ax, iter_arr, Pw16, color = color3, markersize = 8)
    scatter!(ax, iter_arr, Pl16, color = color1, markersize = 8)
    # scatter!(ax, iter_arr, Pb16, color = color2, markersize = 12)
end

function plot_energies()
    iter = Observable(201)
    iter_arr = @lift([$iter])

    color1 = :deepskyblue
    color2 = :orange1
    color3 = :firebrick2

    key8weno = "weno5dd"
    key8leith = "leith"
    key8bilap = "bilap"
    key16weno = "weno9"
    key16leith = "smag"
    key16bilap = "lapleith"

    en4  = jldopen("energies_quarter.jld2")
    en16 = jldopen("energies_quarter.jld2")

    Ew4  = @lift(  en4[key8weno].KE[1:$iter] .* 1e3)
    El4  = @lift( en4[key8leith].KE[1:$iter] .* 1e3)
    Eb4  = @lift(    en4[key8bilap].KE[1:$iter] .* 1e3)
    Ew16 = @lift( en16[key16weno].KE[1:$iter] .* 1e3)
    El16 = @lift(en16[key16leith].KE[1:$iter] .* 1e3)
    Eb16 = @lift(   en16[key16bilap].KE[1:$iter] .* 1e3)
    Pw4  = @lift(  en4[key8weno].KE[$iter:$iter] .* 1e3)
    Pl4  = @lift( en4[key8leith].KE[$iter:$iter] .* 1e3)
    Pb4  = @lift(    en4[key8bilap].KE[$iter:$iter] .* 1e3)
    Pw16 = @lift( en16[key16weno].KE[$iter:$iter] .* 1e3)
    Pl16 = @lift(en16[key16leith].KE[$iter:$iter] .* 1e3)
    Pb16 = @lift(   en16[key16bilap].KE[$iter:$iter] .* 1e3)

    figE = Figure(resolution = (400, 300))
    ax1 = Axis(figE[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              yticks = ([], []),
              title = "Integrated Kinetic Energy")
        
    plot_stuff!(ax1, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                    Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                    iter_arr, color1, color2, color3)
    axislegend(ax1, position = :lt, framevisible = false)
    record_video!("energy_video2", figE, iter, 201)

    Ew4  = @lift(  en4[key8weno].RPE[1:$iter] .-  en4[key8weno].RPE[1]) 
    El4  = @lift( en4[key8leith].RPE[1:$iter] .- en4[key8leith].RPE[1]) 
    Eb4  = @lift(    en4[key8bilap].RPE[1:$iter] .-    en4[key8bilap].RPE[1]) 
    Ew16 = @lift( en16[key16weno].RPE[1:$iter] .- en16[key16weno].RPE[1]) 
    El16 = @lift(en16[key16leith].RPE[1:$iter] .-en16[key16leith].RPE[1]) 
    Eb16 = @lift(   en16[key16bilap].RPE[1:$iter] .-   en16[key16bilap].RPE[1]) 
    Pw4  = @lift(  en4[key8weno].RPE[$iter:$iter] .-  en4[key8weno].RPE[1]) 
    Pl4  = @lift( en4[key8leith].RPE[$iter:$iter] .- en4[key8leith].RPE[1]) 
    Pb4  = @lift(    en4[key8bilap].RPE[$iter:$iter] .-    en4[key8bilap].RPE[1]) 
    Pw16 = @lift( en16[key16weno].RPE[$iter:$iter] .- en16[key16weno].RPE[1]) 
    Pl16 = @lift(en16[key16leith].RPE[$iter:$iter] .-en16[key16leith].RPE[1]) 
    Pb16 = @lift(   en16[key16bilap].RPE[$iter:$iter] .-   en16[key16bilap].RPE[1]) 

    figR = Figure(resolution = (400, 300))
    ax2 = Axis(figR[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              yticks = ([], []),
              title = "Integrated Reference PE")
        
    plot_stuff!(ax2, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    record_video!("RPE_video2", figR, iter, 201)

    Ew4  = @lift(  en4[key8weno].APE[1:$iter] .-  en4[key8weno].APE[1]) 
    El4  = @lift( en4[key8leith].APE[1:$iter] .- en4[key8leith].APE[1]) 
    Eb4  = @lift(    en4[key8bilap].APE[1:$iter] .-    en4[key8bilap].APE[1]) 
    Ew16 = @lift( en16[key16weno].APE[1:$iter] .- en16[key16weno].APE[1]) 
    El16 = @lift(en16[key16leith].APE[1:$iter] .-en16[key16leith].APE[1]) 
    Eb16 = @lift(   en16[key16bilap].APE[1:$iter] .-   en16[key16bilap].APE[1]) 
    Pw4  = @lift(  en4[key8weno].APE[$iter:$iter] .-  en4[key8weno].APE[1]) 
    Pl4  = @lift( en4[key8leith].APE[$iter:$iter] .- en4[key8leith].APE[1]) 
    Pb4  = @lift(    en4[key8bilap].APE[$iter:$iter] .-    en4[key8bilap].APE[1]) 
    Pw16 = @lift( en16[key16weno].APE[$iter:$iter] .- en16[key16weno].APE[1]) 
    Pl16 = @lift(en16[key16leith].APE[$iter:$iter] .-en16[key16leith].APE[1]) 
    Pb16 = @lift(   en16[key16bilap].APE[$iter:$iter] .-   en16[key16bilap].APE[1]) 

    figR = Figure(resolution = (400, 300))
    ax2 = Axis(figR[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              yticks = ([], []),
              title = "Integrated Available PE")
        
    plot_stuff!(ax2, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    record_video!("APE_video2", figR, iter, 201)

    close(en4)
    close(en16)
end

function plot_energies_together()
    iter = Observable(201)
    iter_arr = @lift([$iter])

    color1 = :deepskyblue
    color2 = :orange1
    color3 = :firebrick2

    en4  = jldopen("energies_quarter.jld2")
    en16 = jldopen("energies_sixteen.jld2")

    KEw4  = @lift(  en4[key8weno].KE[1:$iter] .* 1e3)
    KEl4  = @lift( en4[key8leith].KE[1:$iter] .* 1e3)
    KEb4  = @lift(    en4[key8bilap].KE[1:$iter] .* 1e3)
    KEw16 = @lift( en16[key8weno].KE[1:$iter] .* 1e3)
    KEl16 = @lift(en16[key8leith].KE[1:$iter] .* 1e3)
    KEb16 = @lift(   en16[key8bilap].KE[1:$iter] .* 1e3)
    KPw4  = @lift(  en4[key8weno].KE[$iter:$iter] .* 1e3)
    KPl4  = @lift( en4[key8leith].KE[$iter:$iter] .* 1e3)
    KPb4  = @lift(    en4[key8bilap].KE[$iter:$iter] .* 1e3)
    KPw16 = @lift( en16[key8weno].KE[$iter:$iter] .* 1e3)
    KPl16 = @lift(en16[key8leith].KE[$iter:$iter] .* 1e3)
    KPb16 = @lift(   en16[key8bilap].KE[$iter:$iter] .* 1e3)

    figE = Figure(resolution = (400, 300))
    ax1 = Axis(figE[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Kinetic Energy")
        
    plot_stuff!(ax1, KEw4, KEl4, KEb4, KEw16, KEl16, KEb16, 
                    KPw4, KPl4, KPb4, KPw16, KPl16, KPb16, 
                    iter_arr, color1, color2, color3)
    axislegend(ax1, position = :rc, framevisible = false)

    REw4  = @lift(  en4[key8weno].RPE[1:$iter] .-  en4[key8weno].RPE[1] ) 
    REl4  = @lift( en4[key8leith].RPE[1:$iter] .- en4[key8leith].RPE[1] ) 
    REb4  = @lift(    en4[key8bilap].RPE[1:$iter] .-    en4[key8bilap].RPE[1] ) 
    REw16 = @lift( en16[key8weno].RPE[1:$iter] .- en16[key8weno].RPE[1] ) 
    REl16 = @lift(en16[key8leith].RPE[1:$iter] .-en16[key8leith].RPE[1] ) 
    REb16 = @lift(   en16[key8bilap].RPE[1:$iter] .-   en16[key8bilap].RPE[1] ) 
    RPw4  = @lift(  en4[key8weno].RPE[$iter:$iter] .-  en4[key8weno].RPE[1]) 
    RPl4  = @lift( en4[key8leith].RPE[$iter:$iter] .- en4[key8leith].RPE[1]) 
    RPb4  = @lift(    en4[key8bilap].RPE[$iter:$iter] .-    en4[key8bilap].RPE[1]) 
    RPw16 = @lift( en16[key8weno].RPE[$iter:$iter] .- en16[key8weno].RPE[1]) 
    RPl16 = @lift(en16[key8leith].RPE[$iter:$iter] .-en16[key8leith].RPE[1]) 
    RPb16 = @lift(   en16[key8bilap].RPE[$iter:$iter] .-   en16[key8bilap].RPE[1]) 
        
    plot_stuff!(ax1, REw4, REl4, REb4, REw16, REl16, REb16, 
                     RPw4, RPl4, RPb4, RPw16, RPl16, RPb16, 
                     iter_arr, color1, color2, color3)

    AEw4  = @lift(  en4[key8weno].APE[1:$iter] .-  en4[key8weno].APE[1]) 
    AEl4  = @lift( en4[key8leith].APE[1:$iter] .- en4[key8leith].APE[1]) 
    AEb4  = @lift(    en4[key8bilap].APE[1:$iter] .-    en4[key8bilap].APE[1])
    AEw16 = @lift( en16[key8weno].APE[1:$iter] .- en16[key8weno].APE[1]) 
    AEl16 = @lift(en16[key8leith].APE[1:$iter] .-en16[key8leith].APE[1]) 
    AEb16 = @lift(   en16[key8bilap].APE[1:$iter] .-   en16[key8bilap].APE[1]) 
    APw4  = @lift(  en4[key8weno].APE[$iter:$iter] .-  en4[key8weno].APE[1]) 
    APl4  = @lift( en4[key8leith].APE[$iter:$iter] .- en4[key8leith].APE[1]) 
    APb4  = @lift(    en4[key8bilap].APE[$iter:$iter] .-    en4[key8bilap].APE[1]) 
    APw16 = @lift( en16[key8weno].APE[$iter:$iter] .-  en16[key8weno].APE[1]) 
    APl16 = @lift(en16[key8leith].APE[$iter:$iter] .- en16[key8leith].APE[1]) 
    APb16 = @lift(   en16[key8bilap].APE[$iter:$iter] .-   en16[key8bilap].APE[1])

    plot_stuff!(ax1, AEw4, AEl4, AEb4, AEw16, AEl16, AEb16, 
                     APw4, APl4, APb4, APw16, APl16, APb16, 
                     iter_arr, color1, color2, color3)
    record_video!("totE_video2", figE, iter, 201)

    close(en4)
    close(en16)
end

function plot_enstrophies()
    iter = Observable(201)
    iter_arr = @lift([$iter])

    color1 = :deepskyblue
    color2 = :orange1
    color3 = :firebrick2

    key8weno = "weno5dd"
    key8leith = "lapleith"
    key8bilap = "bilap"
    key16weno = "weno5dd"
    key16leith = "lapleith"
    key16bilap = "bilap"

    @show key8weno

    z4  = jldopen("enstrophies_quarter.jld2")
    z16 = jldopen("enstrophies_sixteen.jld2")

    Zw4  = @lift(  z4[key8weno][1:$iter] ./ 1e6)
    Zl4  = @lift( z4[key8leith][1:$iter] ./ 1e6)
    Zb4  = @lift(    z4[key8bilap][1:$iter] ./ 1e6)
    Zw16 = @lift( z16[key16weno][1:$iter] ./ 1e6)
    Zl16 = @lift(z16[key16leith][1:$iter] ./ 1e6)
    Zb16 = @lift(   z16[key16bilap][1:$iter] ./ 1e6)
    Rw4  = @lift(  z4[key8weno][$iter:$iter] ./ 1e6)
    Rl4  = @lift( z4[key8leith][$iter:$iter] ./ 1e6)
    Rb4  = @lift(    z4[key8bilap][$iter:$iter] ./ 1e6)
    Rw16 = @lift( z16[key16weno][$iter:$iter] ./ 1e6)
    Rl16 = @lift(z16[key16leith][$iter:$iter] ./ 1e6)
    Rb16 = @lift(   z16[key16bilap][$iter:$iter] ./ 1e6)

    figZ = Figure(resolution = (400, 300))
    ax3 = Axis(figZ[1, 1], xgridvisible = false, ygridvisible = false, 
              yticks = ([], []),
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Enstrophy")
        
    plot_stuff!(ax3, Zw4, Zl4, Zb4, Zw16, Zl16, Zb16, 
                     Rw4, Rl4, Rb4, Rw16, Rl16, Rb16, 
                     iter_arr, color1, color2, color3)
    record_video!("Enstrophy_video2", figZ, iter, 201)

    s4  = jldopen("stratif_quarter.jld2")
    s16 = jldopen("stratif_sixteen.jld2")

    Ew4  = @lift(  s4[key8weno][1:$iter] ./ 1e-6)
    El4  = @lift( s4[key8leith][1:$iter] ./ 1e-6)
    Eb4  = @lift(    s4[key8bilap][1:$iter] ./ 1e-6)
    Ew16 = @lift( s16[key16weno][1:$iter] ./ 1e-6)
    El16 = @lift(s16[key16leith][1:$iter] ./ 1e-6)
    Eb16 = @lift(   s16[key16bilap][1:$iter] ./ 1e-6)
    Pw4  = @lift(  s4[key8weno][$iter:$iter] ./ 1e-6)
    Pl4  = @lift( s4[key8leith][$iter:$iter] ./ 1e-6)
    Pb4  = @lift(    s4[key8bilap][$iter:$iter] ./ 1e-6)
    Pw16 = @lift( s16[key16weno][$iter:$iter] ./ 1e-6)
    Pl16 = @lift(s16[key16leith][$iter:$iter] ./ 1e-6)
    Pb16 = @lift(   s16[key16bilap][$iter:$iter] ./ 1e-6)

    figS = Figure(resolution = (400, 300))
    ax4 = Axis(figS[1, 1], xgridvisible = false, ygridvisible = false, 
               xlabel = "days", 
               ylabel = "",
               title = "Mean Stratification",
               yticks = ([0.01, 100], ["0.01", "100"]))
        
    plot_stuff!(ax4, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    axislegend(ax4, position = :lt, framevisible = false)
    record_video!("Stratif_video2", figS, iter, 201)
end
