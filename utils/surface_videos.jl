using GLMakie
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: VerticalVorticity, KineticEnergy, DeformationRadius

using Oceananigans
using JLD2

add_trailing_name(name) = name * "_snapshots.jld2"

function record_video!(name, fig, iter, Nt) 
    GLMakie.record(fig, name * ".mp4", 1:Nt, framerate = 11) do i
        @info "step $i"; 
        iter[] = i; 
    end
end

function surface_videos(trailing_character = "_weaker")
    file_prefix = ["weno5vd", "leith", "lapleith", "bilap",
                   "smag", "weno5dd", "weno9", "weno9dd", "highres"]

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

            @info "deformation radius video $filename"    
            R = @lift(interior(DeformationRadius(f, $iter), :, :, 1))

            fig = Figure()
            ax  = Axis(fig[1, 1])
            heatmap!(ax, R, colorrange = (1e1, 1e4), colormap = :thermometer)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("Lr_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
        end
    end
end

function plot_stuff!(ax, Ew4, El4, Eb4, Ew16, El16, Eb16, Pw4, 
                     Pl4, Pb4, Pw16, Pl16, Pb16, iter_arr,
                     color1, color2, color3)
    lines!(ax, Ew4, linewidth = 2, color = color3, label = "3D-WENO-VI")
    lines!(ax, El4, linewidth = 2, color = color1, label = "Leith")
    lines!(ax, Eb4, linewidth = 2, color = color2, label = "Biharmonic")
    lines!(ax, Ew16, linewidth = 2, linestyle = :dash, color = color3)
    lines!(ax, El16, linewidth = 2, linestyle = :dash, color = color1)
    lines!(ax, Eb16, linewidth = 2, linestyle = :dash, color = color2)

    scatter!(ax, iter_arr, Pw4, color = color3, markersize = 12)
    scatter!(ax, iter_arr, Pl4, color = color1, markersize = 12)
    scatter!(ax, iter_arr, Pb4, color = color2, markersize = 12)
    scatter!(ax, iter_arr, Pw16, color = color3, markersize = 12)
    scatter!(ax, iter_arr, Pl16, color = color1, markersize = 12)
    scatter!(ax, iter_arr, Pb16, color = color2, markersize = 12)
end

function plot_all()
    iter = Observable(201)
    iter_arr = @lift([$iter])

    color1 = :deepskyblue
    color2 = :orange1
    color3 = :firebrick2

    en4  = jldopen("energies_quarter.jld2")
    en16 = jldopen("energies_sixteen.jld2")

    Ew4  = @lift(  en4["weno5dd"].KE[1:$iter] ./ 1e15)
    El4  = @lift( en4["lapleith"].KE[1:$iter] ./ 1e15)
    Eb4  = @lift(    en4["bilap"].KE[1:$iter] ./ 1e15)
    Ew16 = @lift( en16["weno5dd"].KE[1:$iter] ./ 1e15)
    El16 = @lift(en16["lapleith"].KE[1:$iter] ./ 1e15)
    Eb16 = @lift(   en16["bilap"].KE[1:$iter] ./ 1e15)
    Pw4  = @lift(  en4["weno5dd"].KE[$iter:$iter] ./ 1e15)
    Pl4  = @lift( en4["lapleith"].KE[$iter:$iter] ./ 1e15)
    Pb4  = @lift(    en4["bilap"].KE[$iter:$iter] ./ 1e15)
    Pw16 = @lift( en16["weno5dd"].KE[$iter:$iter] ./ 1e15)
    Pl16 = @lift(en16["lapleith"].KE[$iter:$iter] ./ 1e15)
    Pb16 = @lift(   en16["bilap"].KE[$iter:$iter] ./ 1e15)

    figE = Figure(resolution = (400, 300))
    ax1 = Axis(figE[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Kinetic Energy")
        
    plot_stuff!(ax1, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                    Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                    iter_arr, color1, color2, color3)
    axislegend(ax1, position = :rc, framevisible = false)
    record_video!("energy_video", figE, iter, 201)

    Ew4  = @lift(  en4["weno5dd"].RPE[1:$iter])
    El4  = @lift( en4["lapleith"].RPE[1:$iter])
    Eb4  = @lift(    en4["bilap"].RPE[1:$iter])
    Ew16 = @lift( en16["weno5dd"].RPE[1:$iter])
    El16 = @lift(en16["lapleith"].RPE[1:$iter])
    Eb16 = @lift(   en16["bilap"].RPE[1:$iter])
    Pw4  = @lift(  en4["weno5dd"].RPE[$iter:$iter])
    Pl4  = @lift( en4["lapleith"].RPE[$iter:$iter])
    Pb4  = @lift(    en4["bilap"].RPE[$iter:$iter])
    Pw16 = @lift( en16["weno5dd"].RPE[$iter:$iter])
    Pl16 = @lift(en16["lapleith"].RPE[$iter:$iter])
    Pb16 = @lift(   en16["bilap"].RPE[$iter:$iter])

    figR = Figure(resolution = (400, 300))
    ax2 = Axis(figR[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Reference PE")
        
    plot_stuff!(ax2, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    record_video!("RPE_video", figR, iter, 201)

    close(en4)
    close(en16)

    z4  = jldopen("enstrophies_quarter.jld2")
    z16 = jldopen("enstrophies_sixteen.jld2")

    Zw4  = @lift(  z4["weno5dd"][1:$iter] ./ 1e6)
    Zl4  = @lift( z4["lapleith"][1:$iter] ./ 1e6)
    Zb4  = @lift(    z4["bilap"][1:$iter] ./ 1e6)
    Zw16 = @lift( z16["weno5dd"][1:$iter] ./ 1e6)
    Zl16 = @lift(z16["lapleith"][1:$iter] ./ 1e6)
    Zb16 = @lift(   z16["bilap"][1:$iter] ./ 1e6)
    Rw4  = @lift(  z4["weno5dd"][$iter:$iter] ./ 1e6)
    Rl4  = @lift( z4["lapleith"][$iter:$iter] ./ 1e6)
    Rb4  = @lift(    z4["bilap"][$iter:$iter] ./ 1e6)
    Rw16 = @lift( z16["weno5dd"][$iter:$iter] ./ 1e6)
    Rl16 = @lift(z16["lapleith"][$iter:$iter] ./ 1e6)
    Rb16 = @lift(   z16["bilap"][$iter:$iter] ./ 1e6)

    figZ = Figure(resolution = (400, 300))
    ax3 = Axis(figZ[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Enstrophy")
        
    plot_stuff!(ax3, Zw4, Zl4, Zb4, Zw16, Zl16, Zb16, 
                     Rw4, Rl4, Rb4, Rw16, Rl16, Rb16, 
                     iter_arr, color1, color2, color3)
    # axislegend(ax3, position = :rt, framevisible = false)
    record_video!("Enstrophy_video", figZ, iter, 201)

    s4  = jldopen("stratif_quarter.jld2")
    s16 = jldopen("stratif_sixteen.jld2")

    Ew4  = @lift(  s4["weno5dd"][1:$iter])
    El4  = @lift( s4["lapleith"][1:$iter])
    Eb4  = @lift(    s4["bilap"][1:$iter])
    Ew16 = @lift( s16["weno5dd"][1:$iter])
    El16 = @lift(s16["lapleith"][1:$iter])
    Eb16 = @lift(   s16["bilap"][1:$iter])
    Pw4  = @lift(  s4["weno5dd"][$iter:$iter])
    Pl4  = @lift( s4["lapleith"][$iter:$iter])
    Pb4  = @lift(    s4["bilap"][$iter:$iter])
    Pw16 = @lift( s16["weno5dd"][$iter:$iter])
    Pl16 = @lift(s16["lapleith"][$iter:$iter])
    Pb16 = @lift(   s16["bilap"][$iter:$iter])

    figS = Figure(resolution = (400, 300))
    ax4 = Axis(figS[1, 1], xgridvisible = false, ygridvisible = false, 
               xlabel = "days", 
               ylabel = "",
               title = "Mean Stratification")
        
    plot_stuff!(ax4, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    axislegend(ax4, position = :lb, framevisible = false)
    record_video!("Stratif_video", figS, iter, 201)
end