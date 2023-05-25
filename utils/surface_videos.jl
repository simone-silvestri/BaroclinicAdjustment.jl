using GLMakie
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: VerticalVorticity, KineticEnergy, DeformationRadius

using GeoMakie

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

function geographic2cartesian(λ, φ, r=1)
    Nλ = length(λ)
    Nφ = length(φ)

    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ) 
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

function surface_videos(trailing_character = "_weaker")

    file_prefix = ["bilap", "leith", "lapleith", "smag", "qgleith",
                   "weno5v", "weno5d", "weno7v", "weno7d", "weno9v", "weno9d", "wenoHv", "weno9MD",
                   "weno5F", "weno7F", "weno9F"]

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    iter = Observable(1)

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)
            f = all_fieldtimeseries(filename)

            grid = f[:u].grid
            x, y, z = nodes(grid, Face(), Face(), Center())
            N = min(size(x), size(y))[1]
            x = x[1:N]
            y = y[1:N]

            Nt = length(f[:b].times)
            @info "vorticity video $filename"    
            ζ = @lift(interior(VerticalVorticity(f, $iter), :, :, 50))

            fig = Figure()
            ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            surface!(ax, x, y, ζ, colorrange = (-1e-5, 1e-5), shading = false)
            hidedecorations!(ax)
            hidespines!(ax)
        
            record_video!("vort_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
        
            ζ = @lift(interior(VerticalVorticity(f, $iter), :, :, 30))

            fig = Figure()
            ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            surface!(ax, x, y, ζ, colorrange = (-1e-5, 1e-5), shading = false)
            hidedecorations!(ax)
            hidespines!(ax)
        
            record_video!("lowvort_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            @info "z-vorticity video $filename"    
            ζ = @lift(interior(VerticalVorticity(f, $iter), 160, :, :))
            
            fig = Figure()
            ax = Axis(fig[1, 1])            # 
            heatmap!(ax, ζ, colorrange = (-1e-5, 1e-5))
            hidedecorations!(ax)
            hidespines!(ax)

            record_video!("zvort_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            # @info "kinetic energy video $filename"    
            # E = @lift(interior(KineticEnergy(f, $iter), :, :, 50))

            # x, y, z = nodes(grid, Center(), Center(), Center())
            # x = x[1:N+1]
            # y = y[1:N+1]
            # fig = Figure()
            # ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            # surface(ax, x, y, E, colorrange = (0, 1.0), colormap = :solar)
            # hidedecorations!(ax)
            # hidespines!(ax)
            # record_video!("energy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            @info "buoyancy video $filename"    
            b = @lift(interior(f[:b][$iter], :, :, 50))

            fig = Figure()
            ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            surface!(ax, x, y, b, colorrange = (0, 0.006), colormap = :thermometer, shading = false)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("buoyancy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

        end
    end
end


function surface_videos_rectilinear(trailing_character = "_weaker")
    # file_prefix = ["leith", "lapleith", "bilap", "weno5vd",
    #                "smag", "weno5dd", "weno9", "weno9dd", "highres"]

    file_prefix = ["weno9v4_rect"] #, "lapleith"]

    @show file_prefix
    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)

    iter = Observable(1)

    for (prefix, filename) in zip(file_prefix, filenames)
        if isfile(filename)
            f = all_fieldtimeseries(filename)

            grid = f[:u].grid
            x, y, z = nodes(grid, Face(), Face(), Center())
            N = min(size(x), size(y))[1]
            x = x[1:N]
            y = y[1:N]

            Nt = length(f[:b].times)
            @info "vorticity video $filename"    
            ζ = @lift(interior(VerticalVorticity(f, $iter), :, :, 50))

            fig = Figure()
            ax = Axis(fig[1, 1])
            heatmap!(ax, x, y, ζ, colorrange = (-1e-5, 1e-5))
            hidedecorations!(ax)
            hidespines!(ax)
        
            record_video!("vort_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
        
            @info "kinetic energy video $filename"    
            E = @lift(interior(KineticEnergy(f, $iter), :, :, 50))

            x, y, z = nodes(grid, Center(), Center(), Center())
            x = x[1:N+1]
            y = y[1:N+1]
            # fig = Figure()
            # ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            # surface(ax, x, y, E, colorrange = (0, 1.0), colormap = :solar)
            # hidedecorations!(ax)
            # hidespines!(ax)
            # record_video!("energy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            @info "buoyancy video $filename"    
            b = @lift(interior(f[:b][$iter], :, :, 50))

            fig = Figure()
            ax = Axis(fig[1, 1])
            heatmap!(ax, x, y, b, colorrange = (0, 0.006), colormap = :thermometer)
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

    en4  = jldopen("energies_eight.jld2")
    en16 = jldopen("energies_sixteen.jld2")

    Ew4  = @lift(  en4["weno5dd"].KE[1:$iter] .* 1e3)
    El4  = @lift( en4["lapleith"].KE[1:$iter] .* 1e3)
    Eb4  = @lift(    en4["bilap"].KE[1:$iter] .* 1e3)
    Ew16 = @lift( en16["weno5dd"].KE[1:$iter] .* 1e3)
    El16 = @lift(en16["lapleith"].KE[1:$iter] .* 1e3)
    Eb16 = @lift(   en16["bilap"].KE[1:$iter] .* 1e3)
    Pw4  = @lift(  en4["weno5dd"].KE[$iter:$iter] .* 1e3)
    Pl4  = @lift( en4["lapleith"].KE[$iter:$iter] .* 1e3)
    Pb4  = @lift(    en4["bilap"].KE[$iter:$iter] .* 1e3)
    Pw16 = @lift( en16["weno5dd"].KE[$iter:$iter] .* 1e3)
    Pl16 = @lift(en16["lapleith"].KE[$iter:$iter] .* 1e3)
    Pb16 = @lift(   en16["bilap"].KE[$iter:$iter] .* 1e3)

    figE = Figure(resolution = (400, 300))
    ax1 = Axis(figE[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Kinetic Energy")
        
    plot_stuff!(ax1, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                    Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                    iter_arr, color1, color2, color3)
    axislegend(ax1, position = :rc, framevisible = false)
    record_video!("energy_video2", figE, iter, 201)

    Ew4  = @lift(  en4["weno5dd"].RPE[1:$iter] .-  en4["weno5dd"].RPE[1] ) # ./ 1.538 ./ 1e21)
    El4  = @lift( en4["lapleith"].RPE[1:$iter] .- en4["lapleith"].RPE[1] ) # ./ 1.538 ./ 1e21)
    Eb4  = @lift(    en4["bilap"].RPE[1:$iter] .-    en4["bilap"].RPE[1] ) # ./ 1.538 ./ 1e21)
    Ew16 = @lift( en16["weno5dd"].RPE[1:$iter] .- en16["weno5dd"].RPE[1] ) # ./ 1.538 ./ 1e21)
    El16 = @lift(en16["lapleith"].RPE[1:$iter] .-en16["lapleith"].RPE[1] ) # ./ 1.538 ./ 1e21)
    Eb16 = @lift(   en16["bilap"].RPE[1:$iter] .-   en16["bilap"].RPE[1] ) # ./ 1.538 ./ 1e21)
    Pw4  = @lift(  en4["weno5dd"].RPE[$iter:$iter] .-  en4["weno5dd"].RPE[1]) # ./ 1.538 ./ 1e21)
    Pl4  = @lift( en4["lapleith"].RPE[$iter:$iter] .- en4["lapleith"].RPE[1]) # ./ 1.538 ./ 1e21)
    Pb4  = @lift(    en4["bilap"].RPE[$iter:$iter] .-    en4["bilap"].RPE[1]) # ./ 1.538 ./ 1e21)
    Pw16 = @lift( en16["weno5dd"].RPE[$iter:$iter] .- en16["weno5dd"].RPE[1]) # ./ 1.538 ./ 1e21)
    Pl16 = @lift(en16["lapleith"].RPE[$iter:$iter] .-en16["lapleith"].RPE[1]) # ./ 1.538 ./ 1e21)
    Pb16 = @lift(   en16["bilap"].RPE[$iter:$iter] .-   en16["bilap"].RPE[1]) # ./ 1.538 ./ 1e21)

    figR = Figure(resolution = (400, 300))
    ax2 = Axis(figR[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Reference PE")
        
    plot_stuff!(ax2, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    record_video!("RPE_video2", figR, iter, 201)

    Ew4  = @lift(  en4["weno5dd"].APE[1:$iter] .-  en4["weno5dd"].APE[1]) # ./ 1.538 ./ 1e21)
    El4  = @lift( en4["lapleith"].APE[1:$iter] .- en4["lapleith"].APE[1]) # ./ 1.538 ./ 1e21)
    Eb4  = @lift(    en4["bilap"].APE[1:$iter] .-    en4["bilap"].APE[1]) # ./ 1.538 ./ 1e21)
    Ew16 = @lift( en16["weno5dd"].APE[1:$iter] .- en16["weno5dd"].APE[1]) # ./ 1.538 ./ 1e21)
    El16 = @lift(en16["lapleith"].APE[1:$iter] .-en16["lapleith"].APE[1]) # ./ 1.538 ./ 1e21)
    Eb16 = @lift(   en16["bilap"].APE[1:$iter] .-   en16["bilap"].APE[1]) # ./ 1.538 ./ 1e21)
    Pw4  = @lift(  en4["weno5dd"].APE[$iter:$iter] .-  en4["weno5dd"].APE[1]) # ./ 1.538 ./ 1e21)
    Pl4  = @lift( en4["lapleith"].APE[$iter:$iter] .- en4["lapleith"].APE[1]) # ./ 1.538 ./ 1e21)
    Pb4  = @lift(    en4["bilap"].APE[$iter:$iter] .-    en4["bilap"].APE[1]) # ./ 1.538 ./ 1e21)
    Pw16 = @lift( en16["weno5dd"].APE[$iter:$iter] .- en16["weno5dd"].APE[1]) # ./ 1.538 ./ 1e21)
    Pl16 = @lift(en16["lapleith"].APE[$iter:$iter] .-en16["lapleith"].APE[1]) # ./ 1.538 ./ 1e21)
    Pb16 = @lift(   en16["bilap"].APE[$iter:$iter] .-   en16["bilap"].APE[1]) # ./ 1.538 ./ 1e21)

    figR = Figure(resolution = (400, 300))
    ax2 = Axis(figR[1, 1], xgridvisible = false, ygridvisible = false, 
              xlabel = "days", 
              ylabel = "",
              title = "Integrated Available PE")
        
    plot_stuff!(ax2, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    record_video!("APE_video2", figR, iter, 201)

    close(en4)
    close(en16)

    z4  = jldopen("enstrophies_eight.jld2")
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
    record_video!("Enstrophy_video2", figZ, iter, 201)

    s4  = jldopen("stratif_eight.jld2")
    s16 = jldopen("stratif_sixteen.jld2")

    Ew4  = @lift(  s4["weno5dd"][1:$iter] ./ 1e-6)
    El4  = @lift( s4["lapleith"][1:$iter] ./ 1e-6)
    Eb4  = @lift(    s4["bilap"][1:$iter] ./ 1e-6)
    Ew16 = @lift( s16["weno5dd"][1:$iter] ./ 1e-6)
    El16 = @lift(s16["lapleith"][1:$iter] ./ 1e-6)
    Eb16 = @lift(   s16["bilap"][1:$iter] ./ 1e-6)
    Pw4  = @lift(  s4["weno5dd"][$iter:$iter] ./ 1e-6)
    Pl4  = @lift( s4["lapleith"][$iter:$iter] ./ 1e-6)
    Pb4  = @lift(    s4["bilap"][$iter:$iter] ./ 1e-6)
    Pw16 = @lift( s16["weno5dd"][$iter:$iter] ./ 1e-6)
    Pl16 = @lift(s16["lapleith"][$iter:$iter] ./ 1e-6)
    Pb16 = @lift(   s16["bilap"][$iter:$iter] ./ 1e-6)

    figS = Figure(resolution = (400, 300))
    ax4 = Axis(figS[1, 1], xgridvisible = false, ygridvisible = false, 
               xlabel = "days", 
               ylabel = "",
               title = "Mean Stratification",
               yticks = ([2, 4, 8], ["2.0", "4.0", "8.0"]))
        
    plot_stuff!(ax4, Ew4, El4, Eb4, Ew16, El16, Eb16, 
                     Pw4, Pl4, Pb4, Pw16, Pl16, Pb16, 
                     iter_arr, color1, color2, color3)
    axislegend(ax4, position = :lt, framevisible = false)
    record_video!("Stratif_video2", figS, iter, 201)
end

function geographic2cartesian(λ, φ, r=1)
    Nλ = length(λ)
    Nφ = length(φ)

    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ) 
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

function video_highres()

    file = all_fieldtimeseries("highres_final_snapshots.jld2")
    # Domain
    Lz = 1kilometers     # depth [m]
    Ny = 1000
    Nz = 50


    iter = Observable(1)

    grid = LatitudeLongitudeGrid(topology = (Periodic, Bounded, Bounded),
                                 size = (Ny, Ny, Nz), 
                                 longitude = (-10, 10),
                                 latitude = (-60, -40),
                                 z = (-Lz, 0))
    
    slice_north  = @lift(Array(interior(file[:b][$iter], :, grid.Ny, :)))
    slice_south  = @lift(Array(interior(file[:b][$iter], :, 1, :)))
    slice_top    = @lift(Array(interior(file[:b][$iter], :, :, grid.Nz)))
    slice_bottom = @lift(Array(interior(file[:b][$iter], :, :, 1)))
    slice_east   = @lift(Array(interior(file[:b][$iter], grid.Nx, :, :)))
    slice_west   = @lift(Array(interior(file[:b][$iter], 1, :, :)))


    xb, yb, zb = nodes(grid, Center(), Center(), Center())
    xb = (xb .- xb[1]) .* 200
    yb = (yb .- yb[1]) .* 200

    fig = Figure()
    ax = fig[1:5, 1:5] = LScene(fig, show_axis=false)

    surface!(ax, yb, zb, slice_west;   transformation = (:yz, xb[1]),   colormap=:thermometer, colorrange = (0, 0.006))
    surface!(ax, yb, zb, slice_east;   transformation = (:yz, xb[end]), colormap=:thermometer, colorrange = (0, 0.006))
    surface!(ax, xb, zb, slice_south;  transformation = (:xz, yb[1]),   colormap=:thermometer, colorrange = (0, 0.006))
    surface!(ax, xb, zb, slice_north;  transformation = (:xz, yb[end]), colormap=:thermometer, colorrange = (0, 0.006))
    surface!(ax, xb, yb, slice_bottom; transformation = (:xy, zb[1]),   colormap=:thermometer, colorrange = (0, 0.006))
    surface!(ax, xb, yb, slice_top;    transformation = (:xy, zb[end]), colormap=:thermometer, colorrange = (0, 0.006))

    rotate_cam!(ax.scene, (π/18, -π/6, 0))

    record_video!("final_video", fig, iter, 201)

    return fig
end
