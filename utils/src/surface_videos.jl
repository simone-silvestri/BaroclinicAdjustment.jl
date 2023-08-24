using CairoMakie
using BaroclinicAdjustment
using BaroclinicAdjustment: add_trailing_characters
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: VerticalVorticity, KineticEnergy, DeformationRadius

using GeoMakie

using Oceananigans
using Oceananigans.Units
using JLD2

function record_video!(name, fig, iter, Nt) 
    CairoMakie.record(fig, name * ".mp4", 1:Nt, framerate = 11) do i
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

function surface_videos(trailing_character = "_weaker", file_prefix = generate_names())

    filenames = add_trailing_characters.(file_prefix, trailing_character)
    filenames = add_trailing_name.(filenames)
    @show filenames

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

            @info "buoyancy video $filename"    
            b = @lift(interior(f[:b][$iter], :, :, 50))

            fig = Figure()
            ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = false, lonlims = (-10, 10), latlims = (-60, -40))            # 
            surface!(ax, x, y, b, colorrange = (0, 0.006), colormap = :thermometer, shading = false)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("buoyancy_" * prefix * trailing_character * ".jld2", fig, iter, Nt)

            @info "mean buoyancy video $filename"    
            B = @lift(mean(interior(f[:b][$iter], :, :, :), dims = 1)[1, :, :])

            fig = Figure()
            ax = Axis(fig[1, 1])            
            contourf!(ax, B, levels = range(-0.004, 0.006, length = 10), colormap = :thermometer)
            hidedecorations!(ax)
            hidespines!(ax)
            record_video!("buoyancy_mean_" * prefix * trailing_character * ".jld2", fig, iter, Nt)
        end
    end
end

function surface_videos_rectilinear(trailing_character = "_weaker", file_prefix = generate_names())
    
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
