using FFTW

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

@inline onefunc(args...)  = 1.0
@inline hann_window(n, N) = sin(π * n / N)^2 

function average_spectra(var::FieldTimeSeries, xlim, ylim; k = 69, spectra = power_spectrum_1d_x, windowing = onefunc)

    xdomain = xnodes(var[1])[xlim]
    ydomain = ynodes(var[1])[ylim]

    Nt = length(var.times)

    spec = spectra(interior(var[1], xlim, ylim, k), xdomain, ydomain; windowing) 

    for i in 2:Nt
        spec.spec .+= spectra(interior(var[i], xlim, ylim, k), xdomain, ydomain).spec 
    end

    spec.spec ./= Nt

    return spec
end

function power_cospectrum_1d(var1, var2, x; windowing = onefunc)

    Nx = length(x)
    Nfx = Int64(Nx)
    
    spectra = zeros(ComplexF64, Int(Nfx/2))
    
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0, +ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    windowed_var1 = [var1[i] * windowing(i, Nfx) for i in 1:Nfx]
    windowed_var2 = [var2[i] * windowing(i, Nfx) for i in 1:Nfx]
    fourier1      = fft(windowed_var1) / Nfx
    fourier2      = fft(windowed_var2) / Nfx
    spectra[1]    += fourier1[1] .* conj(fourier2[1]) .+ fourier2[1] .* conj(fourier1[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += fourier1[m] .* conj(fourier2[m]) .+ fourier2[m] .* conj(fourier1[m])
    end
    return Spectrum(spectra, freqs)
end

function power_spectrum_1d_x(var, x; windowing = onefunc)

    Nx = length(x)
    Nfx = Int64(Nx)
    
    spectra = zeros(ComplexF64, Int(Nfx/2))
    
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    windowed_var = [var[i] * windowing(i, Nfx) for i in 1:Nfx]
    fourier      = fft(windowed_var) / Nfx
    spectra[1]  += fourier[1] .* conj(fourier[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += 2.0 * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
    end
    return Spectrum(spectra, freqs)
end

@kernel function _compute_zonal_spectra!(Uspec, Vspec, Ωspec, WBspec, grid, u, v, ζ, w, b)
    j, k = @index(Global, NTuple)

    Nx = size(grid, 1)

    Uspec[j, k]  = power_spectrum_1d_x(interior(u, :, j, k), grid.λᶠᵃᵃ[1:Nx])
    Vspec[j, k]  = power_spectrum_1d_x(interior(v, :, j, k), grid.λᶜᵃᵃ[1:Nx])
    Ωspec[j, k]  = power_spectrum_1d_x(interior(ζ, :, j, k), grid.λᶠᵃᵃ[1:Nx])
    WBspec[j, k] = power_cospectrum_1d(interior(w, :, j, k), interior(b, :, j, k), grid.λᶜᵃᵃ[1:Nx])
end

@kernel function _update_spectra!(Ufinal, Vfinal, Ωfinal, WBfinal, Uspec, Vspec, Ωspec, WBspec)
    j, k = @index(Global, NTuple)

    Ufinal[j, k].spec  .+= Uspec[j, k].spec
    Vfinal[j, k].spec  .+= Vspec[j, k].spec
    Ωfinal[j, k].spec  .+= Ωspec[j, k].spec
    WBfinal[j, k].spec .+= WBspec[j, k].spec
end

function compute_spectra(f::Dict, time)
    grid = f[:u].grid

    Nx, Ny, Nz = size(grid)

    U  = Array{Spectrum}(undef, Ny, Nz)
    V  = Array{Spectrum}(undef, Ny, Nz)
    Ω  = Array{Spectrum}(undef, Ny, Nz)
    WB = Array{Spectrum}(undef, Ny, Nz)

    Uspec  = Array{Spectrum}(undef, Ny, Nz)
    Vspec  = Array{Spectrum}(undef, Ny, Nz)
    Ωspec  = Array{Spectrum}(undef, Ny, Nz)
    WBspec = Array{Spectrum}(undef, Ny, Nz)

    u = Field{Face, Center, Center}(grid)
    v = Field{Center, Face, Center}(grid)

    set!(u, f[:u][time[1]])
    set!(v, f[:v][time[1]])
    fill_halo_regions!((u, v))

    w, b = f[:w][time[1]], f[:b][time[1]]

    ζ = compute!(Field(VerticalVorticityOperation((; u, v))))
    Nx, Ny, Nz = size(grid)
    launch!(CPU(), grid, (Ny, Nz), _compute_zonal_spectra!, U, V, Ω, WB, grid, u, v, ζ, w, b)

    @show length(time)
    if length(time) > 1
        for t in time[2:end]
            @info "doing time $time"
            set!(u, f[:u][t])
            set!(v, f[:v][t])
            fill_halo_regions!((u, v))
            ζ = compute!(Field(VerticalVorticityOperation((; u, v))))

            w, b = f[:w][t], f[:b][t]

            launch!(CPU(), grid, (Ny, Nz), _compute_zonal_spectra!, Uspec, Vspec, Ωspec, WBspec, grid, u, v, ζ, w, b)
            launch!(CPU(), grid, (Ny, Nz), _update_spectra!, U, V, Ω, WB, Uspec, Vspec, Ωspec, WBspec)
        end
    end
    
    return (; U, V, Ω, WB)
end
