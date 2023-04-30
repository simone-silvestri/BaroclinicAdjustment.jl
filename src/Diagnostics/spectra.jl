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

function power_spectrum_1d_x(var, x, y; windowing = onefunc)
    Nx  = length(x)
    Ny  = length(y)
    Nfx = Int64(Nx)
    
    spectra = zeros(Float64, Int(Nfx/2))
    
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    for j in 1:Ny
        windowed_var = [var[i, j] * windowing(i, Nfx) for i in 1:Nfx]
        fourier      = fft(windowed_var) / Nfx
        spectra[1]  += fourier[1] .* conj(fourier[1])

        for m in 2:Int(Nfx/2)
            spectra[m] += 2.0 * fourier[m] * conj(fourier[m])  / Ny # factor 2 for neg freq contribution
        end
    end
    return Spectrum(spectra, freqs)
end

function power_spectrum_1d_x(var, x; windowing = onefunc)
    Nx  = length(x)
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

@kernel function _compute_zonal_spectra!(Uspec, Vspec, Ωspec, grid, time, u, v, ζ)
    j, k = @index(Global, NTuple)

    Nx = size(grid, 1)

    Uspec[time, j, k] = power_spectrum_1d_x(interior(u, :, j, k), grid.λᶠᵃᵃ[1:Nx])
    Vspec[time, j, k] = power_spectrum_1d_x(interior(v, :, j, k), grid.λᶜᵃᵃ[1:Nx])
    Ωspec[time, j, k] = power_spectrum_1d_x(interior(ζ, :, j, k), grid.λᶠᵃᵃ[1:Nx])
end

function compute_spectra(f::Dict)
    grid = f[:u].grid

    Nx, Ny, Nz = size(grid)
    Nt = length(f[:u].times)

    Uspec = Array{Spectrum}(undef, Nt, Ny, Nz)
    Vspec = Array{Spectrum}(undef, Nt, Ny, Nz)
    Ωspec = Array{Spectrum}(undef, Nt, Ny, Nz)

    for time in 1:length(f[:u].times)
        @info "doing time $time"
        u = f[:u][time]
        v = f[:v][time]
        ζ = compute!(Field(VerticalVorticityOperation(f, time)))

        launch!(CPU(), grid, :yz, _compute_zonal_spectra!, Uspec, Vspec, Ωspec, grid, time, u, v, ζ)
    end

    return (; Uspec, Vspec, Ωspec)
end