using FFTW
using Oceananigans.Grids: φnode

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

import Base

Base.:(+)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .+ t.spec, s.freq)
Base.:(*)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .* t.spec, s.freq)
Base.:(/)(s::Spectrum, t::Int)      = Spectrum(s.spec ./ t, s.freq)

Base.real(s::Spectrum) = Spectrum(real.(s.spec), s.freq)
Base.abs(s::Spectrum) = Spectrum(abs.(s.spec), s.freq)

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

@kernel function _compute_zonal_spectra!(Uspec, Vspec, Bspec, Ωspec, WBspec, STspec, PVspec, grid, u, v, ζ, w, b, st, pv)
    j, k = @index(Global, NTuple)

    Nx = size(grid, 1)

    Uspec[j, k]  = power_spectrum_1d_x(Array(interior(u,  :, j, k)), Array(grid.λᶠᵃᵃ.parent)[1:Nx])
    Vspec[j, k]  = power_spectrum_1d_x(Array(interior(v,  :, j, k)), Array(grid.λᶜᵃᵃ.parent)[1:Nx])
    Bspec[j, k]  = power_spectrum_1d_x(Array(interior(b,  :, j, k)), Array(grid.λᶜᵃᵃ.parent)[1:Nx])
    Ωspec[j, k]  = power_spectrum_1d_x(Array(interior(ζ,  :, j, k)), Array(grid.λᶠᵃᵃ.parent)[1:Nx])
    STspec[j, k] = power_spectrum_1d_x(Array(interior(st, :, j, k)), Array(grid.λᶠᵃᵃ.parent)[1:Nx])
    PVspec[j, k] = power_spectrum_1d_x(Array(interior(pv, :, j, k)), Array(grid.λᶠᵃᵃ.parent)[1:Nx])
    WBspec[j, k] = power_cospectrum_1d(Array(interior(w,  :, j, k)), Array(interior(b, :, j, k)), Array(grid.λᶜᵃᵃ.parent)[1:Nx])
end

@kernel function _update_spectra!(Ufinal, Vfinal, Bfinal, Ωfinal, WBfinal, STfinal, PVfinal, Uspec, Vspec, Bspec, Ωspec, WBspec, STspec, PVspec)
    j, k = @index(Global, NTuple)

    Ufinal[j, k].spec  .+= Uspec[j, k].spec
    Vfinal[j, k].spec  .+= Vspec[j, k].spec
    Bfinal[j, k].spec  .+= Bspec[j, k].spec
    Ωfinal[j, k].spec  .+= Ωspec[j, k].spec
    WBfinal[j, k].spec .+= WBspec[j, k].spec
    STfinal[j, k].spec .+= STspec[j, k].spec
    PVfinal[j, k].spec .+= PVspec[j, k].spec
end

@inline function _stretching_term(i, j, k, grid, ∂zb, N²)
    φ = φnode(i, j, k, grid, Center(), Center(), Face())
    f = 2 * 7.292115e-5 * sind(φ)
    return ∂zb[i, j, k] * f / max(1e-10, N²[i, j, k])
end

function compute_spectra(f::Dict, time)
    grid = f[:u].grid

    Nx, Ny, Nz = size(grid)

    U  = Array{Spectrum}(undef, Ny, Nz)
    V  = Array{Spectrum}(undef, Ny, Nz)
    B  = Array{Spectrum}(undef, Ny, Nz)
    Ω  = Array{Spectrum}(undef, Ny, Nz)
    WB = Array{Spectrum}(undef, Ny, Nz)
    ST = Array{Spectrum}(undef, Ny, Nz)
    PV = Array{Spectrum}(undef, Ny, Nz)

    Uspec  = Array{Spectrum}(undef, Ny, Nz)
    Vspec  = Array{Spectrum}(undef, Ny, Nz)
    Bspec  = Array{Spectrum}(undef, Ny, Nz)
    Ωspec  = Array{Spectrum}(undef, Ny, Nz)
    WBspec = Array{Spectrum}(undef, Ny, Nz)
    STspec = Array{Spectrum}(undef, Ny, Nz)
    PVspec = Array{Spectrum}(undef, Ny, Nz)

    u = Field{Face, Center, Center}(grid)
    v = Field{Center, Face, Center}(grid)

    b = Field{Center, Center, Center}(grid)
    set!(u, f[:u][time[1]])
    set!(v, f[:v][time[1]])
    set!(b, f[:b][time[1]])
    fill_halo_regions!((u, v, b))

    w, b = f[:w][time[1]], f[:b][time[1]]

    ζ = compute!(Field(VerticalVorticityOperation((; u, v))))

    ∂zb = compute!(Field(∂z(b)))
    N²  = compute!(Field(Average(∂zb, dims = 1)))

    st_op = KernelFunctionOperation{Center, Center, Face}(_stretching_term, grid, ∂zb, N²)
    st    = compute!(Field(st_op))
    pv    = compute!(Field(st + ζ))

    Nx, Ny, Nz = size(grid)
    launch!(CPU(), grid, (Ny, Nz), _compute_zonal_spectra!, U, V, B, Ω, WB, ST, PV, grid, u, v, ζ, w, b, st, pv)

    @show length(time)
    if length(time) > 1
        for t in time[2:end]
            @info "doing time $time"
            set!(u, f[:u][t])
            set!(v, f[:v][t])
            fill_halo_regions!((u, v))
            ζ = compute!(Field(VerticalVorticityOperation((; u, v))))

            w, b = f[:w][t], f[:b][t]

            launch!(CPU(), grid, (Ny, Nz), _compute_zonal_spectra!, Uspec, Vspec, Bspec, Ωspec, WBspec, STspec, PVspec, grid, u, v, ζ, w, b, st, pv)
            launch!(CPU(), grid, (Ny, Nz), _update_spectra!, U, V, B, Ω, WB, ST, PV, Uspec, Vspec, Bspec, Ωspec, WBspec, STspec, PVspec)
        end
    end

    @info "finished"
    
    return (; U, V, B, Ω, WB, ST, PV)
end
