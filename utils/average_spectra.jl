using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using BaroclinicAdjustment.Diagnostics: Spectrum
using KernelAbstractions: @index, @kernel
# using GLMakie
using JLD2

using Statistics

import Base

Base.:(+)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .+ t.spec, s.freq)
Base.:(*)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .* t.spec, s.freq)
Base.:(/)(s::Spectrum, t::Int)      = Spectrum(s.spec ./ t, s.freq)

Base.real(s::Spectrum) = Spectrum(real.(s.spec), s.freq)
Base.abs(s::Spectrum)  = Spectrum( abs.(s.spec), s.freq)

function add_box_inset(fig, left=65, right=200, bottom=65, top=150)
    inset_box = Axis(fig, bbox=BBox(left, right, bottom, top), 
    xgridvisible = false, ygridvisible = false, xscale = log10, yscale = log10,
    xticklabelsize=12, yticklabelsize=12, 
    xticks = ([], []),
    yticks = ([], [])
        )
    translate!(inset_box.scene, 0, 0, 10)  # bring content upfront
    return inset_box
end

function plot_spectra(trailing_character, jrange, krange; speckeys = nothing, 
                      plot_val = true, ax1 = nothing, ax2 = nothing, ax3 = nothing,
                      inset_ax = nothing)
    file     = jldopen(trailing_character * ".jld2")
    filekeys = keys(file)
    spectra  = Dict()

    if plot_val
        fig = (isnothing(ax1) && isnothing(ax2)) ? Figure() : nothing
        axE  = isnothing(ax1) ? Axis(fig[1, 1], xscale = log10, yscale = log10) : ax1
        axΩ  = isnothing(ax2) ? Axis(fig[1, 2], xscale = log10, yscale = log10) : ax2
        axW  = isnothing(ax3) ? Axis(fig[1, 3], xscale = log10, yscale = log10) : ax3
        iax  = isnothing(inset_ax) ? add_box_inset(fig) : inset_ax
    end

    for key in filekeys
        if isnothing(speckeys) || (key ∈ speckeys)
            spec = Spectrum[]
            for Var in [:U, :V, :Ω, :WB, :ST, :PV]
                @info "doing key $key, variable $Var"
                var = getproperty(file[key][1], Var)
                var = mean(real.(var[jrange, krange]))
                push!(spec, var)
            end
            @show key
            spectra[key] = NamedTuple{(:E, :Ω, :WB, :ST, :PV)}((spec[1] + spec[2], spec[3], spec[4], spec[5], spec[6]))
        end
    end

    if plot_val 
        for (key, spec) in spectra
            freq = spec[1].freq[2:end] * 1 / 100kilometers 
            spe3 = spec[3].spec[2:end]
            spe3[spe3 .<= 0] .= NaN
            if (plot_key(Val(Symbol(key))) && isnothing(speckeys)) || (key ∈ speckeys)
                lines!(iax, freq, spec[1].spec[2:end], color = clr(Val(Symbol(key))), linewidth = width(length(freq)), label = string(key), linestyle = sty(length(freq)))
                lines!(axE, freq, spec[1].spec[2:end], color = clr(Val(Symbol(key))), linewidth = width(length(freq)), label = string(key), linestyle = sty(length(freq)))
                lines!(axΩ, freq, spec[2].spec[2:end], color = clr(Val(Symbol(key))), linewidth = width(length(freq)), label = string(key), linestyle = sty(length(freq)))
                lines!(axW, freq, spe3,                color = clr(Val(Symbol(key))), linewidth = width(length(freq)), label = string(key), linestyle = sty(length(freq)))
            end
            if key == "weno9pV" && length(freq) > 100 && length(freq) < 200
                idx = searchsortedlast(freq, freq[end] / 10)
                x1  = freq[idx]
                s1  = spec[1].spec[idx]

                C = s1 / x1^(-3)
                x = freq[3:end]
                y = C .* x.^(-3)

                lines!(axE, x, y, linestyle = :dash, color = :gray, linewidth = 2)

                s2  = spec[2].spec[idx]

                C = s2 / x1^(-1)
                x = freq[3:end]
                y = C .* x.^(-1)

                lines!(axΩ, x, y, linestyle = :dash, color = :gray, linewidth = 2)
            end
            if key == "weno9pV" && length(freq) > 100 && length(freq) < 200
                smin1 = minimum(spec[1].spec)
                smin2 = minimum(spec[2].spec)
                smin3 = minimum(filter(!isnan, spe3))
                kw = [2π / (3.9 * 7e3), 2π / (3.9 * 7e3)]
                lines!(axE, kw, [smin1 * 0.8, 1e-3],  linestyle = :dashdot, color = :gray, linewidth = 2)
                lines!(axΩ, kw, [smin2 * 0.8, 1e-10], linestyle = :dashdot, color = :gray, linewidth = 2)
                lines!(axW, kw, [smin3 * 0.8, 1e-7],  linestyle = :dashdot, color = :gray, linewidth = 2)
            end
        end
        # axislegend(axE, position = :lb)
        # axislegend(axΩ, position = :lb)
        # axislegend(axW, position = :lb)
        xlims!(iax, (1e-5, 6e-5))
        ylims!(iax, (5e-3, 1))
        
        xlims!(axW, (2.5e-6, 5e-4))
        ylims!(axW, (1e-11,  1e-7))
        xlims!(axE, (2.5e-6, 5e-4))
        ylims!(axE, (1e-6,   1.8))
        xlims!(axΩ, (2.5e-6, 5e-4))
        ylims!(axΩ, (2e-13,  1e-9))
    end 

    close(file)
    return fig
end

res = ceil.(Ref(Int), (1000 * 0.9, 350 * 0.9))

function plot_more_spec(chars, sizes; keys = nothing)
    fig = Figure(resolution = res)

    if isnothing(keys)
        keys = Tuple(nothing for i in 1:length(sizes))
    end

    axE = Axis(fig[1, 1], xscale = log10, yscale = log10,
    # xlabel = L"\text{Wavenumber m}^{-1}",
    # title = L"\text{Energy Spectra}",
    ylabel = "",
    xgridvisible = false, ygridvisible = false,
    # xticks = ([1e-5, 1e-4, 1e-3, 1e-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    # yticks = ([1e-6, 1e-4, 1e-2, 1],    [L"10^{-6}", L"10^{-4}", L"10^{-3}", L"10^{0}"]),
    xtickalign = 1,
    ytickalign = 1)

    axΩ = Axis(fig[1, 2], xscale = log10, yscale = log10,
    # xlabel = L"\text{Wavenumber m}^{-1}",
    # title = L"\text{Enstrophy Spectra}",
    ylabel = "",
    xgridvisible = false, ygridvisible = false,
    # xticks = ([1e-5, 1e-4, 1e-3, 1e-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    # yticks = ([1e-13, 1e-11, 1e-9], [L"10^{-13}", L"10^{-11}", L"10^{-9}"]),
    xtickalign = 1,
    ytickalign = 1)

    axW = Axis(fig[1, 3], xscale = log10, yscale = log10,
    # xlabel = L"\text{Wavenumber m}^{-1}",
    # title = L"wb \text{ cospectra}",
    ylabel = "",
    xgridvisible = false, ygridvisible = false,
    # xticks = ([1e-5, 1e-4, 1e-3, 1e-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    # yticks = ([1e-11, 1e-10, 1e-9, 1e-8, 1e-7], [L"10^{-11}", L"10^{-10}", L"10^{-9}", L"10^{-8}", L"10^{-8}"]),
    xtickalign = 1,
    ytickalign = 1)

    
    inset_ax = add_box_inset(fig)

    for (char, size, key) in zip(chars, sizes, keys)
        plot_spectra(char, size...; speckeys = key, ax1 = axE, ax2 = axΩ, ax3 = axW, inset_ax)
    end
    return fig
end

const color1 = :deepskyblue
const color4 = :deepskyblue2
const color2 = :orange1
const color3 = :firebrick2

clr(::Val{s}) where s = string(s)[1:4] == "smag" ? Symbol("green") :
                        string(s)[1:5] == "leith" ? RGBf(0.8, 0.8, 0.8) :
                        string(s)[1:2] == "qg" ? Symbol("green3") :
                        string(s) == "bilap" ? color2 :
                        string(s) == "weno5pV" ? Symbol("deepskyblue") :
                        string(s) == "weno9pV" ? Symbol("purple") :                        
                        string(s) == "weno9pAllD" ? Symbol("blue") :
                        string(s)[1:7] == "weno5Fl" ? Symbol("gray10") :
                        string(s)[1:7] == "weno9Fl" ? Symbol("orange1") :
                        string(s) == "highres" ? RGBf(0.8, 0.8, 0.8) : color1 

plot_key(::Val{s}) where s = string(s) == "weno5vd" || 
                             string(s) == "lapleith" || 
                             string(s) == "bilap" || 
                             string(s) == "highres" ? true : false
                        
sty(s) = s < 50 ? :solid : s < 100 ? :solid : s < 250 ? :solid : :solid
width(s) = s < 50 ? 0.75 : s < 100 ? 1.8 : s < 250 ? 3.2 : 7