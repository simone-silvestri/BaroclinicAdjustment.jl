module BaroclinicAdjustment

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing, architecture
using Oceananigans.Operators
using Oceananigans.Utils: getnamewrapper, launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Advection: VelocityStencil, DefaultStencil, EnergyConservingScheme

using Oceananigans.Advection: FunctionStencil, divergence_smoothness
using Oceananigans.Advection: CrossAndSelfUpwinding, OnlySelfUpwinding, VelocityUpwinding

using Oceananigans.Advection: VectorInvariantCrossVerticalUpwinding, VectorInvariantSelfVerticalUpwinding, VectorInvariantVelocityVerticalUpwinding

using KernelAbstractions: @kernel, @index
using JLD2
using Random

include("horizontal_visc.jl")
include("qg_leith_viscosity.jl")
include("outputs.jl")

function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.75)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δ  = 1 / sqrt(1 / Δx^2 + 1 / Δy^2)

    return  Base.Int(ceil(2 * Δt / (CFL / wave_speed * Δ)))
end

viscosity_name(clo::ScalarBiharmonicDiffusivity) = typeof(clo.ν)
viscosity_name(clo::ScalarDiffusivity)           = typeof(clo.ν)
viscosity_name(clo) = typeof(clo)
advection_name(adv) = getnamewrapper(adv.vorticity_scheme)

add_trailing_characters(name, trailing_character = "_weaker") = name * trailing_character

getupwindingscheme(::VectorInvariantCrossVerticalUpwinding)    = "f"
getupwindingscheme(::VectorInvariantSelfVerticalUpwinding)     = "p"
getupwindingscheme(::VectorInvariantVelocityVerticalUpwinding) = "s"

getname(s::VectorInvariant{N}) where N = "weno" * string(N * 2 - 1) * getupwindingscheme(s) * "$(s.vorticity_stencil isa VelocityStencil ? "V" : "D")"

function testcases(FT)

    advection_schemes   = []
    horizontal_closures = []
    names               = []

    # First five are the "Explicit" LES closures
    for i in 1:5
        push!(advection_schemes, VectorInvariant(vorticity_scheme = EnergyConservingScheme(), vertical_scheme = EnergyConservingScheme()))
    end

    hi2 = HorizontalScalarBiharmonicDiffusivity(FT; ν = geometric_νhb, discrete_form = true, parameters = 5days)
    hi3 = leith_viscosity(HorizontalFormulation(), FT)
    hi4 = leith_laplacian_viscosity(HorizontalFormulation(), FT)
    hi5 = smagorinski_viscosity(HorizontalFormulation(), FT)
    hi6 = QGLeith(FT)

    push!(horizontal_closures, hi2, hi3, hi4, hi5, hi6)
    push!(names, "bilap", "leith", "lapleith", "smag", "qgleith")

    # Next twelve are the "Implicit" LES closures
    for order in [5, 9]
        for Upwind in (CrossAndSelfUpwinding, OnlySelfUpwinding, VelocityUpwinding)
            for vorticity_stencil in (VelocityStencil(), DefaultStencil())
                push!(advection_schemes, VectorInvariant(; vorticity_scheme = WENO(FT; order), 
                                                           vorticity_stencil,
                                                           vertical_scheme = WENO(FT), 
                                                           upwinding= Upwind(cross_scheme = WENO(FT))))
                push!(horizontal_closures, nothing)
                push!(names, getname(advection_schemes[end]))
            end
        end
    end

    # Incorrect stencil usage
    for order in [5, 9]
        push!(advection_schemes, VectorInvariant(; vorticity_scheme = WENO(FT; order),  
                                                   vorticity_stencil = DefaultStencil(),
                                                   upwinding= OnlySelfUpwinding(;
                                                                                cross_scheme = WENO(FT),
                                                                                δU_stencil  = DefaultStencil(),
                                                                                δV_stencil  = DefaultStencil(),
                                                                                δu²_stencil = DefaultStencil(),
                                                                                δv²_stencil = DefaultStencil()),
                                                   vertical_scheme = WENO(FT)))
                
        push!(horizontal_closures, nothing)
    end

    push!(names, "weno5pAllD", "weno9pAllD")
    
    # Flux form WENO schemes
    push!(advection_schemes, WENO(FT))
    push!(advection_schemes, WENO(FT; order = 9))
    push!(horizontal_closures, nothing, nothing)

    push!(names, "weno5Fl", "weno9Fl")

    push!(advection_schemes, VectorInvariant(vorticity_scheme = WENO(FT),            vertical_scheme = WENO(FT), multi_dimensional_stencil = true))
    push!(advection_schemes, VectorInvariant(vorticity_scheme = WENO(FT; order = 9), vertical_scheme = WENO(FT), multi_dimensional_stencil = true))
    push!(horizontal_closures, nothing, nothing)
    push!(names, "weno5MD", "weno9MD")

    return advection_schemes, horizontal_closures, names
end

include("baroclinic_adjustment_latlon.jl")
include("restoring_baroclinic_adjustment_latlon.jl")
include("run_simulations.jl")
include("Diagnostics/Diagnostics.jl")

using .Diagnostics

end
