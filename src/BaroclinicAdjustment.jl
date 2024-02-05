module BaroclinicAdjustment

export TestCase
export baroclinic_adjustment_simulation

using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing, architecture
using Oceananigans.Operators
using Oceananigans.Utils: getnamewrapper, launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Advection: VelocityStencil, DefaultStencil, EnergyConserving

using Oceananigans.Advection: FunctionStencil, divergence_smoothness
using Oceananigans.Advection: CrossAndSelfUpwinding, OnlySelfUpwinding, VelocityUpwinding

using KernelAbstractions: @kernel, @index
using JLD2
using Random

struct TestCase{A, H, N}
    a :: A # momentum advection scheme
    h :: H # lateral friction
    n :: N # name
end

function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.75)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δ  = 1 / sqrt(1 / Δx^2 + 1 / Δy^2)

    return  Base.Int(ceil(2 * Δt / (CFL / wave_speed * Δ)))
end

include("Parameterizations/Parameterizations.jl")
include("outputs.jl")
include("baroclinic_adjustment.jl")
include("Diagnostics/Diagnostics.jl")
include("Postprocess/Postprocess.jl")

using .Parameterizations
using .Diagnostics
using .Postprocess

end
