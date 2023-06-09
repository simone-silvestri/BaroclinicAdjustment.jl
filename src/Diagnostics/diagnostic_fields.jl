using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, ∇_dot_qᶜ, AbstractScalarBiharmonicDiffusivity, ExplicitTimeDiscretization

import Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

VerticalVorticityField(fields::Dict, i) = VerticalVorticityField((; u = fields[:u][i], v = fields[:v][i]))
KineticEnergyField(fields::Dict, i)     =     KineticEnergyField((; u = fields[:u][i], v = fields[:v][i]))

VerticalVorticityOperation(fields::Dict, i)   =   VerticalVorticityOperation((; u = fields[:u][i], v = fields[:v][i]))
PotentialVorticityOperation(fields::Dict, i)  =  PotentialVorticityOperation((; u = fields[:u][i], v = fields[:v][i], b = fields[:b][i]))
VerticalDissipationOperation(fields::Dict, i) = VerticalDissipationOperation((; u = fields[:u][i], v = fields[:v][i], b = fields[:b][i]))
KineticEnergyOperation(fields::Dict, i)       =       KineticEnergyOperation((; u = fields[:u][i], v = fields[:v][i]))
StratificationOperation(fields::Dict, i)      =      StratificationOperation(fields[:b][i])

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

@inline _density_operation(i, j, k, grid, b, ρ₀, g) = ρ₀ * (1 - b[i, j, k] / g)

DensityOperation(b; ρ₀ = 1000.0, g = 9.80655) = 
    KernelFunctionOperation{Center, Center, Center}(_density_operation, b.grid, b, ρ₀, g)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(DensityOperation(b; ρₒ, g)))

function HeightField(grid, loc = (Center, Center, Center))  

    zf = Field(loc, grid)
    Lz = grid.Lz

    for k in 1:size(zf, 3)
        interior(zf, :, :, k) .= Lz + znode(k, grid, loc[3]())
    end

    return zf
end

function KineticEnergyField(velocities::NamedTuple)
    u = velocities.u
    v = velocities.v

    E_op = @at (Center, Center, Center) 0.5 * (u^2 + v^2)

    return compute!(Field(E_op))
end

function VerticalVorticityOperation(velocities::NamedTuple)

    grid = velocities.u.grid
    computed_dependencies = (velocities.u, velocities.v)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, computed_dependencies...)

    return ζ_op
end

function StratificationOperation(b)
    grid = b.grid

    N2_op = KernelFunctionOperation{Center, Center, Face}(N²ᶜᶜᶠ, grid, b)

    return N2_op
end

@inline N²ᶠᶠᶠ(i, j, k, grid, b) = max(1e-10, ℑxyᶠᶠᵃ(i, j, k, grid, ∂zᶜᶜᶠ, b))
@inline N²ᶜᶜᶠ(i, j, k, grid, b) = ∂zᶜᶜᶠ(i, j, k, grid, b)

@inline b_term(i, j, k, grid, b) = fᶠᶠᵃ(i, j, k, grid, HydrostaticSphericalCoriolis()) / N²ᶠᶠᶠ(i, j, k, grid, b) * ℑxyᶠᶠᵃ(i, j, k, grid, b)
@inline pvᶠᶠᶜ(i, j, k, grid, u, v, b) = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) + ∂zᶠᶠᶜ(i, j, k, grid, b_term, b) 

function PotentialVorticityOperation(fields::NamedTuple)

    grid = fields.u.grid
    computed_dependencies = (fields.u, fields.v, fields.b)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(pvᶠᶠᶜ, grid, computed_dependencies...)

    return ζ_op
end

function KineticEnergyOperation(velocities::NamedTuple)
    u = velocities.u
    v = velocities.v

    E_op = @at (Center, Center, Center) 0.5 * (u^2 + v^2)

    return E_op
end

VerticalVorticity(f::Dict, i) = compute!(Field(VerticalVorticityOperation(f, i)))
KineticEnergy(f::Dict, i)     = compute!(Field(KineticEnergyOperation(f, i)))
Stratification(f::Dict, i)    = compute!(Field(StratificationOperation(f, i)))

@inline _deformation_radius(i, j, k, grid, b) = sqrt(max(0, ∂zᶜᶜᶠ(i, j, k, grid, b))) / π /
                                                abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, HydrostaticSphericalCoriolis()))

function DeformationRadius(f::Dict, i)
    
    Rop = KernelFunctionOperation{Center, Center, Face}(_deformation_radius, f[:b].grid, f[:b][i])
    R   = compute!(Field(Integral(Rop, dims = 3))) 

    return R
end

function HorizontalFriction(model; ClosureType = AbstractScalarBiharmonicDiffusivity)

    grid          = model.grid
    clock         = model.clock
    closure       = filter(x -> x isa ClosureType, model.closure)
    diffusivities = model.diffusivity_fields
    buoyancy      = model.buoyancy
    velocities    = model.velocities
    free_surface  = model.free_surface
    tracers       = model.tracers
    auxiliary_fields = model.auxiliary_fields

    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    computed_dependencies = (closure, diffusivities, clock, model_fields, buoyancy)

    ∂ⱼ_τ₁ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₁ⱼ, grid, computed_dependencies...)
    ∂ⱼ_τ₂ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₂ⱼ, grid, computed_dependencies...)
    ∂ⱼ_τ₃ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₃ⱼ, grid, computed_dependencies...)

    τ₁ = compute!(Field(∂ⱼ_τ₁ⱼ_op))
    τ₂ = compute!(Field(∂ⱼ_τ₂ⱼ_op))
    τ₃ = compute!(Field(∂ⱼ_τ₃ⱼ_op))

    return (; τ₁, τ₂, τ₃)
end

function HorizontalFriction(fields::NamedTuple, closure)

    grid          = fields.u.grid
    clock         = Clock{eltype(grid)}(0, 0, 1)
    diffusivities = nothing
    buoyancy      = BuoyancyTracer()
    velocities    = (; u = fields.u, v = fields.v, w = ZFaceField(grid))
    free_surface  = nothing
    tracers       = (; b = fields.b)
    auxiliary_fields = NamedTuple()

    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    computed_dependencies = (closure, diffusivities, clock, model_fields, buoyancy)

    ∂ⱼ_τ₁ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₁ⱼ, grid, computed_dependencies...)
    ∂ⱼ_τ₂ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₂ⱼ, grid, computed_dependencies...)
    ∂ⱼ_τ₃ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₃ⱼ, grid, computed_dependencies...)

    τ₁ = compute!(Field(∂ⱼ_τ₁ⱼ_op))
    τ₂ = compute!(Field(∂ⱼ_τ₂ⱼ_op))
    τ₃ = compute!(Field(∂ⱼ_τ₃ⱼ_op))

    return (; τ₁, τ₂, τ₃)
end

function VerticalDissipationOperation(fields::NamedTuple)

    grid          = fields.b.grid
    clock         = Clock{eltype(grid)}(0, 0, 1)
    diffusivities = nothing
    closure       = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), κ = 1e-5)
    buoyancy      = BuoyancyTracer()
    velocities    = (; u = fields.u, v = fields.v, w = nothing)
    tracers       = (; b = fields.b)

    model_fields = merge(velocities, tracers)
    computed_dependencies = (closure, diffusivities, Val(1), clock, model_fields, buoyancy)

    ∇_dot_qᶜ_op = KernelFunctionOperation{Face, Center, Center}(∇_dot_qᶜ, grid, computed_dependencies...)

    return ∇_dot_qᶜ_op * fields.b
end

