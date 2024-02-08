using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, ∇_dot_qᶜ, AbstractScalarBiharmonicDiffusivity, ExplicitTimeDiscretization

 VerticalVorticityOperation(fields::Dict, i) = VerticalVorticityOperation((; u = fields[:u][i], v = fields[:v][i]))
PotentialVorticityOperation(fields::Dict, i) = PotentialVorticityOperation((; u = fields[:u][i], v = fields[:v][i], b = fields[:b][i]))
     KineticEnergyOperation(fields::Dict, i) = KineticEnergyOperation((; u = fields[:u][i], v = fields[:v][i]))
    StratificationOperation(fields::Dict, i) = StratificationOperation(fields[:b][i])
          RiNumberOperation(fields::Dict, i) = RiNumberOperation((; b = fields[:b][i], u = fields[:u][i], v = fields[:v][i]))
 
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

function VerticalVorticityOperation(velocities::NamedTuple)

    grid = velocities.u.grid
    computed_dependencies = (velocities.u, velocities.v)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, computed_dependencies...)

    return ζ_op
end

function StratificationOperation(b)
    grid = b.grid
    loc  = location(b)

    N2_op = KernelFunctionOperation{loc[1], loc[2], Face}(N²ᶜᶜᶠ, grid, b)

    return N2_op
end

@inline N²ᶜᶜᶠ(i, j, k, grid, b)    = ∂zᶜᶜᶠ(i, j, k, grid, b)
@inline fᶜᶜᵃ(i, j, k, grid)        = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, HydrostaticSphericalCoriolis())
@inline ζ₃ᶜᶜᶜ(i, j, k, grid, u, v) = ℑxyᶜᶜᵃ(i, j, k, grid, ζ₃ᶠᶠᶜ, u, v)

@inline pvᶜᶜᶜ(i, j, k, grid, u, v, b) = (ζ₃ᶜᶜᶜ(i, j, k, grid, u, v) + fᶜᶜᵃ(i, j, k, grid)) * ℑzᵃᵃᶜ(i, j, k, grid, N²ᶜᶜᶠ, b) +
                                         ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u) * ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b) + 
                                         ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
                                        

function PotentialVorticityOperation(fields::NamedTuple)

    grid = fields.u.grid
    computed_dependencies = (fields.u, fields.v, fields.b)

    ζ_op = KernelFunctionOperation{Center, Center, Center}(pvᶜᶜᶜ, grid, computed_dependencies...)

    return ζ_op
end

function KineticEnergyOperation(velocities::NamedTuple)
    u = velocities.u
    v = velocities.v

    E_op = @at (Center, Center, Center) 0.5 * (u^2 + v^2)

    return E_op
end

using Oceananigans.TurbulenceClosures: Riᶜᶜᶠ

function RiNumberOperation(fields::NamedTuple)

    grid = fields.u.grid
    computed_dependencies = ((; u = fields.u,  v = fields.v), BuoyancyTracer(), (; b = fields.b))

    Ri_op = KernelFunctionOperation{Center, Center, Center}(Riᶜᶜᶠ, grid, computed_dependencies...)

    return Ri_op
end

PotentialVorticity(f::Dict, i) = compute!(Field(PotentialVorticityOperation(f, i)))
VerticalVorticity(f::Dict, i)  = compute!(Field(VerticalVorticityOperation(f, i)))
KineticEnergy(f::Dict, i)      = compute!(Field(KineticEnergyOperation(f, i)))
Stratification(f::Dict, i)     = compute!(Field(StratificationOperation(f, i)))
RiNumber(f::Dict, i)           = compute!(Field(RiNumberOperation(f, i)))

@inline _deformation_radius(i, j, k, grid, b) = sqrt(max(0, ∂zᶜᶜᶠ(i, j, k, grid, b))) / π /
                                                abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, HydrostaticSphericalCoriolis()))

function DeformationRadius(f::Dict, i)
    
    Rop = KernelFunctionOperation{Center, Center, Face}(_deformation_radius, f[:b].grid, f[:b][i])
    R   = compute!(Field(Integral(Rop, dims = 3))) 

    return R
end