using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalFormulation
using Oceananigans.Operators
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz

@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz) =  2 * (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

@inline function νhb_smagorinski_final(i, j, k, grid, lx, ly, lz, clock, fields, p)

   location = (lx, ly, lz)
   from_Dₛ = (Center(), Center(), Center()) 
   from_Dₜ = (Face(),   Face(),   Center()) 
	
   δ₁ = ℑxyz(i, j, k, grid, from_Dₛ, location, Dₛ, fields.u, fields.v)    
   δ₂ = ℑxyz(i, j, k, grid, from_Dₜ, location, Dₛ, fields.u, fields.v)    

   dynamic_visc = p.C * sqrt(δ₁^2 + δ₂^2)

   A = p.Area(i, j, k, grid, lx, ly, lz)

   return A^2 * dynamic_visc
end

function smagorinski_viscosity(formulation; Cₛₘ = 0.45, Area = Δ²ᵃᵃᵃ)

    @show C = (Cₛₘ / π)^2 / 8

    return ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_smagorinski_final, discrete_form=true,  
				                       parameters = (; C, Area))
end

@inline function νhb_leith_final(i, j, k, grid, lx, ly, lz, clock, fields, p)
    
    location = (lx, ly, lz)
    from_∂xζ = (Center(), Face(), Center()) 
    from_∂yζ = (Face(), Center(), Center()) 
    from_∂xδ = (Face(), Center(), Center()) 
    from_∂yδ = (Center(), Face(), Center()) 
	
    ∂xζ = ℑxyz(i, j, k, grid, from_∂xζ, location, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxyz(i, j, k, grid, from_∂yζ, location, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂xδ = ℑxyz(i, j, k, grid, from_∂xδ, location, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑxyz(i, j, k, grid, from_∂yδ, location, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
   
    dynamic_visc = sqrt((∂xζ^2 + ∂yζ^2) + (∂xδ^2 + ∂yδ^2) ) * p.C
 
    A = p.Area(i, j, k, grid, lx, ly, lz)
    
    return dynamic_visc * A^(5/2)
end

function leith_viscosity(formulation; Cₗ = 1.0, Area = Δ²ᵃᵃᵃ)

    @show C = (Cₗ / π)^3 / 8

    visc = ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_leith_final, discrete_form=true,  
                                       parameters = (; C, Area))

    @show typeof(visc.ν)

    return visc
end

@inline function νhb_leith_laplacian_final(i, j, k, grid, lx, ly, lz, clock, fields, p)
    
    location = (lx, ly, lz)
    from_∂xζ = (Center(), Face(), Center()) 
    from_∂yζ = (Face(), Center(), Center()) 
    from_∂xδ = (Face(), Center(), Center()) 
    from_∂yδ = (Center(), Face(), Center()) 
	
    ∂xζ = ℑxyz(i, j, k, grid, from_∂xζ, location, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxyz(i, j, k, grid, from_∂yζ, location, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂xδ = ℑxyz(i, j, k, grid, from_∂xδ, location, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑxyz(i, j, k, grid, from_∂yδ, location, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
   
    dynamic_visc = sqrt((∂xζ^2 + ∂yζ^2) + (∂xδ^2 + ∂yδ^2)) * p.C
 
    A = p.Area(i, j, k, grid, lx, ly, lz)

    return dynamic_visc * A^(3/2)
end

function leith_laplacian_viscosity(formulation = HorizontalFormulation(); Cₗ = 1.0, Area = Δ²ᵃᵃᵃ)

    @show C = (Cₗ / π)^3 

    visc = ScalarDiffusivity(formulation; 
                             ν=νhb_leith_laplacian_final, discrete_form=true,  
                             parameters = (; C, Area))

    @show typeof(visc.ν)

    return visc
end

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ

