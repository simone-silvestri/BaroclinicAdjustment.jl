using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalFormulation
using Oceananigans.Operators
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz

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

function leith_viscosity(formulation, FT::DataType = Float64; Cₗ = FT(1.0), Area = Δ²ᵃᵃᵃ)

    @show C = (Cₗ / π)^3 / 8

    visc = ScalarBiharmonicDiffusivity(formulation, FT; 
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

function leith_laplacian_viscosity(formulation = HorizontalFormulation(), FT::DataType = Float64; Cₗ = FT(1.0), Area = Δ²ᵃᵃᵃ)

    @show C = (Cₗ / π)^3 

    visc = ScalarDiffusivity(formulation, FT; 
                             ν=νhb_leith_laplacian_final, discrete_form=true,  
                             parameters = (; C, Area))

    @show typeof(visc.ν)

    return visc
end

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ

