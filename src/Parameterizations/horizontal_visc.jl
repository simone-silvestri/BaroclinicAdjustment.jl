using Oceananigans.TurbulenceClosures
using Oceananigans.Operators: Δx, Δy

GeometricBilaplacian(FT; λ = 5days) = 
    HorizontalScalarBiharmonicDiffusivity(FT; 
                                          ν=geometric_νhb, discrete_form=true,  
                                          parameters = λ)

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ

