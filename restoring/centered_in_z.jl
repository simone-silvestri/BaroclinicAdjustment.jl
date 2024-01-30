
struct ThreeDimensionalAdvection{N, FT, A, B, C} <: AbstractAdvectionScheme{N, FT}
    x :: A
    y :: B
    z :: C

    ThreeDimensionalAdvection{N, FT}(x::A, y::B, z::C) where {N, FT, A, B, C} = new{N, FT, A, B, C}(x, y, z)
end

function ThreeDimensionalAdvection(; x, y, z)
    Nx = required_halo_size(x) 
    Ny = required_halo_size(y) 
    Nz = required_halo_size(z) 

    FT = eltype(x)

    return ThreeDimensionalAdvection{max(Nx, Ny, Nz), FT}(x, y, z)
end

import Oceananigans.Advection: div_Uc

@inline function div_Uc(i, j, k, grid, advection::ThreeDimensionalAdvection, U, c)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection.x, U.u, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection.y, U.v, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, _advective_tracer_flux_z, advection.z, U.w, c))
end

momentum_advection = VectorInvariant(vorticity_scheme = WENO(order = 9),
                                     vertical_scheme = WENO(order = 5))
tracer_advection = ThreeDimensionalAdvection(x = WENO(order = 7, buffer_scheme = Centered()), 
                                             y = WENO(order = 7, buffer_scheme = Centered()), 
                                             z = Centered())

horizontal_closure = nothing

res = 1/16
trl = "_sixteen_center_buffer"

name = "weno9pV"

simulation = BaroclinicAdjustment.baroclinic_adjustment_latlong(res, name * trl; 
                                                                    arch = GPU(), momentum_advection,
                                                                    stop_time, horizontal_closure,
                                                                    tracer_advection,
                                                                    buoyancy_forcing_timescale)
run!(simulation)
                    