using BaroclinicAdjustment
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Advection: EnergyConserving, DefaultStencil
using BaroclinicAdjustment.Parameterizations
using BaroclinicAdjustment: baroclinic_adjustment_latlong

using CUDA
CUDA.device!(1)

stop_time = 1000days
buoyancy_forcing_timescale = 50days

energy_conserving_advection = VectorInvariant(vorticity_scheme = EnergyConserving(), 
                                               vertical_scheme = EnergyConserving())

best_weno = VectorInvariant(vorticity_scheme = WENO(; order = 9),
                             vertical_scheme = WENO(),
                           divergence_scheme = WENO())

worst_upwinding = OnlySelfUpwinding(; cross_scheme = WENO(FT),
                                      δU_stencil  = DefaultStencil(),
                                      δV_stencil  = DefaultStencil(),
                                      δu²_stencil = DefaultStencil(),
                                      δv²_stencil = DefaultStencil())                               

worst_weno = VectorInvariant(vorticity_scheme = WENO(; order = 9),
                            vorticity_stencil = DefaultStencil(),
                              vertical_scheme = WENO(),
                            divergence_scheme = WENO(),
                                    upwinding = worst_upwinding)

# Define test cases:

qgleith_test    = TestCase(energy_conserving_advection, QGLeith(), "qgleith")
omp25_test      = TestCase(energy_conserving_advection, OMp25Closure(), "omp25")
best_weno_test  = TestCase(best_weno,  nothing, "weno9pV")
worst_weno_test = TestCase(worst_weno, nothing, "weno9pAllD")
upwind_test     = TestCase(UpwindBiased(order = 3), nothing, "upwind")

tests = (qgleith_test, omp25_test, best_weno_test, worst_weno_test, upwind_test)

resolutions = (1/8, 1/16, 1/32)
resnames    = ("_eight", "_sixteen", "_thirtytwo")

for (resolution, trailing_character) in zip(resolutions, resnames)
   for test in tests    
      simulation = baroclinic_adjustment_latlong(test, resolution, trailing_character; arch = GPU())
      run!(simulation)
   end
end


