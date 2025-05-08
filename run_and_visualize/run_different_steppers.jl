using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: EnergyConserving, DefaultStencil, OnlySelfUpwinding
using BaroclinicAdjustment
using BaroclinicAdjustment.PostProcess
using BaroclinicAdjustment.Parameterizations

####
#### Advection schemes
####

# Classical energy conserving vector-invariant advection
energy_conserving = VectorInvariant(vorticity_scheme = EnergyConserving(), 
                                     vertical_scheme = EnergyConserving())

# W9V advection scheme
best_weno = WENOVectorInvariant()
omp25     = OMp25Closure()

####
#### Construct test cases
####

# OMP25 test case AB2
omp25AB2_test = TestCase(energy_conserving, omp25, "omp25AB2", :QuasiAdamsBashforth2)

# OMP25 test case RK3
omp25RK3_test = TestCase(energy_conserving, omp25, "omp25RK3", :SplitRungeKutta3)

# W9V test case AB2
w9vAB2_test  = TestCase(best_weno, nothing, "wenoAB2", :QuasiAdamsBashforth2)

# W9V test case RK3
w9vRK3_test  = TestCase(best_weno, nothing, "wenoRK3", :SplitRungeKutta3)


# All test cases
all_tests = (omp25AB2_test, 
             omp25RK3_test, 
             w9vAB2_test,
             w9vRK3_test)

####
#### Let's run an postprocess!
####

for test in all_tests    
   # Define the simulation
   # simulation = BaroclinicAdjustment.baroclinic_adjustment_simulation(test, 1/16, "_sixteen"; 
   #                                                                    arch = GPU(),
   #                                                                    buoyancy_forcing_timescale = nothing,
   #                                                                    stop_time = 500days)
   
   # # Let's run
   # run!(simulation)

   # Postprocessing the outputs
   Postprocess.calculate_diagnostics(test, "_sixteen")
end


