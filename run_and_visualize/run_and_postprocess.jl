using Oceananigans
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
best_weno = VectorInvariant(vorticity_scheme = WENO(; order = 9),
                             vertical_scheme = WENO(),
                           divergence_scheme = WENO())

# Upwinding using default stencils
default_stencils = OnlySelfUpwinding(; cross_scheme = WENO(),
                                       δU_stencil  = DefaultStencil(),
                                       δV_stencil  = DefaultStencil(),
                                       δu²_stencil = DefaultStencil(),
                                       δv²_stencil = DefaultStencil())                               

# W9D advection scheme
worst_weno = VectorInvariant(vorticity_scheme = WENO(; order = 9),
                            vorticity_stencil = DefaultStencil(),
                              vertical_scheme = WENO(),
                            divergence_scheme = WENO(),
                                    upwinding = default_stencils)

# Flux form upwind third order
upwind = UpwindBiased(order = 3)

#### 
#### Explicit closures
#### 

qgleith = QGLeith(; C = 2)
omp25   = OMp25Closure()
bileith = BiharmonicLeith(; C = 2)
ebs     = EnergyBackScatter()

####
#### Construct test cases
####

# QG2 test case
qgleith_test = TestCase(energy_conserving, qgleith, "qgleith")

# OMP25 test case
omp25_test = TestCase(energy_conserving, omp25, "omp25")

# BL2 test case
bileith_test = TestCase(energy_conserving, bileith, "bileith")

# EBS test case
ebs_test = TestCase(energy_conserving, ebs, "ebs")

# W9V test case
best_weno_test  = TestCase(best_weno, nothing, "weno9pV")

# W9D test case 
worst_weno_test = TestCase(worst_weno, nothing, "weno9pAllD")

# Upwind test case
upwind_test = TestCase(upwind, nothing, "upwind")

# All test cases
all_tests = (qgleith_test, 
             omp25_test, 
             bileith_test,
             ebs_test,
             best_weno_test, 
             worst_weno_test, 
             upwind_test)

####
#### Let's run an postprocess!
####

# All resolutions
resolutions = (1/8, 1/16, 1/32)
resnames    = ("_eight", "_sixteen", "_thirtytwo")

for (resolution, trailing_character) in zip(resolutions, resnames)
   for test in all_tests    
      # Define the simulation
      simulation = BaroclinicAdjustment.baroclinic_adjustment_latlong(test, resolution, trailing_character; arch = GPU())
      
      # Let's run
      run!(simulation)

      # Postprocessing the outputs
      Postprocess.calculate_diagnostics(test, trailing_character)
   end
end


