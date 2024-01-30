using BaroclinicAdjustment
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Advection: DefaultStencil, FunctionStencil, AUGX, AUGY, AUGZ, HOADV, MPData
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using BaroclinicAdjustment: geometric_νhb

stop_time = 1000days
buoyancy_forcing_timescale = 50days

advs, hors, names = BaroclinicAdjustment.testcases(Float64)

iter = get(ENV, "CASE", "5")

if iter == "Upwind"
   advs = [UpwindBiased(order = 3)]
   hors = [nothing]
   names = ["upwind3"]
elseif iter == "MPData"
   advs = [nothing]
   hors = [nothing]
   names = ["mpdata_full"]
else
   iter = parse(Int, iter)
   iter = [iter] #[14] #, 19, 2]
   advs = advs[iter]
   hors = hors[iter]
   names = names[iter]
end

hors = [HorizontalScalarBiharmonicDiffusivity(Float64; ν = geometric_νhb, discrete_form = true, parameters = 10days)]

for (res, trl) in zip((1/32, ), ("_thirtytwo", )), (momentum_advection, horizontal_closure, name) in zip(advs, hors, names)
    @show name    
    simulation = BaroclinicAdjustment.baroclinic_adjustment_latlong(res, name * trl; 
                                                                    arch = GPU(), momentum_advection,
                                                                    stop_time, horizontal_closure,
                                                                    tracer_advection = WENO(order = 7),
                                                                    buoyancy_forcing_timescale)
    run!(simulation)
end

