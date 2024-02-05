using CairoMakie

include("buoyancy_plots.jl")
CairoMakie.activate!()
CairoMakie.save("buoyancy_fig.eps", figb)

include("buoyancy_contour.jl")
CairoMakie.activate!()
CairoMakie.save("buoyancy_contour.eps", fig)

include("energy_plots.jl")
CairoMakie.activate!()
CairoMakie.save("energies_3d.eps", fig)
