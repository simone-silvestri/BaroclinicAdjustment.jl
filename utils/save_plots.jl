using CairoMakie

include("../utils/plotting.jl")
include("../utils/average_spectra.jl")

# fig1, fig2, fig3 = all_figures()

# CairoMakie.save("energy_plots.eps", fig1)
# CairoMakie.save("energy_stratif.eps", fig2)
# CairoMakie.save("enstrophy_plots.eps", fig3)

# fig1 = all_contours_white()

# CairoMakie.save("white_contours.eps", fig1)

# fig1, fig2, fig3, fig4 = all_contours()

# CairoMakie.save("energy_contours.png", fig1)
# CairoMakie.save("vorticity_contours.png", fig2)
# CairoMakie.save("enstrophy_contours.png", fig3)
# CairoMakie.save("buoyancy_contours.png", fig4)

key1 = ["qgleith", "weno9pV"]

ffig = plot_more_spec(["spectra_fifty", "spectra_quarter", "spectra_eight", "spectra_sixteen"],
                      [(10:990, 1:50), (10:70, 1:50), (10:150, 1:50), (10:310, 1:50)]; 
                       keys = (["leith"], key1, key1, key1))
                    
CairoMakie.save("spectra_plots.eps", ffig)
