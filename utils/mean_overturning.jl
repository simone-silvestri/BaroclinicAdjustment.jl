using Oceananigans
using BaroclinicAdjustment
using KernelAbstractions
using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

wn08 = jldopen("weno9pV_eight_new_postprocess.jld2")
qg08 = jldopen("qgleith2_eight_postprocess.jld2")
up08 = jldopen("upwind_eight_postprocess.jld2")
eb08 = jldopen("ebs_eight_postprocess.jld2")
om08 = jldopen("omp25_eight_postprocess.jld2")
wd08 = jldopen("weno9pAllD_eight_new_postprocess.jld2")
qg16 = jldopen("qgleith2_sixteen_postprocess.jld2")
up16 = jldopen("upwind_sixteen_postprocess.jld2")
om16 = jldopen("omp25_sixteen_postprocess.jld2")
wd16 = jldopen("weno9pAllD_sixteen_new_postprocess.jld2")
wn16 = jldopen("weno9pV_sixteen_new_postprocess.jld2")
qg32 = jldopen("qgleith_thirtytwo_new_postprocess.jld2")
wn32 = jldopen("weno9pV_thirtytwo_new_postprocess.jld2")

@kernel function _calculate_ψ(ψ, V, Δz, ::Val{Nz}) where Nz
    j = @index(Global, Linear)
    
    ψ[1, j, 1] = 0

    @unroll for k in 2:Nz+1
        ψ[1, j, k] = ψ[1, j, k-1] + Δz * V[1, j, k-1]
    end
end

function overturning(mean_file)
    V    = mean_file["mean"].V
    grid = V.grid
    Nz = grid.Nz
    Ny = grid.Ny
    Δz = 20

    ψ = Field((Nothing, Face, Face), grid)

    _calculate_ψ(KernelAbstractions.CPU(), 16, Ny+1)(ψ, V, Δz, Val(Nz))

    return ψ
end

Ψwn08 = overturning(wn08)
Ψqg08 = overturning(qg08)
Ψup08 = overturning(up08)
Ψeb08 = overturning(eb08)
Ψom08 = overturning(om08)
Ψwd08 = overturning(wd08)
Ψqg16 = overturning(qg16)
Ψup16 = overturning(up16)
Ψom16 = overturning(om16)
Ψwd16 = overturning(wd16)
Ψwn16 = overturning(wn16)
Ψqg32 = overturning(qg32)
Ψwn32 = overturning(wn32)

