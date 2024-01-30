using Oceananigans
using BaroclinicAdjustment
using BaroclinicAdjustment.Diagnostics
using Statistics: mean

weno_series = all_fieldtimeseries("weno9pV_sixteen_new_snapshots.jld2")
cent_series = all_fieldtimeseries("weno9pV_sixteen_center_snapshots.jld2")
vert_series = all_fieldtimeseries("weno9pV_sixteen_center_vertical_snapshots.jld2")
buff_series = all_fieldtimeseries("weno9pV_sixteen_center_buffer_snapshots.jld2")

bw = weno_series[:b]
bc = cent_series[:b]
bv = vert_series[:b]
bb = buff_series[:b]

iters = 50:201
Niter = length(iters)

Bw = mean(bw[iters[1]], dims = 1)
Bc = mean(bc[iters[1]], dims = 1)
Bv = mean(bv[iters[1]], dims = 1)
Bb = mean(bb[iters[1]], dims = 1)

for i in iters[2:end]
    @info "doing iteration $i"
    Bw .+= mean(bw[i], dims = 1)
    Bc .+= mean(bc[i], dims = 1)
    Bv .+= mean(bv[i], dims = 1)
    Bb .+= mean(bb[i], dims = 1)
end

Bw ./= Niter
Bc ./= Niter
Bv ./= Niter
Bb ./= Niter

Vccc = BaroclinicAdjustment.Diagnostics.VolumeField(bw.grid)

Vdomain = sum(Vccc)

# b²w = propagate(bw, Vccc; func = (x, y) -> x^2 * y / Vdomain)
# b²c = propagate(bc, Vccc; func = (x, y) -> x^2 * y / Vdomain)
# b²v = propagate(bv, Vccc; func = (x, y) -> x^2 * y / Vdomain)
# b²b = propagate(bb, Vccc; func = (x, y) -> x^2 * y / Vdomain)

# b2w = []
# b2c = []
# b2v = []
# b2b = []

# for i in 1:201
#     @info "doing iteration $i"
#     push!(b2w, sum(interior(b²w[i])))
#     push!(b2c, sum(interior(b²c[i])))
#     push!(b2v, sum(interior(b²v[i])))
#     push!(b2b, sum(interior(b²b[i])))
# end

# b2w = Float64.(b2w)
# b2c = Float64.(b2c)
# b2v = Float64.(b2v)
# b2b = Float64.(b2b)

N²w = propagate(bw, Vccc; func = (x, y) -> ∂z(x) * y / Vdomain)
N²c = propagate(bc, Vccc; func = (x, y) -> ∂z(x) * y / Vdomain)
N²v = propagate(bv, Vccc; func = (x, y) -> ∂z(x) * y / Vdomain)
N²b = propagate(bb, Vccc; func = (x, y) -> ∂z(x) * y / Vdomain)

N2w = Float64[]
N2c = Float64[]
N2v = Float64[]
N2b = Float64[]

for i in 1:201
    @info "doing iteration $i"
    push!(N2w, sum(interior(N²w[i], :, :, 2:49)))
    push!(N2c, sum(interior(N²c[i], :, :, 2:49)))
    push!(N2v, sum(interior(N²v[i], :, :, 2:49)))
    push!(N2b, sum(interior(N²b[i], :, :, 2:49)))
end
