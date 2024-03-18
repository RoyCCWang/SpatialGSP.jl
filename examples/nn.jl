
# using LinearAlgebra
# import Random
# Random.seed!(25)

# import PlotlyLight as PLY
# import ColorSchemes

# PLT.close("all")
# fig_num = 1

# D = 3
# T = Float64
# #N = 10
# N = 1_000_000

# X = collect( randn(T, D) for _ = 1:N )

# X_mat = reshape(collect(Iterators.flatten(X)), D, N)

# # set up.
# knn_metric = Distances.Euclidean()
# println("Timing: NN.KDTree(), N = $N, D = $D")
# @time tree = NN.KDTree(X_mat, knn_metric)
# # 9 sec
# # results are sorted.
# k_knn = 2*D+1
# println("Timing: NN.knn(), k_knn = $k_knn")
# @time idxs, dists = NN.knn(tree, X_mat, k_knn, true)

# """
# Timing: NN.KDTree(), N = 10000000, D = 3
#   9.610241 seconds (31 allocations: 625.612 MiB, 8.75% gc time)
# Timing: NN.knn(), k_knn = 7
#  34.906567 seconds (20.00 M allocations: 2.235 GiB, 10.98% gc time)
# """


include("./helpers/image.jl")
include("./helpers/viz/misc.jl")

Random.seed!(25)

fig_num = 1
PLT.close("all")

T = Float64
#atol = convert(T, 1e-6)

load_folder = joinpath(
    homedir(),
    "work/GSP/"
)

D = 2
project_name = "parrot"
img = loadkodakimage(T, "./data/kodim23.png"; discard_pixels = 1)
x_nD, x_ranges = image2samples(img)
sz_x = size(x_nD)

nbs = GSP.getgridnbs(size(x_nD))
x = vec(x_nD)

# generate the sampling locations.

#N = 10
N = length(img)

X = collect(
    GSP.convertcompactdomain(
        rand(T, D), zeros(T, D), ones(T, D),
        ones(T, D), size(img) .* ones(T, D),
    )
    for _ = 1:N
)

# itp.
import Interpolations
Xrs = (1:size(img,1), 1:size(img,2))
itp = Interpolations.interpolate(
    img,
    Interpolations.BSpline( 
        Interpolations.Cubic(    
            Interpolations.Line(Interpolations.OnGrid()),
        ),
    ),
)
#itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
scaled_itp = Interpolations.scale(itp, Xrs...)
etp = Interpolations.extrapolate(scaled_itp, zero(T))

etp_X = collect( etp(x...) for x in X )

# ## visualize graph signal.
fig_num = VIZ.plotmeshgrid2D(
    PLT,
    x_ranges,
    x_nD,
    [],
    "x",
    fig_num,
    "original image, x";
    cmap = "Greys_r",
    matrix_mode = true
)

# 
k_knn = 2*D+1


nbs, dists = GSP.getnbs_knn(
    X, k_knn;
    single_distance = false
)

knn_config = GSP.KNNConfig{T}(
    k = k_knn,
    kernel_gain = 1.23,
    kernel_σ² = 2.3,
)
kernel = GSP.SqExp(T)
println("Timing: create_knn_graph")
@time G = GSP.create_knn_graph(knn_config, X)

warp_config = GSP.WarpConfig{T}()
println("get_knn_warp_samples()")
@time W, G = GSP.get_knn_warp_samples(
    X, vec(img), warp_config, knn_config,
)

println("get_grid_warp_samples()")
@time W_grid, G_grid = GSP.get_grid_warp_samples(x_nD, warp_config)

# ## visualize graph signal.
fig_num = VIZ.plotmeshgrid2D(
    PLT,
    x_ranges,
    reshape(W_grid, size(x_nD)),
    [],
    "x",
    fig_num,
    "warp samples, grid";
    cmap = "bwr",
    matrix_mode = true
)

nothing