"""
    @kwdef struct WarpConfig{T <: AbstractFloat}
        N_bands::Int = 20
        λ_max::T = 2.0 # in (0, 2]
        iterate_threshold::T = convert(T, 10)
        discount_factor::T = convert(T, 0.5)
        max_iters::Int = 10
        aggregate_option::Symbol = :sum_normalized
        verbose::Bool = false
    end

`aggregate_option` can be `:sum_normalized` or `:sum`
"""
@kwdef struct WarpConfig{T <: AbstractFloat}
    N_bands::Int = 20
    λ_max::T = 2.0 # in (0, 2]
    iterate_threshold::T = convert(T, 10)
    discount_factor::T = convert(T, 0.5)
    max_iters::Int = 10
    aggregate_option::Symbol = :sum_normalized
    verbose::Bool = false
end


# for use with dimensional expansion kernels.
"""
    get_grid_warp_samples(
        x_nD::AbstractArray{T},
        config::WarpConfig,
    ) where T <: Real

Returns:
- `W`, an array of the same size as `x_nD` that contains the warp samples.
- `G`, the unit graph associated with `x_nD`.
"""
function get_grid_warp_samples(
    x_nD::AbstractArray{T},
    config::WarpConfig) where T <: Real
    
    # grid graph, no wrap-around, undirected.
    nbs = getgridnbs(size(x_nD))
    G = UnitGrid(nbs)

    TL = create_rwlaplacian(T, G)
    return get_warp_samples(vec(x_nD), TL, config), G
end

"""
    get_knn_warp_samples(
        X::AbstractVector,
        y::AbstractVector{T},
        config::WarpConfig,
        knn_config::KNNConfig;
        check_finite::Bool = true,
    ) where T <: Real

To get an undirected `G` output, make sure `knn_config.make_symmetric` is `true`.

Returns:
- `W`, warp samples for dataset `(X,y)`.
- `G`, the k-nearest neighbors graph for the dataset.
"""
function get_knn_warp_samples(
    X::AbstractVector,
    y::AbstractVector{T},
    config::WarpConfig,
    knn_config::KNNConfig;
    check_finite::Bool = true,
    ) where T <: Real
    
    G = create_knn_graph(knn_config, X)

    TL = create_rwlaplacian(G)
    if check_finite
        if !isfinite(minimum(TL)) || !isfinite(maximum(TL))
            error("The random-walk Laplacian has non-finite values. One possible solution is to change the graph weights, and try again.")
        end
    end

    return get_warp_samples(y, TL, config), G
end

"""
    get_warp_samples(
        x::Vector{T},
        TL::Union{Matrix, SparseMatrixCSC},
        config::WarpConfig,
    ) where T <: AbstractFloat

`TL` is the fundamental one-hop operator. One choice is the random-walk Laplacian of a graph.

`x` is the graph signal from which we derive the warp samples from.

Returns:
- `W`, warp samples for dataset `(X,y)`.
"""
function get_warp_samples(
    x::Vector{T},
    TL::Union{Matrix, SparseMatrixCSC},
    config::WarpConfig,
    ) where T <: AbstractFloat

    Ws, _ = computewarpsamplesiterated(
        x,
        TL,
        BernsteinFBConfig(
            N_bands = config.N_bands,
            λ_max = config.λ_max,
        );
        iterate_threshold = config.iterate_threshold,
        discount_factor = config.discount_factor,
        max_iters = config.max_iters,
        verbose = config.verbose,
    )

    W_out = sum(Ws)
    if config.aggregate_option == :sum_normalized
        W_out = sumnormalized(Ws)
    end

    return W_out
end