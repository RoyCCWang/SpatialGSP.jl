

# knn results to 1-hop neighbours for each node.
function getnbs_knn(X::AbstractVector{AV}, k::Integer; single_distance = false, atol = 1e-6) where AV <: AbstractVector
    @assert !isempty(X)

    X_mat = reshape(collect(Iterators.flatten(X)), length(X[begin]), length(X))
    return getnbs_knn(X_mat, k; single_distance = single_distance, atol = atol)
end

function getnbs_knn(X_mat::Matrix{T}, k::Integer; single_distance = false, atol = eps(T)*100) where T <: AbstractFloat
    
    knn_metric = Distances.Euclidean()
    tree = NN.KDTree(X_mat, knn_metric)
    
    # results are sorted.
    idxs, dists = NN.knn(tree, X_mat, k, true)

    # prune the current node.
    for n in eachindex(idxs)
        deleteat!(idxs[n], 1)
        deleteat!(dists[n], 1)
    end

    if single_distance
        for n in eachindex(dists)
            keep_inds = findall(
                xx->isapprox(
                    xx, dists[begin]; atol = atol,
                ),
                dists[n]
            )
            dists[n] = dists[n][keep_inds]
            idxs[n] = idxs[n][keep_inds]
        end
    end

    return idxs, dists
end

"""
    @kwdef struct KNNConfig{T <: AbstractFloat} <: GraphConstructionConfig
        single_distance::Bool = false
        atol::T = convert(T, 1e-3)
        normalize_min_weight::Bool = false
        make_symmetric::Bool = true
        k::Int = 4
        kernel_gain = convert(T, Inf)
        kernel_σ² = convert(T, Inf)
    end

If `normalize_min_wieght` is `false`, then the minimum weight is 1.
"""
@kwdef struct KNNConfig{T <: AbstractFloat} <: GraphConstructionConfig
    single_distance::Bool = false
    atol::T = convert(T, 1e-3)
    normalize_min_weight::Bool = false # if true, the minimum weight is 1.
    make_symmetric::Bool = true
    k::Int = 4
    kernel_gain = convert(T, Inf)
    kernel_σ² = convert(T, Inf)
end

function create_knn_graph(
    config::KNNConfig{T},
    X::Vector{AV},
    ) where {T <: AbstractFloat, AV <: AbstractVector}

    k = config.k
    nbs, dists = getnbs_knn(X, k+1; single_distance = config.single_distance, atol = config.atol)

    return create_knn_graph(config, nbs, dists)
end

function create_knn_graph(
    config::KNNConfig{T},
    nbs,
    dists::Vector{Vector{T}},
    ) where {T <: AbstractFloat}
    
    @assert length(nbs) == length(dists)

    gain, σ² = config.kernel_gain, config.kernel_σ²
    kernel = SqExpGain(1/(2*σ²), gain)
    if !isfinite(σ²) || !isfinite(gain)
        min_dist = minimum(Iterators.flatten(dists))
        gain = convert(T, exp(1/2))
        σ² = min_dist
        kernel = SqExpGain(1/(2*σ²), gain)
    end

    edges_tuple, ws = nbs_2_edges(nbs, dists, kernel, config.make_symmetric)
    if config.normalize_min_weight
        Z = minimum( abs(w) for w in ws )
        ws = ws ./ Z
    end

    srcs = map(xx->xx[begin], edges_tuple)
    dests = map(xx->xx[begin+1], edges_tuple)
    edges = EdgeList(srcs, dests, ws)

    #A = sparse(srcs, dests, ws) # version for directed graph
    N_nodes = length(nbs)
    rs, cs, as = srcs, dests, ws
    A = sparse(
        vcat(rs, cs), vcat(cs, rs), vcat(as, as), N_nodes, N_nodes,
    )
    graph = Graphs.SimpleGraph(A)

    return UndirectedGraph(graph, edges)
end

function nbs_2_edges(
    nbs::Vector{Vector{Int}},
    dists::Vector{Vector{T}},
    kernel,
    make_symmetric::Bool,
    ) where T <: AbstractFloat

    @assert length(nbs) == length(dists)

    edges = Vector{Tuple{Int,Int}}(undef, 0)
    ws = Vector{T}(undef, 0)

    for n in eachindex(nbs)
        for i in eachindex(nbs[n])
            push!(edges, (n, nbs[n][i]))
            push!(ws, evalkernel(kernel, dists[n][i]))
        end
    end
    tmp = unique(xx->xx[begin], zip(edges, ws))
    edges = map(xx->xx[begin], tmp)
    ws = map(xx->xx[begin+1], tmp)

    if make_symmetric
        for i in eachindex(tmp)
            edge = tmp[i][begin]
            w = tmp[i][begin+1]
            push!(edges, (edge[begin+1], edge[begin]) )
            push!(ws, w)
        end
        tmp = unique(xx->xx[begin], zip(edges, ws))
        edges = map(xx->xx[begin], tmp)
        ws = map(xx->xx[begin+1], tmp)
    end

    return edges, ws
end
