
"""
    @kwdef struct AxisSearchConfig{T <: AbstractFloat} <: GraphConstructionConfig
        kernel_σ_multiplier::T = convert(T, 3)
        make_undirected::Bool = true
        remove_null_nbs::Bool = true
        w_lb::T = convert(T, 0.1)
    end

If `remove_null_nbs` is set to false, you can infer which neighbors came from which axis direction from the generated graph.
"""
@kwdef struct AxisSearchConfig{T <: AbstractFloat} <: GraphConstructionConfig
    kernel_σ_multiplier::T = convert(T, 3)
    make_undirected::Bool = true
    remove_null_nbs::Bool = true # if set to false, you can infer which neighbors came from which axis direction.
    w_lb::T = convert(T, 0.1)
end

function axis_search(
    d::Integer,
    node::Integer,
    X::Vector{AV},
    kernel_ref::WeightKernel,
    config::AxisSearchConfig{T},
    ) where {AV <: AbstractVector, T}

    x0 = X[node]

    # find the most axis-aligned points; one for positive and one for negtaive axis directions.
    cos_evals, candidates = get_angle_candidates(
        d, X, node,
        convert(T, cos(π/4)),
    )
    #degs = rad2deg.( acos.(cos_evals)) # debug.

    rs = collect(
        X[n][d] - x0[d]
        for n in candidates
    )

    # posiitve
    inds = findall(xx->xx>0, rs)
    cs_positive = cos_evals[inds]
    can_positive = candidates[inds]

    nb_positive = 0
    if !isempty(cs_positive)

        nb_positive = pickcandidate(
            #rs_positive,
            cs_positive,
            can_positive,
            X,
            node,
            kernel_ref,
            config.kernel_σ_multiplier
        )
    end

    # negative
    inds = findall(xx->xx<0, rs) # don't include the case r == 0.
    #rs_negative = rs[inds]
    cs_negative = cos_evals[inds]
    can_negative = candidates[inds]

    nb_negative = 0
    if !isempty(cs_negative)

        nb_negative = pickcandidate(
            #rs_negative,
            cs_negative,
            can_negative,
            X,
            node,
            kernel_ref,
            config.kernel_σ_multiplier
        )
    end

    return nb_positive, nb_negative
end

function pickcandidate(
    cos_evals, candidates, X, node,
    kernel_ref, kernel_σ_multiplier,
    )

    @assert !isempty(cos_evals)
    
    _, ind0 = findmax(cos_evals)
    min_angle_node = candidates[ind0]

    dist = norm(X[min_angle_node] - X[node])
    #a = solveparam(kernel_ref, dist, kernel_σ_multiplier)
    #kernel = createkernel(kernel_ref, a)
    kernel = createkernelσ(
        kernel_ref, dist, kernel_σ_multiplier,
    )

    score = collect(
        #cos_evals[i]*evalkernel(kernel, abs(rs[i]))
        cos_evals[i]*evalkernel(
            kernel,
            norm(X[node]-X[candidates[i]]),
        )
        for i in eachindex(candidates)
    )

    _, ind = findmax(score)

    # # debug
    # k_evals = collect(
    #     evalkernel(kernel, norm(X[node]-X[n]))
    #     for n in candidates
    # )
    # println("[candidates cos_evals k_evals score]")
    # display([candidates cos_evals k_evals score])
    # println("Initial: ", candidates[ind0])
    # println("Picked: ", candidates[ind])
    # println()

    return candidates[ind]
end

function axiscosineangle(d, xn, x0)
    p = copy(x0)
    p[d] = xn[d]

    a = p-x0
    b = xn-x0
    return dot(a,b)/(norm(a)*norm(b))
end

function get_angle_candidates(
    d::Integer, X::Vector{AV}, node::Integer, lb_threshold
    ) where {T, AV <: AbstractVector{T}}

    @assert 1 <= node <= length(X)

    cos_evals = collect(
        axiscosineangle(d, X[n], X[node])
        for n in eachindex(X)
    )
    inds = findall(
        xx->(isfinite(xx) && xx >= lb_threshold),
        cos_evals,
    )

    return cos_evals[inds], inds
end

function rad2deg(x::T)::T where T <: AbstractFloat
    return x * 180 / convert(T, pi)
end

#### graph.

# find neighbours for one node. # assumes X is non-empty.
function findaxisnbs(
    X::Vector{AV},
    node::Integer,
    kernel_ref::WeightKernel,
    config::AxisSearchConfig,
    )  where AV <: AbstractVector

    return collect(
        Iterators.flatten(
            axis_search(
                d, node, X, kernel_ref, config,
            ) for d in eachindex(X[begin])
        )
    )
end

# find neighbours for all nodes.
function findaxisnbs(
    config::AxisSearchConfig,
    X::Vector{AV};
    kernel_ref = SqExp(T),
    ) where {T, AV <: AbstractVector{T}}

    @assert !isempty(X)

    nbs = collect(
        findaxisnbs(X, node, kernel_ref, config)
        for node in eachindex(X)
    )

    if config.remove_null_nbs
        for n in eachindex(nbs)
            filter!(xx->(xx != 0), nbs[n])
        end
    end

    return nbs
end

# make the graph symmetric. nbs is a list of 1-hop neighbours for each node.
function symconnections(nbs::Vector{Vector{Int}})
    
    M0 = sum( length(nbs[n]) for n in eachindex(nbs) )
    N_edges_ub = M0*2
    edges = Vector{Tuple{Int,Int}}(undef, N_edges_ub)

    k = 0
    for i in eachindex(nbs)
        for j in nbs[i]
            
            k += 1
            edges[k] = (i,j)

            k += 1
            edges[k] = (j,i)
        end
    end
    resize!(edges, k)
    unique!(xx->lexorder(xx), edges) # O(N^2) run-time. TODO see if we can do this without unique().
    # there are also non-unique entries for 45 degree stations, so unique() here is a simple implementation.

    srcs = map(xx->xx[begin], edges)
    dests = map(xx->xx[begin+1], edges)
    return srcs, dests
end

function lexorder(s::Tuple{Int,Int})::Tuple{Int,Int}
    if s[begin+1] > s[begin]
        return (s[begin+1], s[begin])
    end
    return s
end

function nbs2connections(nbs::Vector{Vector{Int}})
    dests = collect( Iterators.flatten(nbs) )
    srcs = Vector{Int}(undef, length(dests))
    
    k = 0
    for n in eachindex(nbs)
        for _ in eachindex(nbs)
            k += 1
            srcs[k] = n
        end
    end
    return srcs, dests
end

function computedists(
    metric::Distances.Metric,
    srcs::Vector{Int},
    dests::Vector{Int},
    X,
    )

    @assert length(srcs) == length(dests)
    return collect(
        Distances.evaluate(metric, X[srcs[m]], X[dests[m]])
        for m in eachindex(srcs)
    )
end


"""
    create_axis_graph(
        config::AxisSearchConfig{T},
        X::Vector{AV},
    ) where {T, AV <: AbstractVector{T}}

Constructs an undirected graph via the axis-search algorithm.
"""
function create_axis_graph(
    config::AxisSearchConfig{T},
    X::Vector{AV},
    ) where {T, AV <: AbstractVector{T}}

    return create_axis_graph(
        config,
        Distances.Euclidean(),
        X,
        SqExp(T),
    )
end

function create_axis_graph(
    config::AxisSearchConfig,
    metric::Distances.Metric,
    X::Vector{AV},
    weight_kernel_ref,
    ) where AV <: AbstractVector

    w_lb = config.w_lb
    
    nbs = findaxisnbs(config, X)
    srcs, dests = symconnections(nbs)
    dists = computedists(metric, srcs, dests, X)
    
    max_dist = maximum(dists)
    kernel = solve_createkernel(weight_kernel_ref, max_dist, w_lb)
    ws = map(xx->evalkernel(kernel, xx), dists)

    edges = EdgeList(srcs, dests, ws)
    #A = sparse(srcs, dests, ws) # version for directed graph
    N_nodes = length(X)
    rs, cs, as = srcs, dests, ws
    A = sparse(
        vcat(rs, cs), vcat(cs, rs), vcat(as, as), N_nodes, N_nodes,
    )
    graph = Graphs.SimpleGraph(A)

    return UndirectedGraph(graph, edges)
end
