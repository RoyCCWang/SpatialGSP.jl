# for 2D spatial graphs.

# for coordinate d.
function getedgecoord(
    G::UndirectedGraph, X::Vector{AV}, d::Integer,
    ) where {T, AV <: AbstractVector{T}}
    
    edges = getedgelist(G)
    srcs, dests = edges.srcs, edges.dests
    return getedgecoord(srcs, dests, X, d)
end

function getedgecoord(
    srcs::Vector{Int},
    dests::Vector{Int},
    X::Vector{AV},
    d::Integer,
    ) where {T, AV <: AbstractVector{T}}
    
    # parse & checks.
    @assert !isempty(X)
    @assert 1 <= d <= length(X[begin])
    @assert length(srcs) == length(dests)
    
    # 
    k = 0
    out = Vector{Any}(undef, length(srcs)*3)
    for j in eachindex(srcs)

        k += 1
        out[k] = X[srcs[j]][begin+d-1]

        k += 1
        out[k] = X[dests[j]][begin+d-1]

        # we need a `null` JavaScript value to break up each line segment.
        k+=1
        out[k] = nothing # forces the generated HTML or JSON to use a null value.
    end
    resize!(out, k)

    return out
end

