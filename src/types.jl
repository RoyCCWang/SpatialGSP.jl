
# Unit pre-fix means unweighted, i.e. all edges have unit weights, i.e. a value of 1.

abstract type WeightKernel end

function createkernel(::KT, args...) where KT <: WeightKernel
    return KT(args...)
end

struct SqExpGain{T <: AbstractFloat} <: WeightKernel
    a::T # 1/(2*π^2)
    gain::T
end

function evalkernel(θ::SqExpGain{T}, x::T)::T where T
    return θ.gain*exp(-θ.a*x^2)
end

struct SqExp{T <: AbstractFloat} <: WeightKernel
    a::T # 1/(2*π^2)
end

function SqExp(::Type{T}) where T <: AbstractFloat
    return SqExp(one(T))
end

function evalkernel(θ::SqExp{T}, x::T)::T where T
    return exp(-θ.a*x^2)
end

# y = exp(-a*x^2), given (x,y), solves for a.
function solveparam(::SqExp, x::T, y::T)::T where T
    return -log(y)/x^2
end

function solve_createkernel(θ::WeightKernel, x, y)
    #
    ps = solveparam(θ, x, y)
    return createkernel(θ, ps...)
end

function createkernelσ(::SqExp, x, k)
    #
    #ps = solveparam(θ, x, y)
    σ = x/k
    a = 1/(2*σ^2)
    return SqExp(a)
end

# struct RQ32Weight{T} <: WeightKernel
#     a::T
# end

# function evalkernel(θ::RQ32Weight{T}, x::T)::T where T
#     tmp = 1 + θ.a*x^2
#     return 1/(tmp^(3/2))
# end


# function solveparam(::RQ32Weight, x::T, y::T)::T where T
#     B = (1/y)^(3/2)
#     return (B-1)/(x^2)
# end

#####

abstract type GraphContainer end

function getNnodes(G::GraphContainer)
    return Graphs.nv(getgraph(G))
end

function getgraph(G::GraphContainer)
    return G.graph
end

function getNnbs(G::GraphContainer)
    graph = getgraph(G)
    N_nodes = Graphs.nv(graph)
    return [
        length(Graphs.neighbors(graph, n)) for n = 1:N_nodes
    ]
end

# # Weighted graphs.
# see knn.jl, nnk.jl for construction.

struct EdgeList{T}
    # all of the same length.
    srcs::Vector{Int}
    dests::Vector{Int}
    ws::Vector{T}
end

function getsrcs(E::EdgeList)
    return E.srcs
end

function getdests(E::EdgeList)
    return E.dests
end

function getweights(E::EdgeList)
    return E.ws
end


# The edges only contain one direction, not both, and create_adjacency() will double up the edges when forming the adjacency matrix.
struct UndirectedGraph{T <: AbstractFloat, GT} <: GraphContainer
    #
    graph::GT
    #nbs::Vector{Vector{Int}} # 1-hop neighbors.

    edges::EdgeList{T}  # note the edge iterator for graph and the ordering in `edges` are different in general.
end

function getnbs(G::UndirectedGraph)
    graph = getgraph(G)
    N_nodes = Graphs.nv(graph)

    return collect(
        Graphs.neighbors(graph, n)
        for n = 1:N_nodes
    )
end

function getedgelist(G::UndirectedGraph)
    return G.edges
end


abstract type GraphConstructionConfig end
