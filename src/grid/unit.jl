@kwdef struct UnitGraphConfig{UT <: Tuple}
    sz::UT
end

abstract type GraphDirectionality end
struct Undirected <: GraphDirectionality end
struct Directed <: GraphDirectionality end

function getshiftCIs(::Directed, ::NTuple{D,Int}) where D
    shifts_up = collect( CartesianIndex( ntuple(xx->(xx==d ? 1 : 0), D) ) for d = 1:D )
    return shifts_up
end

function getshiftCIs(::Undirected, ::NTuple{D,Int}) where D
    shifts_down = collect( CartesianIndex( ntuple(xx->(xx==d ? -1 : 0), D) ) for d = 1:D )
    shifts_up = collect( CartesianIndex( ntuple(xx->(xx==d ? 1 : 0), D) ) for d = 1:D )
    shifts = vcat(shifts_down, shifts_up)

    return shifts
end

# unweighted.
struct UnitGrid{D, GT} <: GraphContainer
    graph::GT

    # 1-hop neighbors only! not 0-hop (self-loops) nor k-hops, k >= 2.
    nbs::Array{Vector{CartesianIndex{D}}, D}
end

function UnitGrid(nbs)
    return UnitGrid(unitgraph_grid(nbs), nbs)
end

function getnbs(G::UnitGrid)
    return G.nbs
end

# get 1-hop neighbours. undirected, no wrap-around.
function getgridnbs(
    sz_A::NTuple{D,Int},
    ) where D

    # no degenerate dimensions allowed. This avoids self-loops.
    for d in eachindex(sz_A)
        @assert sz_A[d] > 1
    end
    
    N_nodes = prod(sz_A)
    shifts = getshiftCIs(Undirected(), sz_A)

    # one-hop neighbors.
    CIs = CartesianIndices(sz_A)
    @assert length(CIs) == N_nodes

    # no wrap-around means the neighbhours must lie in CIs.
    nb1s = collect(
    collect(
            Iterators.filter(xx->(xx in CIs), n_CI .+ shifts)
        )
        for n_CI in CIs
    )

    return nb1s
end

# Does not check for self-loops.
function unitgraph_grid(
    nbs::Array{Vector{CartesianIndex{D}},D}
    ) where D

    g = Graphs.SimpleGraph(length(nbs))

    #edge_list = Vector{Tuple{Int,Int}}(undef, 0)
    LIs = LinearIndices(nbs)

    for src in LIs

        for ci in nbs[src]
            dest = LIs[ci]

            Graphs.add_edge!(g, src, dest)

            #push!(edge_list, (src, dest)) # debug
            # if src != dest # avoid self-loop
            #     Graphs.add_edge!(g, src, dest)
            # end
        end
    end

    return g#, edge_list
end

# front end.

function creategraph(config::UnitGraphConfig, args...)
    nbs = getgridnbs(config.sz)
    return UnitGrid(nbs)
end