# sparse matrix routines

##### for unit weight, undirected grid graphs.

"""
    create_degree(g::UndirectedGraph)

Returns the degree matrix of graph `g`.
"""
function create_degree(g::UnitGrid)
    graph = getgraph(g)
    return spdiagm(Graphs.degree(graph))
end

function create_invdegree(::Type{T}, g::UnitGrid) where T
    
    deg = create_degree(g)
    return spdiagm( collect( one(T) ./ diag(deg)) )
end

function create_sqrtinvdegree(::Type{T}, g::UnitGrid) where T
    
    inv_deg = create_invdegree(T, g)
    map!(sqrt, inv_deg, inv_deg)
    return inv_deg
end

# random walk Laplacian
"""
    create_rwlaplacian(::Type{T}, g::UnitGrid) where T

Returns the random-walk Laplacian of a unit grid graph, with scalar data type `T`.
"""
function create_rwlaplacian(::Type{T}, g::UnitGrid) where T

    inv_deg = create_invdegree(T, g)
    L = create_laplacian(g)
    
    return inv_deg*L # pg 99.
end

"""
    create_snlaplacian(::Type{T}, g::UnitGrid) where T

Returns the symmetric normalized Laplacian of a unit grid graph, with scalar data type `T`.
"""
function create_snlaplacian(::Type{T}, g::UnitGrid) where T

    deg_inv_sqrt = create_sqrtinvdegree(T, g)
    L = create_laplacian(g)
    
    return deg_inv_sqrt*L*deg_inv_sqrt
end

"""
    create_adjacency(g::UnitGrid)

Returns the adjacency matrix of `g`.
"""
function create_adjacency(g::UnitGrid)
    graph = getgraph(g)
    return Graphs.adjacency_matrix(graph)
end


##### for weighted undirect graphs

# look at mats.jl from RealGSP.jl
"""
    create_adjacency(g::UndirectedGraph)
"""
function create_adjacency(g::UndirectedGraph)
    
    P = getedgelist(g)
    as, rs, cs = P.ws, P.srcs, P.dests

    N = getNnodes(g)
    return sparse(
        vcat(rs, cs), vcat(cs, rs), vcat(as, as), N, N,
    )
    #return sparse(rs, cs, as, N, N) # version for directed graph
end

function create_degreediag(g::UndirectedGraph)
    A = create_adjacency(g)
    
    return vec(reduce(+, A, dims = 2)) # same as A*ones(T, N)
end

"""
    create_degree(g::UndirectedGraph)
"""
function create_degree(g::UndirectedGraph{T}) where T <: AbstractFloat
        return spdiagm(create_degreediag(g))
end

function create_invdegree(g::UndirectedGraph)

    tmp = create_degreediag(g)
    v = collect( 1/tmp[i] for i in eachindex(tmp) ) 
    return spdiagm(v)
end

function create_sqrtinvdegree(g::UndirectedGraph)

    tmp = create_degreediag(g)
    v = collect( sqrt(1/tmp[i]) for i in eachindex(tmp) ) 
    return spdiagm(v)
end

# random walk Laplacian
"""
create_rwlaplacian(g::UndirectedGraph)
"""
function create_rwlaplacian(g::UndirectedGraph)

    inv_deg = create_invdegree(g)
    L = create_laplacian(g)
    
    return inv_deg*L # pg 99.
end

"""
create_snlaplacian(g::UndirectedGraph)
"""
function create_snlaplacian(g::UndirectedGraph)

    deg_inv_sqrt = create_sqrtinvdegree(g)
    L = create_laplacian(g)
    
    return deg_inv_sqrt*L*deg_inv_sqrt
end

##### common graph routines.

"""
    create_laplacian(g::UndirectedGraph)

Returns the combinatorial Laplacian of `g`.
"""
function create_laplacian(g::GraphContainer)

    A = create_adjacency(g)
    deg = create_degree(g)
    L = deg - A
    return L
end

