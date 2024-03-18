module SpatialGSP

using LinearAlgebra
using SparseArrays
using Statistics

import Graphs
import Distances
import NearestNeighbors as NN

# constant values.
function pidiv4(::Type{T})::T where T <: AbstractFloat
    return convert(T, π/4)
end

function pidiv2(::Type{T})::T where T <: AbstractFloat
    return convert(T, π/2)
end

function twopi(::Type{T})::T where T <: AbstractFloat
    return convert(T, 2*π)
end

function onepi(::Type{T})::T where T <: AbstractFloat
    return convert(T, π)
end

include("types.jl")

include("grid/unit.jl")
include("spatial/knn.jl")
#include("spatial/LLE.jl")
include("spatial/axis_search.jl")
include("spatial/viz.jl")

include("matrices.jl")

include("filterbank/bernstein.jl")

include("frontend/warpmap.jl")

include("utils.jl")

export WarpConfig,
get_warp_samples,
get_grid_warp_samples,

create_adjacency,
create_rwlaplacian,
create_snlaplacian,
create_degree,
create_laplacian,

AxisSearchConfig,
create_axis_graph,

KNNConfig,
get_knn_warp_samples,

convertcompactdomain

end # module SpatialGSP
