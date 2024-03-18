



function dcKDtree(
    X::Vector{Vector{T}},
    metric,
    batch_ranges,
    save_path::String,
    ) where T <: Number
    
    @assert batch_ranges[end][end] == length(X)
    @assert batch_ranges[begin][begin] == 1

    tree = NN.KDTree(hcat(X...), metric)
    serialize(save_path, (tree, X, batch_ranges))

    return nothing
end


function dcknn(
    batch_ind::Integer,
    k::Integer,
    load_path::String,
    save_path_prefix::String,
    )

    tree, X, batch_ranges = deserialize(load_path)
    pts_mat = hcat(X[batch_ranges[batch_ind]]...)
    inds, dists = NN.knn(tree, pts_mat, k, false)

    serialize(
        "$(save_path_prefix)$(batch_ind)",
        (inds, dists),
    )
    return nothing
end

function getknngc(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    batch_ranges,
    metric,
    k::Integer,
    tmp_folder::String;
    delete_files = false,
    ) where T <: AbstractFloat

    inds, dists = getknn(X, batch_ranges, metric, k, tmp_folder; delete_files = delete_files)
    Base.GC.gc()

    return inds, dists
end

function getknn(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    batch_ranges,
    metric,
    k::Integer,
    tmp_folder::String;
    delete_files = false,
    ) where T <: AbstractFloat

    # set up.
    tree_path = joinpath(tmp_folder, "kd_tree")
    dcKDtree(X, metric, batch_ranges, tree_path)

    # distributed computing.
    save_path_prefix = joinpath(tmp_folder, "knn_")
    processfunc = nn->dcknn(nn, k, tree_path, save_path_prefix)
    
    # processfunc(1)
    # @assert 4==3
    N_batches = length(batch_ranges)
    @sync begin
        pmap(processfunc, 1:N_batches)
    end

    # assemble the results.
    inds, dists = assembledcknn(T, save_path_prefix, N_batches; delete_files = delete_files)

    if delete_files
        rm(tree_path)
    end

    return inds, dists
end



# not type stable, since we don't know what the type of ps and f_ps.
function assembledcknn(
    ::Type{T},
    load_prefix::String,
    N_batches::Int;
    delete_files = false) where T <: AbstractFloat

    dics = collect(
        deserialize("$(load_prefix)$(n)")
        for n = 1:N_batches
    )

    # solutions.
    inds_raw = collect( Iterators.flatten( collect( dics[n][begin] for n in eachindex(dics) ) ) )
    dists_raw = collect( Iterators.flatten( collect( dics[n][begin+1] for n in eachindex(dics) ) ) )
    
    inds = convert(Vector{Vector{Int}}, inds_raw)
    dists = convert(Vector{Vector{T}}, dists_raw)

    if delete_files
        # delete files.
        for n = 1:N_batches
            t = @task begin
                rm(joinpath(load_prefix, "$(load_prefix)$(n)"))
            end
            schedule(t)
            wait(t)
        end
    end
    
    return inds, dists
end


# post process NN.knn()
function removeselfloop!(knn_inds::Vector{Vector{Int}}, dists)

    for n in eachindex(knn_inds)
        
        del_inds = findall(xx->xx==n, knn_inds[n])
        sort!(del_inds)
        deleteat!(knn_inds[n], del_inds)
        deleteat!(dists[n], del_inds)
    end

    return nothing
end