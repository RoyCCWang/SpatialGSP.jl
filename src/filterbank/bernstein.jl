# the usual filters from Unser/Van de Ville group, as presented in Leonardi 2013 (DOI: 10.1109/TSP.2013.2259825).
# These filters instantiates a tight frame, if they are implemented exactly in the frequency domain.
# for node-domain implementation via polynomials, some approximation error will occur for reconstruction.
# For now, we just implement the forward transform.
# The implementation here is based on Leronardi's article, with minor modification on parameterization.


##### evalate polynomial.

# have numerical issues.
function evalbernsteinpolynomialdirect(i::Integer, n::Integer, x::T) where T <: Real
    @assert 0 <= i <= n
    return binomial(n,i)*x^i * (1-x)^(n-i)
end

# based on de Casteligau's algorithm.
function evalallbernsteins(N::Integer, x::T) where T <: Real

    if N == 0
        return collect( ones(T,1) for _ = 1:1)
    end

    one_minus_x = 1-x

    Bs_x = collect(
        zeros(T, n+1)
        for n = 0:N
    )
    
    Bs_x[begin][begin] = one(T)

    Bs_x[begin+1][begin] = one_minus_x
    Bs_x[begin+1][end] = x

    for n = 2:N
        
        B = Bs_x[n-1+begin]

        # special cases.
        Bs_x[n+begin][begin] = B[begin]*one_minus_x
        Bs_x[n+begin][end] = B[end]*x

        # use recurrance
        for k = 1:n-1
            Bs_x[n+begin][k+begin] = B[k-1+begin]*x + B[k+begin]*one_minus_x
        end
    end

    return Bs_x
end

# x ∈ [lb, ub]. Affine map to [0,1], then evaluate the Bernstein polynomials.
function evalallbernsteins(N::Integer, x::T, λ_max::T) where T <: Real
    t = clamp(x/λ_max, zero(T), one(T))
    return evalallbernsteins(N, t)
end

#### 

# This is the version used for node-domain filtering.
function evalallbernsteins(
    N::Integer,
    x::Union{Matrix, SparseMatrixCSC},
    v::Vector{T},
    ) where T <: Real

    if N == 0
        return v
    end

    buffer = similar(v)

    one_minus_x = LinearAlgebra.I - x

    Bs_x = collect(
        collect(
            Vector{T}(undef, length(v)) for _ = 1:n+1
        ) for n = 0:N
    ) # [Bernstein degree][Berstein index][filtered output index]
    
    Bs_x[begin][begin] = v

    #Bs_x[begin+1][begin] = one_minus_x
    mul!(Bs_x[begin+1][begin], one_minus_x, v)
    
    #Bs_x[begin+1][end] = x
    mul!(Bs_x[begin+1][end], x, v)

    for n = 2:N
        
        B = Bs_x[n-1+begin]

        # special cases.
        #Bs_x[n+begin][begin] = B[begin]*one_minus_x
        mul!(Bs_x[n+begin][begin], one_minus_x, B[begin])

        #Bs_x[n+begin][end] = B[end]*x
        mul!(Bs_x[n+begin][end], x, B[end])

        # use recurrance
        for k = 1:n-1
            # #Bs_x[n+begin][k+begin] = B[k-1+begin]*x + B[k+begin]*one_minus_x
            # mul!(
            #     Bs_x[n+begin][k+begin],
            #     x,
            #     B[k-1+begin],
            # )     
            # tmp = one_minus_x*B[k+begin]
            # Bs_x[n+begin][k+begin] += tmp

            innerprocess!(
                Bs_x[n+begin][k+begin],
                buffer,
                one_minus_x,
                x,
                B[k+begin],
                B[k-1+begin],
            )
        end
    end

    return Bs_x[end]
end

# modified buffer.
function innerprocess!(
    out, buffer, one_minus_x, x, B_k, B_km1,
    )

    mul!(buffer, one_minus_x, B_k)
    mul!(out, x, B_km1)
    
    for i in eachindex(out)
        out[i] = out[i] + buffer[i]
    end
    #axpy!(1, buffer, out) # stores to out.

    return nothing
end

# function myadd!(y, x)
#     for i in eachindex(y)
#         y[i] = y[i] + x[i]
#     end
#     return nothing
# end

#### front end

@kwdef struct BernsteinFBConfig{T <: AbstractFloat}
    N_bands::Int = 10
    λ_max::T
end


function applybersteinfilterbank(
    v::Vector{T},
    L::Union{Matrix, SparseMatrixCSC},
    config::BernsteinFBConfig,
    ) where T <: AbstractFloat

    #
    λ_max = config.λ_max
    Bernstein_degree = config.N_bands - 1
    Z = L ./λ_max
    ys_ND = evalallbernsteins(Bernstein_degree, Z, v)

    return ys_ND
end

function computewarpsamples(
    v::Vector{T},
    L::Union{Matrix, SparseMatrixCSC},
    #config::BernsteinFBConfig,
    λ_max::T,
    N_bands::Integer,
    ) where T <: AbstractFloat

    Bernstein_degree = N_bands - 1
    Z = L ./λ_max
    ys_ND = evalallbernsteins(Bernstein_degree, Z, v)

    R_norm, ind = findmax( norm.(ys_ND) )
    R = ys_ND[ind]

    W = zeros(T, length(ys_ND[begin]))
    for i in eachindex(ys_ND)

        if i != ind
            W += ys_ND[i]
        end
    end

    return W, R, R_norm
end

#### iterated.

function computewarpsamplesiterated(
    v::Vector{T},
    L::Union{Matrix, SparseMatrixCSC},
    config::BernsteinFBConfig;
    iterate_threshold::T = convert(T, 10),
    discount_factor::T = convert(T, 0.5),
    max_iters::Integer = 5,
    verbose::Bool = false,
    ) where T <: AbstractFloat

    λ_max = config.λ_max
    N_bands = config.N_bands - 1

    #
    Ws = Vector{Vector{T}}(undef, 1)
    Rs = Vector{Vector{T}}(undef, 1)
    Ws[begin], Rs[begin], R_norm = computewarpsamples(
        v, L, λ_max, N_bands,
    )
    W_norm = norm(Ws[begin])
    v_next = Rs[begin]
    R_prev_norm = R_norm
    #W_prev_norm = W_norm

    λs = collect( λ_max for _ = 1:1 )

    λ = λ_max
    iter = 2
    while R_norm / W_norm > iterate_threshold && 
        R_norm <= R_prev_norm && #&& W_norm <= W_prev_norm &&
        iter <= max_iters

        R_prev_norm = R_norm
        #W_prev_norm = W_norm

        λ = λ*discount_factor
        W, R, R_norm = computewarpsamples(
            v_next, L, λ, N_bands,
        )
        W_norm = norm(W)

        if R_norm <= R_prev_norm #&& W_norm >= W_prev_norm
            push!(Ws, W)
            push!(Rs, R)
            push!(λs, λ)
        end

        v_next = R

        if verbose
            @show iter, R_norm, W_norm
            @show R_norm / W_norm > iterate_threshold
            @show R_norm <= R_prev_norm
            @show R_norm, R_prev_norm, λ
        end

        iter += 1
    end
    #@show R_norm, norm(Ws[begin])

    return Ws, Rs, λs
end

# one way to aggregate the response.
function sumnormalized(xs::Vector{Vector{T}}) where T <: AbstractFloat
    return sum(
        xs[l] ./ norm(xs[l]) for l in eachindex(xs)
    )
end