"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function convertcompactdomain(x::Vector{T}, a, b, c, d)::Vector{T} where T <: Real

    return collect( convertcompactdomain(x[i], a[i], b[i], c[i], d[i]) for i = 1:length(x) )
end

