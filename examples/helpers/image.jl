
function loadkodakimage(::Type{T}, file_path::String; discard_pixels::Integer = 0) where T <: AbstractFloat
    
    img = Images.load(file_path)
    gray_img = Images.Gray.(img)
    y_nD = convert(Array{T}, gray_img)
    
    #some kodak images have border artefacts. remove pixels.
    a = discard_pixels
    y_nD = y_nD[begin+a:end-a, begin+a:end-a]    

    return y_nD
end

function image2samples(y_nD::Matrix{T}) where T

    Nr, Nc = size(y_nD)
    x_ranges = getmatrixranges(Nr, Nc)

    y = y_nD

    return y, x_ranges
end

function getmatrixranges(Nr::Integer, Nc::Integer)

    v_range = LinRange(1, Nr, Nr)
    h_range = LinRange(1, Nc, Nc)
    

    x_ranges = Vector{LinRange{T,Int}}(undef, 2)
    x_ranges[1] = v_range
    x_ranges[2] = h_range

    return x_ranges
end
