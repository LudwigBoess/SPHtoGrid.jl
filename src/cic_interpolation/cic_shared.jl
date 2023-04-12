"""
    function pix_index_min_max(x::T, hsml::T, n_pixels::Int64) where T

Calculates the minimum and maximum pixel to which a particle contributes.
"""
function pix_index_min_max(x::T, hsml::T, n_pixels::Int64) where T

    iMin = max(floor(Integer, x - hsml), 0)
    iMax = min(floor(Integer, x + hsml), n_pixels-1)

    return iMin, iMax
end


"""
    function get_dxyz(x::Real, hsml::Real, i::Integer)

Calculates the extent of the current particle size in units of pixels.
"""
function get_dxyz(x::Real, hsml::Real, i::Integer)
    min(x + hsml, i + 1) - max(x - hsml, i)
end


"""
    function get_x_dx(x, hsml, i)


"""
@inline function get_x_dx(x, hsml, i)

    dx = get_dxyz(x, hsml, i)
    x_dist = x - i - 0.5

    return x_dist, dx
end


"""
    function get_xyz( pos, hsml, k::Integer,
                      par::mappingParameters)

Calculates `x, y, z` position in units of pixels and performs periodic mapping, if required.
"""
function get_xyz( pos::Array{T},
                  par::mappingParameters) where T

    x, y, z = pos
    # convert from code units to pixels
    x *= par.len2pix
    y *= par.len2pix
    z *= par.len2pix

    x += 0.5 * par.Npixels[1]
    y += 0.5 * par.Npixels[2]
    z += 0.5 * par.Npixels[3]
    
    return x, y, z

end



"""
    function update_image( image::Real, w_image::Real, 
                           wk::Real, bin_q::Real, 
                           geometry_norm::Real )

Applies the different contributions to the image and the weight image.
"""
function update_image( image::Float64, w_image::Float64, 
                        wk::Float64, A::Float64, geometry_norm::Float64,
                        bin_q::Float64 )

    if wk > 0
        pix_weight  = geometry_norm * wk * A
        image      += bin_q * pix_weight
        w_image    += pix_weight
    end

    return image, w_image
end




function free_memory(x, hsml, m, rho, bin_q, weights)
    x = hsml = m = rho = bin_q = weights = nothing
    GC.gc()
end