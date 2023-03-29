"""
    contributing_pixels(pos, hsml, r)

Computes the pixels the particle contributes to.
"""
function contributing_pixels(pos::Vector{T}, hsml::T, r::T, res::Resolution, allsky_map) where T

    # transform position vector to spherical coordinates 
    (theta, phi) = vec2ang(pos...)

    # get radius around particle position in radians
    radius = atan(hsml / r)

    # collect indices of all healpix pixels this particles contributes to 
    pixidx = queryDiscRing(res, theta, phi, radius)

    # also add the index of the pixel which contains the particle center
    push!(pixidx, ang2pix(allsky_map, theta, phi))

    # paranoia check for index uniqueness
    unique!(pixidx)

    return pixidx
end
