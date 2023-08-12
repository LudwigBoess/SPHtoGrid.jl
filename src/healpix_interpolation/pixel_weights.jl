"""
    contributing_area(dx, hsml, ang_pix)

Calculates the area of the pixel the particle contributes to.
"""
function contributing_area(r::T, proj_hsml::T, ang_pix::T) where {T}
    max(0.0, min(ang_pix, abs(proj_hsml - (r - 0.5*ang_pix)))) / ang_pix
end


"""
    distance_to_pixel_center(r::T, pos::Vector{T}, pixel_center::Vector{T}) where T<:Real

Compute the distance between particle vector and pixel center in radians.
"""
function distance_to_pixel_center(r::T, pos::Vector{T}, pixel_center::Tuple{T,T,T}) where T<:Real
    d = 0.0
    @inbounds for i = eachindex(pos)
        d += pos[i] * pixel_center[i]
    end
    acos( min(d/r, 1.0) )
end

"""
    weight_per_index(_Δx, _pos, _hsml, _hsml_inv,
                    _pixidx,
                    distr_area, distr_weight, 
                    n_tot_pix, n_distr_pix, 
                    res, ang_pix,
                    kernel)

Helper function to compute the weights for one pixel
"""
function weight_per_index(_Δx::T, _pos::Vector{T}, _hsml::T, _hsml_inv::T,
                          _pixidx::Integer,
                          distr_area::T, distr_weight::T, 
                          n_tot_pix::Integer, n_distr_pix::Integer, 
                          res, ang_pix::T,
                          kernel::AbstractSPHKernel) where T

    # get vector to pixel center at horizon of unit sphere
    pixel_center = pix2vecRing(res, _pixidx)

    # compute distance to pixel center 
    dx = distance_to_pixel_center(_Δx, _pos, pixel_center)
    # convert in units of hsml
    u = dx * _hsml_inv

    # fraction of particle area in the pixel
    _A = contributing_area(dx, _hsml, ang_pix)
    # account for size of pixel on particle horizon 
    _A /= (ang_pix * _Δx)^2

    # total area over which particle is distributed 
    distr_area += _A

    # count up total pixels 
    n_tot_pix += 1

    # if the pixel center is within the particle
    if u <= 1
        # evaluate kernel
        _wk = 𝒲(kernel, u, _hsml_inv)

        # total distributed weight
        distr_weight += _wk * _A
        # number of pixels the particle contributes to
        n_distr_pix += 1
    else
        _wk = 0.0
    end

    return _A, _wk, 
        distr_area, distr_weight, 
        n_tot_pix, n_distr_pix
end

"""
    calculate_weights(wk::Vector{Float64}, A::Vector{Float64}, 
                    _pos::Vector{Float64}, _hsml::Float64,
                    _Δx::Float64, res::Resolution,
                    pixidx::Vector{Int64}, ang_pix::Float64,
                    kernel::AbstractSPHKernel)


"""
function calculate_weights(wk::Vector{T}, A::Vector{T}, 
                            _pos::Vector{T}, _hsml::T,
                            _Δx::T, res::Resolution,
                            pixidx::Vector{Int64}, ang_pix::T,
                            kernel::AbstractSPHKernel) where T

    # number of pixels to which the particle contributes
    Npixels = length(pixidx)

    # compute here once
    hsml_inv = 1 / _hsml

    # storage variables for count operations
    n_distr_pix  = 0
    n_tot_pix    = 0
    distr_weight = 0.0
    distr_area   = 0.0

    # loop over all relevant pixels
    @inbounds for ipixel = 1:Npixels

        A[ipixel], wk[ipixel], 
        distr_area, distr_weight, 
        n_tot_pix, n_distr_pix = weight_per_index(_Δx, _pos, _hsml, hsml_inv,
                                                  pixidx[ipixel],
                                                  distr_area, distr_weight, 
                                                  n_tot_pix, n_distr_pix, 
                                                  res, ang_pix,
                                                  kernel)

    end # loop over all relevant pixels

    # if particle contributes to pixels
    # but does not overlap with any pixel center
    if iszero(distr_weight)
        
        n_distr_pix = n_tot_pix

        # write full particle quantity into the pixel
        wk[1:Npixels] .= 1.0
        
        # the weight is normalized by the pixel area
        if !iszero(distr_area)
            weight_per_pix = n_distr_pix / distr_area
        else
            weight_per_pix = 1
        end
    else 
        weight_per_pix = n_distr_pix / distr_weight
    end

    return wk, A, n_distr_pix, weight_per_pix

end