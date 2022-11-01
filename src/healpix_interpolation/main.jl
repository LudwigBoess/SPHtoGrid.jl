using Statistics

"""
    contributing_area(dx, hsml, pix_radius)

Calculates the area of the pixel the particle contributes to.
"""
function contributing_area(dx, hsml, pix_radius)
    max(0.0, min((2pix_radius)^2, 2pix_radius * ( pix_radius + hsml - dx)))
end

"""
    weight_per_index(_Δx, _pos, _hsml, _hsml_inv,
                    _pixidx,
                    distr_area, distr_weight, 
                    n_tot_pix, n_distr_pix, 
                    res, pix_radius,
                    kernel)

Helper function to compute the weights for one pixel
"""
function weight_per_index(_Δx, _pos, _hsml, _hsml_inv,
                          _pixidx,
                          distr_area, distr_weight, 
                          n_tot_pix, n_distr_pix, 
                          res, pix_radius,
                          kernel)

    # get vector to pixel center at horizon of particle
    pixel_center = _Δx .* pix2vecRing(res, _pixidx)

    # compute distance to pixel center 
    dx = get_norm(_pos .- pixel_center)
    # convert in units of hsml
    u = dx * _hsml_inv

    # fraction of particle area in the pixel
    _A = contributing_area(dx, _hsml, pix_radius)
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
                    pixidx::Vector{Int64}, pix_radius::Float64,
                    kernel::AbstractSPHKernel)


"""
function calculate_weights(wk::Vector{Float64}, A::Vector{Float64}, 
                            _pos::Vector{Float64}, _hsml::Float64,
                            _Δx::Float64, res::Resolution,
                            pixidx::Vector{Int64}, pix_radius::Float64,
                            kernel::AbstractSPHKernel)

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
                                                  res, pix_radius,
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


"""
    find_in_shell(Δx, radius_limits)

Checks if a particle is contained in a shell around 
"""
function find_in_shell(Δx, radius_limits)
    @. ( radius_limits[1] <= Δx <= radius_limits[2] )
end


"""
    particle_area_and_depth(d, hsml, pix_radian)

Computes the particle area and depth.
Caution: This has to be represented as a cylinder instead of a sphere, which introduces an error by design.
Also computes the pixel radius at the particle horizon.
"""
function particle_area_and_depth(hsml)

    # general particle quantities
    area = π * hsml^2
    dz = 2*hsml

    return area, dz
end


"""
    contributing_pixels(pos, hsml, r)

Computes the pixels the particle contributes to.
"""
function contributing_pixels(pos, hsml, r, res, allsky_map)

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


"""
    update_progress!(P, min_worker, show_progress)

Updates the progress bar if conditions are met.
"""
function update_progress!(P, min_worker, show_progress)
    if show_progress && myid() == min_worker
        next!(P)
    end
end


"""
    update_image!(allsky_map, weight_map, wk,
                pixidx,
                area, dz,
                A, N, weight_per_pix,
                quantity_weight, bin_quantity)

Updates the images.
"""
function update_image!(allsky_map, weight_map, wk,
                       pixidx,
                       area, dz,
                       A, N, weight_per_pix,
                       quantity_weight, bin_quantity)

    # normalisation factors for pixel contribution
    kernel_norm = area / N
    area_norm = kernel_norm * weight_per_pix * quantity_weight * dz

    # loop over all relevant pixels
    @inbounds for ipixel = 1:length(pixidx)
        # combine the different weighting contributions
        pix_weight = area_norm * wk[ipixel] * A[ipixel]

        # update image and weight image
        allsky_map[pixidx[ipixel]] += bin_quantity * pix_weight
        weight_map[pixidx[ipixel]] += pix_weight
    end # loop over all relevant pixels

end

"""
    healpix_map(pos, hsml, m, rho, bin_q, weights;
                center::Vector{<:Real}=[0.0, 0.0, 0.0],
                Nside::Integer=1024,
                kernel::AbstractSPHKernel,
                show_progress::Bool=true,
                radius_limits::Vector{<:Real}=[0.0, Inf])

Calculate an allsky map from SPH particles.
"""
function healpix_map(pos, hsml, bin_q, weights;
                    center::Vector{<:Real}=[0.0, 0.0, 0.0],
                    Nside::Integer=1024,
                    kernel::AbstractSPHKernel,
                    show_progress::Bool=true,
                    radius_limits::Vector{<:Real}=[0.0, Inf])

    # worker ID for output
    min_worker = minimum(workers())

    # subtract center
    pos .-= center

    # allocate arrays
    allsky_map = HealpixMap{Float64,RingOrder}(Nside)
    weight_map = HealpixMap{Float64,RingOrder}(Nside)
    res = Healpix.Resolution(Nside)

    # maximum angular distance (in radians) between a pixel center 
    # and any of its corners
    tan_pix_radian = tan(0.5√(4π/length(allsky_map)))

    # storage array for kernel weights
    wk = Vector{Float64}(undef, length(allsky_map))
    # storage array for mapped area
    A  = Vector{Float64}(undef, length(allsky_map))

    # calculate radii of all particles
    _Δx = @. √( pos[1,:]^2 + pos[2,:]^2 + pos[3,:]^2 )

    # select contributing particles
    sel = find_in_shell(_Δx, radius_limits)

    # define progressmeter
    P = Progress(length(sel))
    
    if show_progress && myid() == min_worker
        println("min, max dist: $(minimum(_Δx)), $(maximum(_Δx))")
        println("radius limits: $radius_limits")
        println("$(length(findall(sel))) / $(length(_Δx)) in image")
    end

    @inbounds for ipart ∈ 1:length(sel)

        # if the particle is not contained in the desired shell -> skip it
        # this seems to be faster than taking the new index list
        if !sel[ipart]
            update_progress!(P, min_worker, show_progress)
            continue 
        end

        # get distance to particle
        Δx = get_norm(pos[:, ipart])

        # if the particle is closer than its smoothing length
        # we get too much noise in the map
        if Δx < hsml[ipart]
            update_progress!(P, min_worker, show_progress)
            continue
        end

        # pixel radius at particle horizon
        pix_radius = Δx * tan_pix_radian

        # find pixels to which particle contributes
        pixidx = contributing_pixels(pos[:, ipart], hsml[ipart], Δx, res, allsky_map)

        # area of particle and length along line of sight
        area, dz = particle_area_and_depth(hsml[ipart])
        
        # calculate kernel weights and mapped area
        wk, A, N, weight_per_pix = calculate_weights(wk, A, pos[:, ipart], hsml[ipart],
                                                     Δx, res, pixidx, pix_radius, kernel)

        # update the actual images
        update_image!(allsky_map, weight_map, wk,
                      pixidx,
                      area, dz,
                      A, N, weight_per_pix,
                      weights[ipart], bin_q[ipart])

        # update the progress meter
        update_progress!(P, min_worker, show_progress)
    end # loop over particles

    if show_progress && myid() == min_worker
        println()
        @info "Number of zero pixels in image: $(length(findall(iszero.(allsky_map))))"
        println()
    end

    return allsky_map, weight_map
end