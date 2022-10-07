mutable struct PixelAtomic
    weight::Float64 
    N::Int64
end

using Statistics

using Base.Threads

function calculate_weights(wk::Vector{Float64}, A::Vector{Float64}, 
    _pos::Vector{Float64}, _hsml::Float64,
    _Δx::Float64, res::Resolution,
    pixidx::Vector{Int64}, pix_radius::Float64,
    kernel::AbstractSPHKernel)

    # number of pixels to which the particle contributes
    Npixels = length(pixidx)

    # set kernel weight array to zero 
    wk[1:Npixels] .= 0.0

    # set contributing area weight array to zero 
    A[1:Npixels] .= 0.0

    # compute here once
    hsml_inv = 1 / _hsml

    # storage variable for atomic operations
    #a = PixelAtomic(0.0, 0)
    N      = 0
    weight = 0.0

    # loop over all relevant pixels
    @inbounds for ipixel = 1:Npixels

        # get vector to pixel center at horizon of particle
        pixel_center = _Δx .* pix2vecRing(res, pixidx[ipixel])

        # compute distance to pixel center 
        dx = get_norm(_pos .- pixel_center)
        # convert in units of hsml
        u = dx * hsml_inv

        # fraction of particle area in the pixel
        A[ipixel] = max(0.0, min((2*pix_radius)^2, 2*pix_radius * ( dx + pix_radius - _hsml)))

        # evaluate kernel if pixel center is within particle
        if u <= 1
            # evaluate kernel
            wk[ipixel] = 𝒲(kernel, u, hsml_inv)

            # total distributed weight
            weight += wk[ipixel] * A[ipixel]
            # number of pixels the particle contributes to
            N      += 1
        end

    end # loop over all relevant pixels

    # if particle contributes to pixels
    # but does not overlap with any pixel center
    if iszero(weight)
        
        # write full particle quantity into the pixel
        wk[1:Npixels] .= 1.0
        
        # the weight is normalized by the pixel area
        area_sum = sum(A[1:Npixels])
        if !iszero(area_sum)
            weight_per_pix = 1 / area_sum
        else
            weight_per_pix = 1
        end
        
        # set N by hand
        N = Npixels
    else 
        weight_per_pix = N / weight
    end

    return wk, A, N, weight_per_pix

end


function find_in_shell(Δx, radius_limits)
    findall( @. ( radius_limits[1] <= Δx <= radius_limits[2] ) )
end


"""
    allsky_map(pos::T, hsml::T, m::T, rho::T, bin_q::T, weights::T;
               center::Vector{<:Real}=[0.0, 0.0, 0.0],
               Nside::Integer=1024,
               kernel::AbstractSPHKernel,
               show_progress::Bool=true) where {T<:Array{<:Real}}

Calculate an allsky map from SPH particles.
"""
function allsky_map(pos, hsml, m, rho, bin_q, weights;
    center::Vector{<:Real}=[0.0, 0.0, 0.0],
    Nside::Integer=1024,
    kernel::AbstractSPHKernel,
    show_progress::Bool=true,
    radius_limits::Vector{<:Real}=[0.0, Inf])

    # worker ID for output
    min_worker = minimum(workers())

    # subtract center
    pos .-= center

    allsky_map = HealpixMap{Float64,RingOrder}(Nside)
    weight_map = HealpixMap{Float64,RingOrder}(Nside)
    res = Healpix.Resolution(Nside)

    # maximum angular distance (in radians) between a pixel center 
    # and any of its corners
    pix_radian = sqrt(4π/length(allsky_map)) #max_pixrad(res)

    # storage array for kernel weights
    wk = Vector{Float64}(undef, length(allsky_map))
    # storage array for mapped area
    A  = Vector{Float64}(undef, length(allsky_map))

    # calculate radii of all particles
    _Δx = @. √( pos[1,:]^2 + pos[2,:]^2 + pos[3,:]^2 )

    # select contributing particles
    sel = find_in_shell(_Δx, radius_limits)

    if show_progress
        P = Progress(length(sel))
        println("min, max dist: $(minimum(_Δx)), $(maximum(_Δx))")
        println("radius limits: $radius_limits")
        println("$(length(sel)) / $(length(m)) in image")
    end

    gt1_pixel = 0

    @inbounds for ipart ∈ sel#[9297073]# sel

        # get distance to particle
        Δx = get_norm(pos[:, ipart])

        # if the particle is closer than its smoothing length
        # we get too much noise in the map
        if Δx < hsml[ipart]
            continue
        end

        # maximum radius of pixel at particle distance
        pix_radius   = 0.5 * Δx * tan(pix_radian)
        # area of one pixel at particle distance
        pix_area_inv = 1 / (2*pix_radius)^2

        # general particle quantities
        #area = π * hsml[ipart]^2
        dz = 2*hsml[ipart]
        area = m[ipart] / rho[ipart] / dz

        # transform position vector to spherical coordinates 
        (theta, phi) = vec2ang(pos[:, ipart]...)

        # get radius around particle position in radians
        radius = atan(hsml[ipart] / Δx)

        # collect indices of all healpix pixels this particles contributes to 
        pixidx = queryDiscRing(res, theta, phi, radius)

        # also add the index of the pixel which contains the particle center
        push!(pixidx, ang2pix(allsky_map, theta, phi))

        # paranoia check for index uniqueness
        unique!(pixidx)

        # calculate kernel weights and mapped area
        wk, A, N, weight_per_pix = calculate_weights(wk, A, pos[:, ipart], hsml[ipart],
                                            Δx, res, pixidx, pix_radius, kernel)

        # normalisation factors for pixel contribution
        kernel_norm = area * pix_area_inv / N
        area_norm = kernel_norm * weight_per_pix * weights[ipart] * dz

        # loop over all relevant pixels
        @inbounds for ipixel = 1:length(pixidx)
            pix_weight = area_norm * wk[ipixel] * A[ipixel]
            allsky_map[pixidx[ipixel]] += bin_q[ipart] * pix_weight
            weight_map[pixidx[ipixel]] += pix_weight
        end # loop over all relevant pixels

        if show_progress
            next!(P)
        end
    end # loop over particles

    if show_progress
        println()
        @info "Number of zero pixels in image: $(length(findall(iszero.(allsky_map))))"
        println()
    end

    # # construct 2D array by deprojecting
    # image, mask, maskflag = project(mollweideprojinv, allsky_map, 2Nside, Nside)
    # w_image, mask, maskflag = project(mollweideprojinv, weight_map, 2Nside, Nside)

    # return [image[1:end] w_image[1:end]]

    return allsky_map, weight_map
end