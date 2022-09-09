mutable struct PixelAtomic
    weight::Float64 
    N::Int64
end


function calculate_weights(wk::Vector{Float64}, _pos::Vector{Float64}, _hsml::Float64,
    _Δx::Float64, res::Resolution,
    pixidx::Vector{Int64}, pix_radian::Float64,
    kernel::AbstractSPHKernel)

    Npixels = length(pixidx)

    # set kernel weight array to zero 
    wk[1:Npixels] .= 0.0

    # compute here once
    hsml_inv = 1 / _hsml

    # storage variable for atomic operations
    a = PixelAtomic(0.0, 0)

    # loop over all relevant pixels
    @inbounds for ipixel = 1:Npixels

        # get pixel center at horizon of particle
        pixel_center = _Δx .* pix2vecRing(res, pixidx[ipixel])

        # compute distance to pixel center in units of hsml
        u = get_norm(_pos .- pixel_center) * hsml_inv

        @inbounds if u <= 1

            # maximum radius of pixel at particle distance
            pix_radius = _Δx * sin(pix_radian)

            # fraction of particle in the image
            dxdy = π * (_hsml - pix_radius)^2

            # evaluate kernel
            wk[ipixel] = 𝒲(kernel, u, hsml_inv)
            wk[ipixel] *= dxdy

            # atomic operations
            a.weight += wk[ipixel]
            a.N += 1
        else
            wk[ipixel] = 0.0
        end

    end # loop over all relevant pixels

    return wk, a

end


"""
    allsky_map(pos::T, hsml::T, m::T, rho::T, bin_q::T, weights::T;
               center::Vector{<:Real}=[0.0, 0.0, 0.0],
               Nside::Integer=1024,
               kernel::AbstractSPHKernel,
               show_progress::Bool=true) where {T<:Array{<:Real}}

Calculate an allsky map from SPH particles.
"""
function allsky_map(pos::T, hsml::T, m::T, rho::T, bin_q::T, weights::T;
    center::Vector{<:Real}=[0.0, 0.0, 0.0],
    Nside::Integer=1024,
    kernel::AbstractSPHKernel,
    show_progress::Bool=true) where {T<:Array{<:Real}}

    # subtract center
    pos .-= center

    map = HealpixMap{Float64,RingOrder}(Nside)
    weight_map = HealpixMap{Float64,RingOrder}(Nside)
    res = Healpix.Resolution(Nside)

    # maximum angular distance (in radians) between a pixel center 
    # and any of its corners
    pix_radian = max_pixrad(res)

    wk = Vector{Float64}(undef, length(map))

    #Δx = @. √( pos[1,:]^2 + pos[2,:]^2 + pos[3,:]^2 )

    if show_progress
        P = Progress(length(m))
    end

    @inbounds for ipart = 1:length(m)

        # general particle quantities
        area = π * hsml[ipart]^2
        dz = m[ipart] / rho[ipart] / area

        # transform position vector to spherical coordinates 
        (theta, phi) = vec2ang(pos[:, ipart]...)

        # get distance to particle
        Δx = get_norm(pos[:, ipart])

        # get radius around particle position in radians
        radius = atan(hsml[ipart] / Δx)

        # collect indices of all healpix pixels this particles contributes to 
        pixidx = queryDiscRing(res, theta, phi, radius)

        # particle only contributes to one pixel, but does not overlap with pixel center
        # -> write the full particle into one pixel
        if iszero(length(pixidx))

            # get index of contributing pixel
            ipixel = ang2pix(map, theta, phi)

            # maximum radius of pixel at particle distance
            pix_radius = Δx * sin(pix_radian)

            # fraction of particle in the image
            dxdy = π * pix_radius^2

            area_norm = area / dxdy * weights[ipart] * dz

            pix_weight = area_norm * dxdy

            map[ipixel] += bin_q[ipart] * pix_weight
            weight_map[ipixel] += pix_weight

            continue
        end

        wk, pix = calculate_weights(wk, pos[:, ipart], hsml[ipart],
            Δx, res, pixidx, pix_radian, kernel)

        weight_per_pix = pix.N / pix.weight
        kernel_norm = area / pix.N
        area_norm = kernel_norm * weight_per_pix * weights[ipart] * dz

        # loop over all relevant pixels
        @inbounds for ipixel = 1:length(pixidx)
            pix_weight = area_norm * wk[ipixel]
            map[pixidx[ipixel]] += bin_q[ipart] * pix_weight
            weight_map[pixidx[ipixel]] += pix_weight
        end # loop over all relevant pixels


        if show_progress
            next!(P)
        end
    end # loop over particles

    # construct 2D array by deprojecting
    image, mask, maskflag = project(mollweideprojinv, map, 2Nside, Nside)
    w_image, mask, maskflag = project(mollweideprojinv, weight_map, 2Nside, Nside)

    return [image[1:end] w_image[1:end]]
end