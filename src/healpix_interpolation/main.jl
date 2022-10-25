mutable struct PixelAtomic
    @atomic distr_weight::Float64 
    @atomic distr_area::Float64 
    @atomic n_distr_pix::Int64
    @atomic n_tot_pix::Int64
end

using Statistics

using Base.Threads

function get_area(dx, hsml, pix_radius)
    max(0.0, min((2pix_radius)^2, 2pix_radius * ( pix_radius + hsml - dx)))
    #1.0
end

function calculations_per_index(_Δx, _pos, _hsml, _hsml_inv,
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
    _A = get_area(dx, _hsml, pix_radius)
    # total area over which particle is distributed 
    distr_area += _A

    # count up total pixels 
    n_tot_pix += 1

    # number of pixels the particle contributes to
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

function calculate_weights(wk::Vector{Float64}, A::Vector{Float64}, 
    _pos::Vector{Float64}, _hsml::Float64,
    _Δx::Float64, res::Resolution,
    pixidx::Vector{Int64}, pix_radius::Float64,
    kernel::AbstractSPHKernel)

    # number of pixels to which the particle contributes
    Npixels = length(pixidx)

    # compute here once
    hsml_inv = 1 / _hsml

    # storage variable for atomic operations
    # a = PixelAtomic(0.0, 0.0, 0, 0)
    n_distr_pix  = 0
    n_tot_pix    = 0
    distr_weight = 0.0
    distr_area   = 0.0

    # loop over all relevant pixels
    @inbounds for ipixel = 1:Npixels

        A[ipixel], wk[ipixel], 
        distr_area, distr_weight, 
        n_tot_pix, n_distr_pix = calculations_per_index(_Δx, _pos, _hsml, hsml_inv,
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


function find_in_shell(Δx, radius_limits)
    @. ( radius_limits[1] <= Δx <= radius_limits[2] )
end


function get_pixel_quantities(r, _pos, _hsml, pix_radian, allsky_map, res)

    # maximum radius of pixel at particle distance
    pix_radius   = r * tan(pix_radian)
    # area of one pixel at particle distance
    pix_area_inv = 1 / (2*pix_radius)^2

    # general particle quantities
    area = π * _hsml^2
    dz = 2*_hsml
    #area = m[ipart] / rho[ipart] / dz

    # transform position vector to spherical coordinates 
    (theta, phi) = vec2ang(_pos...)

    # get radius around particle position in radians
    radius = atan(_hsml / r)

    # collect indices of all healpix pixels this particles contributes to 
    pixidx = queryDiscRing(res, theta, phi, radius)

    # also add the index of the pixel which contains the particle center
    push!(pixidx, ang2pix(allsky_map, theta, phi))

    # paranoia check for index uniqueness
    unique!(pixidx)

    return pixidx, pix_radius, area, dz
end

"""
    healpix_map(pos::T, hsml::T, m::T, rho::T, bin_q::T, weights::T;
               center::Vector{<:Real}=[0.0, 0.0, 0.0],
               Nside::Integer=1024,
               kernel::AbstractSPHKernel,
               show_progress::Bool=true) where {T<:Array{<:Real}}

Calculate an allsky map from SPH particles.
"""
function healpix_map(pos, hsml, m, rho, bin_q, weights;
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

    #sel = sel[sortperm(_Δx[sel])]

    if show_progress && myid() == min_worker
        P = Progress(length(sel))
        println("min, max dist: $(minimum(_Δx)), $(maximum(_Δx))")
        println("radius limits: $radius_limits")
        println("$(length(findall(sel))) / $(length(_Δx)) in image")
    end

    @inbounds for ipart ∈ 1:length(sel)

        if !sel[ipart]
            continue 
        end

        # get distance to particle
        Δx = get_norm(pos[:, ipart])

        # if the particle is closer than its smoothing length
        # we get too much noise in the map
        if Δx < hsml[ipart]
            continue
        end

        pixidx, pix_radius, area, dz = get_pixel_quantities(Δx, pos[:, ipart], hsml[ipart], pix_radian, allsky_map, res)
        
        # calculate kernel weights and mapped area
        wk, A, N, weight_per_pix = calculate_weights(wk, A, pos[:, ipart], hsml[ipart],
                                                        Δx, res, pixidx, pix_radius, kernel)

        # normalisation factors for pixel contribution
        kernel_norm = area / N #* pix_area_inv 
        area_norm = kernel_norm * weight_per_pix * weights[ipart] * dz

        # loop over all relevant pixels
        @inbounds for ipixel = 1:length(pixidx)
            pix_weight = area_norm * wk[ipixel] * A[ipixel]
            allsky_map[pixidx[ipixel]] += bin_q[ipart] * pix_weight
            weight_map[pixidx[ipixel]] += pix_weight
        end # loop over all relevant pixels

        if show_progress && myid() == min_worker
            next!(P)
        end
    end # loop over particles

    if show_progress && myid() == min_worker
        println()
        @info "Number of zero pixels in image: $(length(findall(iszero.(allsky_map))))"
        println()
    end

    return allsky_map, weight_map
end