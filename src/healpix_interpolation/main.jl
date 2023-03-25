using Statistics

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
    particle_area_and_depth(d, hsml, pix_radian)

Computes the particle area and depth.
Caution: This has to be represented as a cylinder instead of a sphere, which introduces an error by design.
Also computes the pixel radius at the particle horizon.
"""
function particle_area_and_depth(hsml::T, m::T, rho::T) where {T}

    # general particle quantities
    dz = 2hsml
    area = (m / rho) / dz

    return area, dz
end

"""
    particle_area_and_depth(d, hsml, pix_radian)

Computes the particle area and depth.
Caution: This has to be represented as a cylinder instead of a sphere, which introduces an error by design.
Also computes the pixel radius at the particle horizon.
"""
function particle_area_and_depth(hsml::T) where {T}

    # general particle quantities
    dz = 2hsml
    area = π * hsml^2

    return area, dz
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
function healpix_map(Pos, Hsml, M, Rho, Bin_q, Weights;
    center::Vector{<:Real}=[0.0, 0.0, 0.0],
    Nside::Integer=1024,
    kernel::AbstractSPHKernel,
    show_progress::Bool=true,
    radius_limits::Vector{<:Real}=[0.0, Inf],
    calc_mean::Bool=true)

    # worker ID for output
    min_worker = minimum(workers())

    # Npart before filtering
    Npart_in = length(Hsml)

    if show_progress && myid() == min_worker
        println("filtering particles")
    end

    # filter and sort the particles, return arrays with only relevant particles
    pos, hsml, m, rho, bin_q, weights = filter_sort_particles(Pos, Hsml,  M, Rho, Bin_q, Weights, center, radius_limits)

    # Npart after filtering 
    Npart = length(hsml)

    if show_progress && myid() == min_worker
        println("$Npart / $Npart_in in image")
    end

    # allocate arrays
    allsky_map = HealpixMap{Float64,RingOrder}(Nside)
    weight_map = HealpixMap{Float64,RingOrder}(Nside)
    res = Healpix.Resolution(Nside)

    # maximum angular distance (in radians) between a pixel center 
    # and any of its corners
    #tan_pix_radian = tan(0.5√(4π/length(allsky_map)))
    tan_pix_radian = tan(max_pixrad(res))

    # storage array for kernel weights
    wk = Vector{Float64}(undef, length(allsky_map))
    # storage array for mapped area
    A = Vector{Float64}(undef, length(allsky_map))

    # storage for grid and particle masses 
    grid_mass = 0.0
    part_mass = 0.0

    # define progressmeter
    P = Progress(Npart)

    @inbounds for ipart ∈ 1:Npart

        if !calc_mean
            if iszero(bin_q[ipart])
                update_progress!(P, min_worker, show_progress)
                continue
            end
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
        area, dz = particle_area_and_depth(hsml[ipart], m[ipart], rho[ipart])
        #area, dz = particle_area_and_depth(hsml[ipart])

        # calculate kernel weights and mapped area
        wk, A, N, weight_per_pix = calculate_weights(wk, A, pos[:, ipart], hsml[ipart],
            Δx, res, pixidx, pix_radius, kernel)

        # update the actual images
        update_image!(allsky_map, weight_map, wk,
            pixidx,
            area, dz,
            A, N, weight_per_pix,
            weights[ipart], bin_q[ipart])

        # store mass on grid and in particles 
        grid_mass += rho[ipart] * sum(wk[1:length(pixidx)]) * sum(A[1:length(pixidx)]) * dz
        part_mass += m[ipart]

        # update the progress meter
        update_progress!(P, min_worker, show_progress)
    end # loop over particles

    if show_progress && myid() == min_worker
        println()
        @info "Number of zero pixels in image: $(length(findall(iszero.(allsky_map))))"
        println()
        @info "Mass conservation:"
        @info "\tMass on grid:      $(grid_mass*1.e10) Msun"
        @info "\tMass in particles: $(part_mass*1.e10) Msun"
        @info "\tRel. Error:        $(abs(part_mass-grid_mass)/part_mass)"
    end

    return allsky_map, weight_map
end