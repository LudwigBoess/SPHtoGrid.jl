"""
    function calculate_weights_2D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                    kernel::AbstractSPHKernel,
                                    y_pixels::Integer )

Calculates the kernel- and geometric weights of the pixels a particle contributes to.
`y_pixels` is the number of pixels along the (inner-loop) y axis, i.e. the row
stride of the flattened image.
"""
@fastmath function calculate_weights( wk::Vector{Float64}, A::Vector{Float64},
                                      iMin::Integer, iMax::Integer,
                                      jMin::Integer, jMax::Integer,
                                      x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                      kernel::AbstractSPHKernel,
                                      y_pixels::Integer )

    # storage variables for count operations
    n_distr_pix  = 0
    n_tot_pix    = 0
    distr_weight = 0.0
    distr_area   = 0.0


    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            # projected distance to pixel center in units of hsml
            u = get_d_hsml(x_dist, y_dist, hsml_inv)

            # pixel area 
            dxdy = dx * dy

            # index in flattened 2D array
            idx = calculate_index(i, j, y_pixels)

            A[idx], wk[idx], 
            distr_area, distr_weight, 
            n_tot_pix, n_distr_pix = get_weight_per_pixel( distr_area, distr_weight, 
                                                        n_tot_pix, n_distr_pix, 
                                                        dxdy, u, hsml_inv, kernel)

        end # j
    end # i

    # if particle contributes to pixels
    # but does not overlap with any pixel center
    if iszero(distr_weight)
        
        n_distr_pix = n_tot_pix

        # write full particle quantity into the pixel
        @inbounds for i = iMin:iMax, j = jMin:jMax
            idx = calculate_index(i, j, y_pixels)
            wk[idx] = 1.0
        end
        
        # the weight is normalized by the pixel area
        if !iszero(distr_area)
            weight_per_pix = n_distr_pix / distr_area
        else
            # particle touches no in-grid pixel (clamped empty range);
            # use 1.0 (not Int 1) for type stability. The caller guards
            # against n_distr_pix == 0 so this particle is skipped and
            # does not corrupt the image with Inf/NaN.
            weight_per_pix = 1.0
        end
    else
        weight_per_pix = n_distr_pix / distr_weight
    end

    return wk, A, n_distr_pix, weight_per_pix
end

"""
    get_quantities_2D( pos, weight, hsml, 
                       rho, m, len2pix::T) where T

Helper function to convert quantities to pixel units and the correct data type.
"""
function get_quantities_2D( pos, weight, hsml, 
                            rho, m, len2pix::T) where T
        
    hsml    *= T(len2pix)
    hsml_inv = T(1/hsml)
    area     = (2hsml)^2 # Effective area of squared particle [pix^2]

    rho *= T(1/(len2pix*len2pix*len2pix)) # [10^10 Msun * pix^-3]
    dz   = m / rho / area # [pix]

    return T.(pos), T(weight), hsml, hsml_inv, area, dz
end


"""
    _cic_map_particle!( image, wk, A, p, Pos, HSML, M, Rho, Bin_Q, Weights, RM,
                        touched_pixel, param, kernel, N_images, calc_mean, stokes )

Deposits a single particle `p` onto `image`, using the reusable scratch
buffers `wk`/`A`. Factored out of [`cic_mapping_2D`](@ref) so the particle
loop can run either serially or split across threads, each thread owning its
own `image`/`wk`/`A`.

Returns `true` when the caller's progress meter should advance (i.e. the
particle was considered for mapping) and `false` for a particle skipped early
because its quantity is zero (matching the original progress semantics).
"""
function _cic_map_particle!(image::Array{Float64}, wk::Vector{Float64}, A::Vector{Float64},
                            p::Integer,
                            Pos, HSML, M, Rho, Bin_Q, Weights, RM,
                            touched_pixel,
                            param::mappingParameters, kernel::AbstractSPHKernel,
                            N_images::Integer, calc_mean::Bool, stokes::Bool)

    # assign bin quantity and convert to Float64
    if N_images == 1
        bin_q = Float64(Bin_Q[p])
    else
        bin_q = Float64.(Bin_Q[:,p])

        if bin_q == zeros(length(bin_q))
            bin_q = 0.0
        end
    end

    # if the quantity is zero we can skip the particle
    if iszero(bin_q) && !calc_mean
        return false
    end

    _pos, los_weight, hsml, hsml_inv, area, dz = get_quantities_2D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

    # simplify position quantities for performance
    x, y, z = get_xyz( _pos, param)

    # calculate relevant pixel range
    iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
    jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )

    # calculate all relevant quantities
    wk, A, n_distr_pix, weight_per_pix  = calculate_weights(wk, A,
                                                  iMin, iMax, jMin, jMax,
                                                  x, y, hsml, hsml_inv,
                                                  kernel,
                                                  param.Npixels[2])

    # particle touches no in-grid pixel (e.g. clamped out at the edge or
    # off-grid): skip it instead of dividing area / 0 = Inf and writing
    # NaNs into the image. The off-grid fraction is lost by design for a
    # sub-region map.
    if iszero(n_distr_pix)
        return true
    end

    # normalisation factors for pixel contribution
    kernel_norm = area / n_distr_pix
    area_norm = kernel_norm * weight_per_pix * los_weight * dz

    # loop over all contributing pixels
    @inbounds for i = iMin:iMax, j = jMin:jMax

        # get the current index in the image array
        idx = calculate_index(i, j, param.Npixels[2])

        # compute pixel weight
        pix_weight = wk[idx] * A[idx] * area_norm

        # if we map faraday rotation we need to rotate the previous emission
        if !isnothing(RM)

            # only rotate pixel if it has been processed
            if touched_pixel[idx]
                # apply faraday rotation to pixel
                faraday_rotate_pixel!(image, idx, RM[p], pix_weight, stokes)
            end
        end

        if !iszero(pix_weight)
            update_image!(image, idx, pix_weight, bin_q)

            # pixel has been processed now
            if !isnothing(RM)
                touched_pixel[idx] = true
            end
        end

    end # i, j

    return true
end

"""
   cic_mapping_2D( Pos, HSML, M, Rho, Bin_Q, Weights;
                   param::mappingParameters,
                   kernel::AbstractSPHKernel,
                   show_progress::Bool=false,
                   calc_mean::Bool=true,
                   threaded::Bool=false )

Underlying function to map SPH data to a 2D grid.

Set `threaded=true` to run the particle loop across `Threads.nthreads()`
threads. Particles are split into chunks (via `domain_decomposition`), each
chunk deposits into its own private `image`, and the partials are summed at
the end — the same race-free reduction the distributed (`@spawnat`) path
uses, without inter-process copies. The Faraday-rotation path (`RM !==
nothing`) is order-dependent and always runs serially.

When `show_progress=true`, a grid-vs-particle mass-conservation check is
logged at the end (see [`mass_conservation_report`](@ref)): the fraction of
`Σ Weights·M/Rho` captured on the grid, a *coverage* diagnostic (off-grid /
too-small particles), not an interpolation-error metric.
"""
function cic_mapping_2D( Pos, HSML,
                        M, Rho,
                        Bin_Q, Weights,
                        RM=nothing;
                        param::mappingParameters,
                        kernel::AbstractSPHKernel,
                        show_progress::Bool=false,
                        calc_mean::Bool=true,
                        stokes::Bool=false,
                        threaded::Bool=false )

    N = size(M,1)  # number of particles

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2]

    # check if Bin_Q is only one quantity or an array
    if ndims(Bin_Q) == 1
        N_images = 1
    else
        # if Bin_Q is an array we need to allocate more images
        N_images = size(Bin_Q, 1)
    end

    # NOTE on periodic boundaries (Correctness issue C3 / mass-conservation B3):
    # The CIC path does NOT deposit periodic ghost images of a particle. With
    # `boxsize` set, `center_particles` wraps particle *centers* into the box
    # (minimum-image, full-boxsize shift) and the hsml-aware particle filter
    # keeps edge-straddling particles. The per-particle renormalisation
    # (1/distr_weight over the covered pixels) then deposits the *full*
    # `m·w/ρ` of every kept particle, so total mass is conserved exactly
    # (verified to machine precision by the mass-conservation tests, incl. the
    # periodic-straddler case). The only approximation is *placement*: the
    # kernel tail of a particle within `hsml` of the box edge is piled onto
    # the near edge instead of being wrapped to the opposite edge. Use the
    # SPLASH mapping path if exact periodic-image placement is required.

    # Threaded shared-memory path: split the particles into nthreads() chunks,
    # give each chunk its own image + scratch buffers, then reduce by summing.
    # This mirrors the distributed (@spawnat) path but avoids inter-process
    # copies. The Faraday-rotation (RM) path is order-dependent — each particle
    # reads and rewrites the accumulated pixel state — so it cannot be reduced
    # and always falls through to the serial loop below.
    if threaded && isnothing(RM) && nthreads() > 1

        nt    = nthreads()
        batch = domain_decomposition(N, nt)

        # per-chunk (== per-thread) accumulators and scratch buffers
        images = [zeros(Float64, N_distr, N_images+1) for _ in 1:nt]
        wks    = [zeros(Float64, N_distr) for _ in 1:nt]
        As     = [Vector{Float64}(undef, N_distr) for _ in 1:nt]

        # index the buffers by the chunk id `t`, never by threadid(), so the
        # result is correct even if a task migrates between threads.
        @threads for t = 1:nt
            img  = images[t]
            wk_t = wks[t]
            A_t  = As[t]
            @inbounds for p in batch[t]
                _cic_map_particle!(img, wk_t, A_t, p,
                                   Pos, HSML, M, Rho, Bin_Q, Weights, nothing,
                                   nothing, param, kernel, N_images, calc_mean, stokes)
            end
        end

        image = sum(images)

    else

        # ---- serial path (also handles the RM / Faraday-rotation case) ----

        # allocate images and weight_image
        image = zeros(Float64, N_distr, N_images+1)

        touched_pixel = isnothing(RM) ? nothing : falses(N_distr)

        # allocate arrays for weights
        wk = zeros(Float64, N_distr)
        # storage array for mapped area
        A  = Vector{Float64}(undef, N_distr)

        if show_progress
            P = Progress(N)
        end

        # loop over all particles
        @inbounds for p = 1:N

            mapped = _cic_map_particle!(image, wk, A, p,
                                        Pos, HSML, M, Rho, Bin_Q, Weights, RM,
                                        touched_pixel, param, kernel,
                                        N_images, calc_mean, stokes)

            # update for ProgressMeter
            if show_progress && mapped
                next!(P)
            end
        end # p
    end

    # report grid-vs-particle mass conservation (coverage diagnostic)
    if show_progress
        mass_conservation_report(image, M, Rho, Weights, param.len2pix; ndim=2)
    end

    return image

end # function
