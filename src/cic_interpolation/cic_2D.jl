"""
    function calculate_weights_2D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                    kernel::AbstractSPHKernel,
                                    x_pixels::Integer )

Calculates the kernel- and geometric weights of the pixels a particle contributes to.
"""
@fastmath function calculate_weights( wk::Vector{Float64}, A::Vector{Float64},
                                      iMin::Integer, iMax::Integer, 
                                      jMin::Integer, jMax::Integer,
                                      x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                      kernel::AbstractSPHKernel,
                                      x_pixels::Integer )

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
            idx = calculate_index(i, j, x_pixels)

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
        wk[1:n_tot_pix] .= 1.0
        
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
   cic_mapping_2D( Pos, HSML, M, Rho, Bin_Q, Weights;
                   param::mappingParameters, 
                   kernel::AbstractSPHKernel,
                   show_progress::Bool=false,
                   calc_mean::Bool=true )

Underlying function to map SPH data to a 2D grid.
"""
function cic_mapping_2D( Pos, HSML, 
                        M, Rho, 
                        Bin_Q, Weights, 
                        RM=nothing;
                        param::mappingParameters, 
                        kernel::AbstractSPHKernel,
                        show_progress::Bool=false,
                        calc_mean::Bool=true,
                        stokes::Bool=false )

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
           
    # allocate images and weight_image
    image = zeros(Float64, N_distr, N_images+1)

    if !isnothing(RM)
        touched_pixel = falses(N_distr)
    end

    if param.periodic
        k_start = 0
    else
        k_start = 7
    end

    # allocate arrays for weights
    wk = zeros(Float64, N_distr)
    # storage array for mapped area
    A  = Vector{Float64}(undef, N_distr)

    if show_progress
        P = Progress(N)
    end

    # grid_mass = 0.0
    # particle_mass = 0.0

    # loop over all particles
    @inbounds for p = 1:N

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
            continue
        end

        _pos, los_weight, hsml, hsml_inv, area, dz = get_quantities_2D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        # simplify position quantities for performance
        x, y, z = get_xyz( _pos, param)
        
        # calculate relevant pixel range
        iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
        jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )

        # calculate all relevant quantities
        wk, A, N, weight_per_pix  = calculate_weights(wk, A, 
                                                      iMin, iMax, jMin, jMax,
                                                      x, y, hsml, hsml_inv, 
                                                      kernel,
                                                      param.Npixels[1])

        # normalisation factors for pixel contribution
        kernel_norm = area / N
        area_norm = kernel_norm * weight_per_pix * los_weight * dz

        #grid_vol = 0.0

        # loop over all contributing pixels
        @inbounds for i = iMin:iMax, j = jMin:jMax

            # get the current index in the image array
            idx = calculate_index(i, j, param.Npixels[1])

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

            #grid_vol += wk[idx] * A[idx] * dz #pix_weight
            
        end # i, j    

        #grid_mass += Rho[p] * grid_vol / param.len2pix^3

        # store mass of contributing particle 
        # particle_mass += M[p]

         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    # if show_progress
    #     @info "Mass conservation:"
    #     @info "\tMass on grid:      $(grid_mass*1.e10) Msun"
    #     @info "\tMass in particles: $(particle_mass*1.e10) Msun"
    #     @info "\tRel. Error:        $(abs(particle_mass-grid_mass)/particle_mass)"
    # end

    return image

end # function