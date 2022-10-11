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

    area_sum     = 0.0
    distr_weight = 0.0
    n_distr_pix  = 0
    n_tot_pix    = 0

    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            # index in flattened 2D array
            idx = calculate_index(i, j, x_pixels)

            dxdy = dx * dy

            # store are of pixel the particle contributes to
            A[idx]    = dxdy
            area_sum +=  dxdy

            u = get_d_hsml(x_dist, y_dist, hsml_inv)

            # evaluate kernel
            wk[idx]       = 𝒲(kernel, u, hsml_inv)

            # count up distributed weigh
            distr_weight += wk[idx] * A[idx]

            # count up total pixels 
            n_tot_pix += 1

            # if pixel center is within the kernel
            # count that pixel contributes to pixel
            if !iszero(wk[idx])
                n_distr_pix  += 1
            end

        end # j
    end # i

    # check if particle intersects with any pixel center
    if !iszero(distr_weight)
        weight_per_pixel = n_distr_pix / distr_weight
    else # particle contributes to pixel but does not overlap with any pixel center
        n_distr_pix = n_tot_pix

        # check if particle contributes to more than one pixel 
        if !iszero(area_sum)
            weight_per_pixel = n_distr_pix / area_sum
        else 
            # particle only contributes to one pixel
            weight_per_pixel = 1
        end
    end

    return wk, A, n_distr_pix, weight_per_pixel
end

"""
    get_quantities_2D( pos, weight, hsml, 
                       rho, m, len2pix::T) where T

Helper function to convert quantities to pixel units and the correct data type.
"""
function get_quantities_2D( pos, weight, hsml, 
                            rho, m, len2pix::T) where T
        
    hsml    *= T(len2pix)
    hsml_inv = T(1.0)/hsml
    area     = (2hsml)^2 # Effective area of squared particle [pix^2]

    rho *= T(1.0/(len2pix*len2pix*len2pix)) # [10^10 Msun * pix^-3]
    dz   = m / rho / area # [pix]

    return T.(pos), T(weight), hsml, hsml_inv, area, dz
end


"""
   cic_mapping_2D( Pos, HSML, 
                  M, Rho, 
                  Bin_Q, Weights;
                  param::mappingParameters, kernel::AbstractSPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 2D grid.
"""
function cic_mapping_2D( Pos, HSML, 
                        M, Rho, 
                        Bin_Q, Weights;
                        param::mappingParameters, kernel::AbstractSPHKernel,
                        show_progress::Bool=false,
                        calc_mean=false )

    N = size(M,1)  # number of particles
    
    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2]

    image = zeros(Float64, N_distr, 2)

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

    grid_mass = 0.0
    particle_mass = 0.0

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Float64(Bin_Q[p])

        if iszero(bin_q) # && !calc_mean
            continue
        end

        _pos, los_weight, hsml, hsml_inv, particle_area, dz = get_quantities_2D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        for k = k_start:7

            # simplify position quantities for performance
            x, y, z, skip_k = get_xyz( _pos, HSML[p], k, param)

            if skip_k
                continue
            end
            
            # calculate relevant pixel range
            iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
            jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )

            # calculate all relevant quantities
            wk, A, n_distr_pix, weight_per_pixel = calculate_weights(wk, A, 
                                                              iMin, iMax, jMin, jMax,
                                                              x, y, hsml, hsml_inv, 
                                                              kernel,
                                                              param.Npixels[1])
           

            part_area_per_pix  = particle_area / n_distr_pix
            kernel_norm = part_area_per_pix * weight_per_pixel * los_weight * dz

            @inbounds for i = iMin:iMax, j = jMin:jMax

                idx = calculate_index(i, j, param.Npixels[1])

                pixel_contribution = A[idx] * wk[idx] * kernel_norm

                # actual image
                image[idx,1] += bin_q * pixel_contribution
                # weight image
                image[idx,2] += pixel_contribution

                grid_mass += Rho[p] * wk[idx] * A[idx] * dz / param.len2pix^3
                
            end # i, j
            
        end # k

        # store mass of contributing particle 
        particle_mass += M[p]

         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    if show_progress
        @info "Mass conservation:"
        @info "\tMass on grid:      $(grid_mass)"
        @info "\tMass in particles: $(particle_mass)"
        @info "\tRel. Error:        $(abs(particle_mass-grid_mass)/particle_mass)"
    end

    return image

end # function