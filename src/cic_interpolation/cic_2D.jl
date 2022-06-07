

"""
    function calculate_weights_2D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                    kernel::AbstractSPHKernel,
                                    x_pixels::Integer )

Calculates the kernel- and geometric weights of the pixels a particle contributes to.
"""
@fastmath function calculate_weights( wk::Vector{Float64}, 
                                      iMin::Integer, iMax::Integer, 
                                      jMin::Integer, jMax::Integer,
                                      x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                      kernel::AbstractSPHKernel,
                                      x_pixels::Integer )

    is_undersampled = false

    if hsml <= 1
        is_undersampled = true
    end

    distr_weight = 0.0
    n_distr_pix  = 0


    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            # index in flattened 2D array
            idx = calculate_index(i, j, x_pixels)

            dxdy = dx * dy

            if is_undersampled

                wk[idx]       = dxdy
                distr_weight += dxdy
                n_distr_pix  += 1

            else # is_undersampled

                u = get_d_hsml(x_dist, y_dist, hsml_inv)

                # check if inside the kernel
                if u <= 1
                    # evaluate kernel
                    wk[idx]       = 𝒲(kernel, u, hsml_inv)
                    wk[idx]      *= dxdy
                    distr_weight += wk[idx]
                    n_distr_pix  += 1
                else
                    wk[idx] = 0.0
                end # u < 1.0

            end # is_undersampled

        end # j
    end # i

    return wk, n_distr_pix, distr_weight
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
    area     = π * hsml^2

    rho *= T(1.0/(len2pix*len2pix*len2pix))
    dz   = m / rho / area

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
                        calc_mean=true )

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

    if show_progress
        P = Progress(N)
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Float64(Bin_Q[p])

        if iszero(bin_q) && !calc_mean
            continue
        end

        _pos, weight, hsml, hsml_inv, area, dz = get_quantities_2D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        for k = k_start:7

            # simplify position quantities for performance
            x, y, z, skip_k = get_xyz( _pos, HSML[p], k, param)

            if skip_k
                continue
            end
            
            # calculate relevant pixel range
            iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
            jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )

            wk, n_distr_pix, distr_weight = calculate_weights(wk, iMin, iMax, jMin, jMax,
                                                              x, y, hsml, hsml_inv, 
                                                              kernel,
                                                              param.Npixels[1])
           
            # skip if the particle is not contained in the image 
            # ( should only happen with periodic boundary conditions )
            if n_distr_pix == 0
                continue
            end

            weight_per_pix  = n_distr_pix / distr_weight
            kernel_norm     = area / n_distr_pix
            area_norm       = kernel_norm * weight_per_pix * weight * dz

            @inbounds for i = iMin:iMax, j = jMin:jMax
                idx = calculate_index(i, j, param.Npixels[1])

                image[idx,1], image[idx,2] = update_image( image[idx,1], image[idx,2], 
                                                           wk[idx], bin_q, area_norm)
                
            end # i, j 
        end # k

         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    return image

end # function