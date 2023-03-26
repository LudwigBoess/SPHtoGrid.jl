function get_weight_per_pixel( _A, _wk, distr_area, distr_weight, 
                            n_tot_pix, n_distr_pix, 
                            dxdy, u, hsml_inv,
                            kernel)

    # store are of pixel the particle contributes to
    _A          = dxdy
    distr_area +=  dxdy

    # count up total pixels 
    n_tot_pix += 1

    # if pixel center is within kernel
    if u <= 1
        # evaluate kernel
        _wk = 𝒲(kernel, u, hsml_inv)

        # count up distributed weigh
        distr_weight += _wk * dxdy

        # if pixel center is within the kernel
        # count that pixel contributes to pixel
        n_distr_pix  += 1
    else
        _wk = 0.0
    end

    return _A, _wk, 
        distr_area, distr_weight, 
        n_tot_pix, n_distr_pix
end


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
            n_tot_pix, n_distr_pix = get_weight_per_pixel(  A[idx], wk[idx], distr_area, distr_weight, 
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
                        calc_mean::Bool=true )

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

        grid_vol = 0.0

        @inbounds for i = iMin:iMax, j = jMin:jMax

            idx = calculate_index(i, j, param.Npixels[1])

            image[idx,1], image[idx,2] = update_image(image[idx,1], image[idx,2], wk[idx], A[idx], area_norm, bin_q)

            grid_vol += wk[idx] * A[idx] * dz #/ param.len2pix^3
            
        end # i, j    

        grid_mass += Rho[p] * grid_vol

        # store mass of contributing particle 
        particle_mass += M[p]

         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    if show_progress
        @info "Mass conservation:"
        @info "\tMass on grid:      $(grid_mass*1.e10) Msun"
        @info "\tMass in particles: $(particle_mass*1.e10) Msun"
        @info "\tRel. Error:        $(abs(particle_mass-grid_mass)/particle_mass)"
    end

    return image

end # function