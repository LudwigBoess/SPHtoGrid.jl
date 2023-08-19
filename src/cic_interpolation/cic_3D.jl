"""
    function calculate_weights_3D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    kMin::Integer, kMax::Integer,
                                    x::Real, y::Real, z::Real, 
                                    hsml::Real, hsml_inv::Real,
                                    kernel::AbstractSPHKernel,
                                    x_pixels::Integer, y_pixels::Integer )
                                                
Calculates the kernel- and geometric weights of the pixels a particle contributes to.
"""
function calculate_weights( wk::Vector{Float64}, V::Vector{Float64},
                            iMin::Integer, iMax::Integer, 
                            jMin::Integer, jMax::Integer,
                            kMin::Integer, kMax::Integer,
                            x::T, y::T, z::T, 
                            hsml::T, hsml_inv::T,
                            kernel::AbstractSPHKernel,
                            x_pixels::Integer, y_pixels::Integer ) where T

    # storage variables for count operations
    n_distr_pix  = 0
    n_tot_pix    = 0
    distr_weight = 0.0
    distr_volume = 0.0

    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            for k = kMin:kMax
                z_dist, dz = get_x_dx(z, hsml, k)

                # current (flattened) index
                idx = calculate_index(i, j, k, x_pixels, y_pixels)

                # contributing volume
                dxdydz = dx * dy * dz

                # distance from pixel center in units of hsml
                u = get_d_hsml(x_dist, y_dist, z_dist, hsml_inv)

                V[idx], wk[idx], 
                distr_volume, distr_weight, 
                n_tot_pix, n_distr_pix = get_weight_per_pixel(distr_volume, distr_weight, 
                                                              n_tot_pix, n_distr_pix, 
                                                              dxdydz, u, hsml_inv, kernel)
            end # k
        end # j
    end # i

    # if particle contributes to pixels
    # but does not overlap with any pixel center
    if iszero(distr_weight)
        
        n_distr_pix = n_tot_pix

        # write full particle quantity into the pixel
        wk[1:n_tot_pix] .= 1.0
        
        # the weight is normalized by the pixel area
        if !iszero(distr_volume)
            weight_per_pix = n_distr_pix / distr_volume
        else
            weight_per_pix = 1
        end
    else 
        weight_per_pix = n_distr_pix / distr_weight
    end

    return wk, V, n_distr_pix, weight_per_pix
end


"""
    get_quantities_3D( pos, weight, hsml, 
                            rho, m, len2pix::T) where T

Helper function to get 
"""
function get_quantities_3D( pos, weight, hsml, 
                            rho, m, len2pix::T) where T

        hsml     = T(hsml * len2pix)
        hsml_inv = T(1.0/hsml)

        rho    /= T(len2pix)^3
        volume  = T(m / rho)

    return T.(pos), T(weight), hsml, hsml_inv, volume
end


"""
   cic_mapping_3D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::AbstractSPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 3D grid.
"""

function cic_mapping_3D( Pos, HSML, 
        M, Rho, 
        Bin_Q, Weights;
        param::mappingParameters, kernel::AbstractSPHKernel,
        show_progress::Bool=false,
        calc_mean=false )

    N = size(M,1)  # number of particles

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2] * param.Npixels[3]

    # allocate image array
    image = zeros(Float64, N_distr, 2)

    # allocate arrays for weights
    wk = zeros(Float64, N_distr)
    # storage array for mapped volume
    V  = Vector{Float64}(undef, N_distr)

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

        # convert to pixel units
        _pos, los_weight, hsml, hsml_inv, vol = get_quantities_3D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        # simplify position quantities for performance
        x, y, z = get_xyz( _pos, param)
        
        # calculate relevant pixel range
        iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
        jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )
        kMin, kMax = pix_index_min_max( z, hsml, param.Npixels[3] )
        
        # calculate all relevant quantities
        wk, V, N, weight_per_pix  = calculate_weights(wk, V, 
                                                            iMin, iMax, jMin, jMax,
                                                            kMin, kMax,
                                                            x, y, z, 
                                                            hsml, hsml_inv, 
                                                            kernel,
                                                            param.Npixels[1], 
                                                            param.Npixels[2])

        # normalisation factors for pixel contribution
        kernel_norm = vol / N
        volume_norm = kernel_norm * weight_per_pix * los_weight

        # loop over all contributing pixels
        @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax

            idx = calculate_index( i, j, k, 
                                    param.Npixels[1],
                                    param.Npixels[2] )

            # compute pixel weight 
            pix_weight = wk[idx] * V[idx] * volume_norm

            if !iszero(pix_weight)
                update_image!(image, idx, pix_weight, bin_q)
            end

            # store mass computed from grid cells
            grid_mass += Rho[p] * wk[idx] * V[idx] / param.len2pix^3
            
        end # i, j, k


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