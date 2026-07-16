"""
    function calculate_weights_3D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    kMin::Integer, kMax::Integer,
                                    x::Real, y::Real, z::Real, 
                                    hsml::Real, hsml_inv::Real,
                                    kernel::AbstractSPHKernel,
                                    y_pixels::Integer, z_pixels::Integer )

Calculates the kernel- and geometric weights of the pixels a particle contributes to.
`y_pixels`/`z_pixels` are the pixel counts along the two inner-loop axes, i.e.
the strides of the flattened 3D cube.
"""
function calculate_weights( wk::Vector{Float64}, V::Vector{Float64},
                            iMin::Integer, iMax::Integer,
                            jMin::Integer, jMax::Integer,
                            kMin::Integer, kMax::Integer,
                            x::T, y::T, z::T,
                            hsml::T, hsml_inv::T,
                            kernel::AbstractSPHKernel,
                            y_pixels::Integer, z_pixels::Integer ) where T

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
                idx = calculate_index(i, j, k, y_pixels, z_pixels)

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
        @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax
            idx = calculate_index(i, j, k, y_pixels, z_pixels)
            wk[idx] = 1.0
        end
        
        # the weight is normalized by the pixel volume
        if !iszero(distr_volume)
            weight_per_pix = n_distr_pix / distr_volume
        else
            # particle touches no in-grid pixel (clamped empty range);
            # use 1.0 (not Int 1) for type stability. The caller guards
            # against n_distr_pix == 0 so this particle is skipped and
            # does not corrupt the cube with Inf/NaN.
            weight_per_pix = 1.0
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
        wk, V, n_distr_pix, weight_per_pix  = calculate_weights(wk, V,
                                                            iMin, iMax, jMin, jMax,
                                                            kMin, kMax,
                                                            x, y, z,
                                                            hsml, hsml_inv,
                                                            kernel,
                                                            param.Npixels[2],
                                                            param.Npixels[3])

        # particle touches no in-grid pixel (e.g. clamped out at the edge or
        # off-grid): skip it instead of dividing vol / 0 = Inf and writing
        # NaNs into the cube. The off-grid fraction is lost by design for a
        # sub-region map.
        if iszero(n_distr_pix)
            if show_progress
                next!(P)
            end
            continue
        end

        # normalisation factors for pixel contribution
        kernel_norm = vol / n_distr_pix
        volume_norm = kernel_norm * weight_per_pix * los_weight * param.len2pix

        # loop over all contributing pixels
        @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax

            idx = calculate_index( i, j, k,
                                    param.Npixels[2],
                                    param.Npixels[3] )

            # compute pixel weight 
            pix_weight = wk[idx] * V[idx] * volume_norm

            if !iszero(pix_weight)
                update_image!(image, idx, pix_weight, bin_q)
            end

        end # i, j, k

        # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    # report grid-vs-particle mass conservation (coverage diagnostic)
    if show_progress
        mass_conservation_report(image, M, Rho, Weights, param.len2pix; ndim=3)
    end

    return image

end # function 