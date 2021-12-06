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
function calculate_weights_splash(wk::Vector{Float64},
    iMin::Integer, iMax::Integer,
    jMin::Integer, jMax::Integer,
    kMin::Integer, kMax::Integer,
    x::T, y::T, z::T,
    hsml::T, hsml_inv::T,
    kernel::AbstractSPHKernel,
    x_pixels::Integer, y_pixels::Integer) where {T}

    is_undersampled = false

    if hsml <= 1
        is_undersampled = true
    end

    n_distr_pix = 0

    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            for k = kMin:kMax
                z_dist, dz = get_x_dx(z, hsml, k)

                idx = calculate_index(i, j, k, x_pixels, y_pixels)

                if !is_undersampled

                    u = get_d_hsml(x_dist, y_dist, z_dist, hsml_inv)

                    # check if inside the kernel
                    if u <= 1
                        # evaluate kernel
                        wk[idx] = 𝒲(kernel, u, hsml_inv)
                        n_distr_pix += 1
                    else
                        wk[idx] = 0
                    end # u < 1.0

                end # is_undersampled
            end # k
        end # j
    end # i

    return wk, n_distr_pix
end



"""
   splash_mapping_3D( Pos, HSML, Bin_Q;
                  param::mappingParameters, kernel::AbstractSPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 3D grid.
"""
function splash_mapping_3D(Pos, HSML, Bin_Q;
    param::mappingParameters, kernel::AbstractSPHKernel,
    show_progress::Bool = false)

    N = size(HSML, 1)  # number of particles

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2] * param.Npixels[3]

    image = zeros(Float64, N_distr)

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
    
        if bin_q == 0.0
            continue
        end
    
        _pos, hsml, hsml_inv = get_quantities_splash(Pos[:, p], HSML[p], param.len2pix)
    
        for k_periodic = k_start:7
        
            x, y, z, skip_k = get_xyz(_pos, HSML[p], k_periodic, param)
        
            if skip_k
                continue
            end
        
            # calculate relevant pixel range
            iMin, iMax = pix_index_min_max(x, hsml, param.Npixels[1])
            jMin, jMax = pix_index_min_max(y, hsml, param.Npixels[2])
            kMin, kMax = pix_index_min_max(z, hsml, param.Npixels[3])
        
            wk, n_distr_pix = calculate_weights_splash(wk,
                iMin, iMax,
                jMin, jMax,
                kMin, kMax,
                x, y, z,
                hsml, hsml_inv,
                kernel,
                param.Npixels[1],
                param.Npixels[2])
        
            # skip if the particle is not contained in the image 
            # ( should only happen with periodic boundary conditions )
            if n_distr_pix == 0
                continue
            end
        
            @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax
        
                idx = calculate_index(i, j, k,
                    param.Npixels[1],
                    param.Npixels[2])
        
                image[idx] = wk[idx] * bin_q
        
            end # i, j, k
        end # k_periodic
    
        # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    return image

end # function 