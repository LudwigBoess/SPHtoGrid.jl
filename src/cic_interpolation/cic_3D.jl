"""
    get_d_hsml_3D( dx::Real, dy::Real, dz::Real,
                   hsml_inv::Real )

Computes the distance in 3D to the pixel center in units of the kernel support.
"""
function get_d_hsml( dx::T, dy::T, dz::T,
                     hsml_inv::T) where T
    √( dx*dx + dy*dy + dz*dz ) * hsml_inv
end



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
@fastmath function calculate_weights( wk::Vector{Float64}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    kMin::Integer, kMax::Integer,
                                    x::T, y::T, z::T, 
                                    hsml::T, hsml_inv::T,
                                    kernel::AbstractSPHKernel,
                                    x_pixels::Integer, y_pixels::Integer ) where T

    is_undersampled = false

    if hsml <= 1.0
        is_undersampled = true
    end

    distr_weight = 0.0
    n_distr_pix  = 0

    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            for k = kMin:kMax
                z_dist, dz = get_x_dx(z, hsml, k)

                idx = calculate_index(i, j, k, x_pixels, y_pixels)

                dxdydz        = dx * dy * dz

                if is_undersampled

                    wk[idx]       = dxdydz
                    distr_weight += dxdydz
                    n_distr_pix  += 1

                    continue

                else # is_undersampled

                    u = get_d_hsml(x_dist, y_dist, z_dist, hsml_inv)

                    if u <= 1.0

                        wk[idx]       =  𝒲(kernel, u, hsml_inv)
                        wk[idx]      *= dxdydz
                        distr_weight += wk[idx]
                        n_distr_pix  += 1
                    
                    else
                        wk[idx] = 0.0
                    end # u < 1.0

                end # is_undersampled
            end # k
        end # j
    end # i

    return wk, n_distr_pix, distr_weight
end


"""
    function calculate_index(i::Integer, j::Integer, x_pixels::Integer)

Calculates the index of a flattened 3D image array.
"""
function calculate_index(  i::T, j::T, k::T, x_pixels::T, y_pixels::T) where T
    return floor(T, i * x_pixels + j * y_pixels + k) + 1
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

        volume  = T(4π/3) * hsml^3
        rho    /= T(len2pix)^3
        dz      = T(m / rho / volume)

    return T.(pos), T(weight), hsml, hsml_inv, volume, dz
end


"""
   sphMapping_3D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::AbstractSPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 3D grid.
"""
function sphMapping_3D( Pos::Array{T}, HSML::Array{T}, 
                        M::Array{T}, Rho::Array{T}, 
                        Bin_Q::Array{T}, Weights::Array{T};
                        param::mappingParameters, kernel::AbstractSPHKernel,
                        show_progress::Bool=false ) where T

    N = size(M,1)  # number of particles

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2] * param.Npixels[3]

    image = zeros(Float64, N_distr,2)

    halfXsize = param.halfsize[1]
    halfYsize = param.halfsize[2]
    halfZsize = param.halfsize[3]

    if param.periodic
        k_start = 0
    else
        k_start = 7
    end

    # allocate arrays for weights
    wk = zeros(Float64, N_distr)

    # allocate array for positions
    #pos = Vector{eltype(Pos[1,1])}(undef, 3)

    if show_progress
        P = Progress(N)
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Float64(Bin_Q[p])

        if bin_q == 0.0
            continue
        end

        _pos, weight, hsml, hsml_inv, volume, dz = get_quantities_3D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        for k_periodic = k_start:7

            x, y, z, skip_k = get_xyz( _pos, HSML[p], k_periodic, param)

            if skip_k
                continue
            end
            
            # calculate relevant pixel range
            iMin, iMax = pix_index_min_max( x, hsml, param.Npixels[1] )
            jMin, jMax = pix_index_min_max( y, hsml, param.Npixels[2] )
            kMin, kMax = pix_index_min_max( z, hsml, param.Npixels[3] )

            wk, n_distr_pix, distr_weight = calculate_weights( wk, 
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

            weight_per_pix = n_distr_pix / distr_weight
            kernel_norm    = volume / n_distr_pix
            volume_norm    = kernel_norm * weight_per_pix * weight * dz

            @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax

                idx = calculate_index( i, j, k, 
                                       param.Npixels[1],
                                       param.Npixels[2] )

                image[idx,1], image[idx,2] = update_image( image[idx,1], image[idx,2], 
                                                           wk[idx], bin_q, volume_norm)
                
            end # i, j, k

        end # k_periodic


         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    return image

end # function 