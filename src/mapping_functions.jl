"""
            Functions for sph mapping to grid.

            sphMapping_2D without property conservation works like SPLASH by
            Daniel Price: http://users.monash.edu.au/~dprice/splash/.

            The property conservation is based on Smac by Dolag et. al. 2005:
            https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...29D/abstract

    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

using ProgressMeter
using Unitful
using SPHKernels

"""
    get_d_hsml_2D( dx::Real, dy::Real,
                   hsml_inv::Real )

Computes the distance in 2D to the pixel center in units of the kernel support.
"""
@inline function get_d_hsml_2D(dx::Real, dy::Real, hsml_inv::Real)
    @fastmath sqrt( dx*dx + dy*dy ) * hsml_inv
end

"""
    get_d_hsml_3D( dx::Real, dy::Real, dz::Real,
                   hsml_inv::Real )

Computes the distance in 3D to the pixel center in units of the kernel support.
"""
@inline function get_d_hsml_3D(dx::Real, dy::Real, dz::Real,
                               hsml_inv::Real)
    @fastmath sqrt( dx*dx + dy*dy + dz*dz ) * hsml_inv
end

"""
    check_in_image( pos::Real, hsml::Real,
                    minCoords::Real, maxCoords::Real )

Checks if a particle is in the image frame.
"""
@inline function check_in_image(pos::Real, hsml::Real,
                                minCoords::Real, maxCoords::Real)

    if ( (minCoords - hsml) < pos < (maxCoords + hsml) )
        return true
    else
        return false
    end
end

"""
    find_min_pixel( pos::Real, hsml::Real,
                    minCoords::Real,
                    pixsize_inv::Real )

Computes the minimum pixel to which the particle contributes.
"""
@inline function find_min_pixel(pos::Real, hsml::Real,
                                minCoords::Real,
                                pixsize_inv::Real)

    @fastmath pix = floor(Int64, (pos - hsml - minCoords) * pixsize_inv )

    return max(pix, 1)
end

"""
    find_max_pixel( pos::Real, hsml::Real,
                    minCoords::AbstractFloat, pixsize_inv::Real, 
                    max_pixel::Integer )

Computes the maximum pixel to which the particle contributes.
"""
@inline function find_max_pixel(pos::Real, hsml::Real,
                                minCoords::AbstractFloat, pixsize_inv::Real, 
                                max_pixel::Integer)

    @fastmath pix = floor(Int64, (pos + hsml  - minCoords) * pixsize_inv )

    return min(pix, max_pixel)
end


@inline function find_position_periodic( pos::Array{<:Real}, k::Integer, boxsize::Real)
    
    x = (k & 0x1) == 0 ? pos[1] : (pos[1] > 0 ? 
			pos[1] - boxsize : pos[1] + boxsize)
	
	y = (k & 0x2) == 0 ? pos[2] : (pos[2] > 0 ? 
			pos[2] - boxsize : pos[2] + boxsize)

	z = (k & 0x4) == 0 ? pos[3] : (pos[3] > 0 ? 
            pos[3] - boxsize : pos[3] + boxsize)
            
    return x, y, z
end


@inline @fastmath function get_dxyz(x::Real, hsml::Real, i::Integer)
    return min(x + hsml, i + 1) - max(x - hsml, i)
end

@inline @fastmath function calculate_index_2D(i::Integer, j::Integer, x_pixels::Integer)
    return floor(Integer, i * x_pixels + j ) + 1
end

@inline @fastmath function calculate_index_3D(  i::Integer, j::Integer, k::Integer,
                                                x_pixels::Integer, y_pixels::Integer)
    return floor(Integer, i * x_pixels + j * y_pixels + k) + 1
end

@inline @fastmath function calculate_weights_2D(  wk::Array{<:Real,1}, 
                                                iMin::Integer, iMax::Integer, 
                                                jMin::Integer, jMax::Integer,
                                                x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                                kernel::SPHKernel,
                                                x_pixels::Integer )

    is_undersampled = false

    if hsml <= 1.0
        is_undersampled = true
    end

    distr_weight = 0.0
    distr_area   = 0.0
    n_distr_pix = 0


    @inbounds for i = iMin:iMax
        dx = get_dxyz(x, hsml, i)
        x_dist = x - i - 0.5

        @inbounds for j = jMin:jMax
            dy = get_dxyz(y, hsml, j)

            idx = calculate_index_2D(i, j, x_pixels)

            dxdy = dx * dy
            distr_area += dxdy

            if is_undersampled

                wk[idx] = dxdy
                distr_weight += dxdy
                n_distr_pix  += 1

                continue

            else # is_undersampled

                y_dist = y - j - 0.5

                u = get_d_hsml_2D(x_dist, y_dist, hsml_inv)

                if u <= 1.0

                    wk[idx] = kernel_value_2D(kernel, u, hsml_inv)
                    wk[idx] *= dxdy
                    distr_weight += wk[idx]
                    n_distr_pix += 1
                
                else
                    wk[idx] = 0.0
                    continue
                end # u < 1.0

            end # is_undersampled
        end # j
    end # i

    return wk, n_distr_pix, distr_weight, distr_area
end

@inline @fastmath function calculate_weights_3D(  wk::Array{<:Real,1}, 
                                                iMin::Integer, iMax::Integer, 
                                                jMin::Integer, jMax::Integer,
                                                kMin::Integer, kMax::Integer,
                                                x::Real, y::Real, z::Real, 
                                                hsml::Real, hsml_inv::Real,
                                                kernel::SPHKernel,
                                                x_pixels::Integer, y_pixels::Integer )

    is_undersampled = false

    if hsml <= 1.0
        is_undersampled = true
    end

    distr_weight = 0.0
    distr_volume = 0.0
    n_distr_pix  = 0


    @inbounds for i = iMin:iMax
        dx = get_dxyz(x, hsml, i)
        x_dist = x - i - 0.5

        @inbounds for j = jMin:jMax
            dy = get_dxyz(y, hsml, j)
            y_dist = y - j - 0.5

            @inbounds for k = kMin:kMax
                dz = get_dxyz(z, hsml, k)

                idx = calculate_index_3D(i, j, k, x_pixels, y_pixels)

                dxdydz = dx * dy * dz
                distr_volume += dxdydz

                if is_undersampled

                    wk[idx] = dxdydz
                    distr_weight += dxdydz
                    n_distr_pix  += 1

                    continue

                else # is_undersampled

                    z_dist = z - k - 0.5

                    u = get_d_hsml_3D(x_dist, y_dist, z_dist, hsml_inv)

                    if u <= 1.0

                        wk[idx] = kernel_value_3D(kernel, u, hsml_inv)
                        wk[idx] *= dxdydz
                        distr_weight += wk[idx]
                        n_distr_pix += 1
                    
                    else
                        wk[idx] = 0.0
                        continue
                    end # u < 1.0

                end # is_undersampled
            end # k
        end # j
    end # i

    return wk, n_distr_pix, distr_weight, distr_volume
end

@inline @fastmath function get_ijk_min_max(x::Real, hsml::Real,
                                          x_pixels::Integer)

    iMin = max(floor(Integer, x - hsml), 0)
    iMax = min(floor(Integer, x + hsml) + 1, x_pixels-1)

    return iMin, iMax
end

@inline @fastmath function update_image( image::Real, w_image::Real, 
                                         wk::Real, bin_q::Real, 
                                         geometry_norm::Real )

    if wk > 0.0
        pix_weight = geometry_norm * wk
        image   += bin_q * pix_weight
        w_image += pix_weight
    end

    return image, w_image
end

@inline @fastmath function get_quantities_2D(pos::Array{<:Real}, weight::Real, hsml::Real, 
                                          rho::Real, m::Real, len2pix::Real)
        hsml *= len2pix
        hsml_inv = 1.0/hsml

        area  = π * hsml^2
        rho /= len2pix^3
        dz  = m / rho / area

    return pos, weight, hsml, hsml_inv, area, m, rho, dz
end

@inline @fastmath function get_quantities_3D(pos::Array{<:Real}, weight::Real, hsml::Real, 
                                          rho::Real, m::Real, len2pix::Real)
        hsml *= len2pix
        hsml_inv = 1.0/hsml

        volume  = 4π/3 * hsml^3
        rho /= len2pix^3
        dz  = m / rho / volume

    return pos, weight, hsml, hsml_inv, volume, m, rho, dz
end

@inline @fastmath function reduce_image_2D( image::Array{<:Real}, w_image::Array{<:Real},
                                         x_pixels::Integer, y_pixels::Integer)


    im_plot = zeros(x_pixels, y_pixels)
    k = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels
        im_plot[i,j] = image[k]
        if w_image[k] > 0.0
            im_plot[i, j] /= w_image[k]
        end
        k += 1
    end
    return im_plot
end

@inline @fastmath function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},
                                         x_pixels::Integer, y_pixels::Integer, z_pixels::Real)


    im_plot = zeros(x_pixels, y_pixels, z_pixels)
    m = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels, k = 1:z_pixels
        im_plot[i,j,k] = image[m]
        if w_image[m] > 0.0
            im_plot[i, j, k] /= w_image[m]
        end
        m += 1
    end
    return im_plot
end

@inline @fastmath function get_xyz( pos::Vector{<:Real}, hsml::Real, k::Integer,
                                   len2pix::Real, x_pixels::Integer, y_pixels::Integer, z_pixels::Integer,
                                   boxsize::Real, periodic::Bool,
                                   halfXsize::Real, halfYsize::Real, halfZsize::Real,
                                   Ndim::Integer)
    
    if periodic
        x, y, z = find_position_periodic(pos, k, boxsize)

        # check if the particle is still in the image
        if  ( x - hsml ) > halfXsize ||  ( x + hsml ) < -halfXsize ||
            ( y - hsml ) > halfYsize ||  ( y + hsml ) < -halfYsize ||
            ( z - hsml ) > halfZsize ||  ( z + hsml ) < -halfZsize

            if Ndim == 2
                return x, y, true
            elseif Ndim == 3
                return x, y, z, true
            end
        end
    else
        x = pos[1]
        y = pos[2]
        z = pos[3]
    end

    x *= len2pix
    y *= len2pix
    x += 0.5 * x_pixels  
    y += 0.5 * y_pixels

    if Ndim == 2
        return x, y, false
    elseif Ndim == 3
        z *= len2pix
        z += 0.5 * z_pixels
        return x, y, z, false
    end
end


"""
   sphMapping_2D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::SPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 2D grid.
"""
function sphMapping_2D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                        M::Array{<:Real}, Rho::Array{<:Real}, 
                        Bin_Q::Array{<:Real}, Weights::Array{<:Real};
                        param::mappingParameters, kernel::SPHKernel,
                        show_progress::Bool=false )

    N = length(M)  # number of particles
    image = zeros(param.Npixels[1] * param.Npixels[2])
    w_image = zeros(param.Npixels[1] * param.Npixels[2])

    # store this here for performance increase
    len2pix = 1.0/param.pixelSideLength

    halfXsize = 0.5 * param.x_size
    halfYsize = 0.5 * param.y_size
    halfZsize = 0.5 * param.z_size

    if param.periodic
        k_start = 0
    else
        k_start = 8
    end

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[2] * param.Npixels[1]

    # allocate arrays for weights
    wk = zeros(N_distr)

    # allocate array for positions
    #pos = Vector{eltype(Pos[1,1])}(undef, 3)

    if show_progress
        P = Progress(N)
        idx_p = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Bin_Q[p,1]

        if bin_q == 0.0
            continue
        end

        _pos, weight, hsml, hsml_inv, area, m, rho, dz = get_quantities_2D(Pos[p,:], Weights[p], HSML[p], Rho[p], M[p], len2pix)

        for k = k_start:8


            x, y, skip_k = get_xyz( _pos, HSML[p], k, 
                                   len2pix, param.Npixels[1], param.Npixels[2], param.Npixels[3], 
                                   param.boxsize, param.periodic,
                                   halfXsize, halfYsize, halfZsize, 2)

            if skip_k
                continue
            end
            
            # calculate relevant pixel range
            iMin, iMax = get_ijk_min_max( x, hsml, param.Npixels[1] )
            jMin, jMax = get_ijk_min_max( y, hsml, param.Npixels[2] )

            wk, n_distr_pix, distr_weight, distr_area = calculate_weights_2D(wk, iMin, iMax, jMin, jMax,
                                                                x, y, hsml, hsml_inv, 
                                                                kernel,
                                                                param.Npixels[1])
           
            # skip if the particle is not contained in the image 
            # ( should only happen with periodic boundary conditions )
            if n_distr_pix == 0
                continue
            end

            weight_per_pix = n_distr_pix / distr_weight
            kernel_norm = area / n_distr_pix
            area_norm = kernel_norm * weight_per_pix * weight * dz

            @fastmath @inbounds for i = iMin:iMax, j = jMin:jMax
                idx = calculate_index_2D(i, j, param.Npixels[1])

                image[idx], w_image[idx] = update_image( image[idx], w_image[idx], 
                                                         wk[idx], bin_q, area_norm)
                

            end # i, j 

        end # k


         # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
            #unlock(P_lock)
        end
    end # p

    return reduce_image_2D( image, w_image, 
                         param.Npixels[1], param.Npixels[2] )

end # function 


"""
   sphMapping_2D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::SPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 3D grid.
"""
function sphMapping_3D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                        M::Array{<:Real}, Rho::Array{<:Real}, 
                        Bin_Q::Array{<:Real}, Weights::Array{<:Real};
                        param::mappingParameters, kernel::SPHKernel,
                        show_progress::Bool=false )

    N = length(M)  # number of particles

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2] * param.Npixels[3]

    image = zeros(N_distr)
    w_image = zeros(N_distr)

    # store this here for performance increase
    len2pix = 1.0/param.pixelSideLength

    if param.periodic
        halfXsize = 0.5 * param.x_size
        halfYsize = 0.5 * param.y_size
        halfZsize = 0.5 * param.z_size
        k_start = 0
    else
        k_start = 8
    end

    # allocate arrays for weights
    wk = zeros(N_distr)

    # allocate array for positions
    #pos = Vector{eltype(Pos[1,1])}(undef, 3)

    if show_progress
        P = Progress(N)
        idx_p = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Bin_Q[p,1]

        if bin_q == 0.0
            continue
        end

        _pos, weight, hsml, hsml_inv, volume, m, rho, dz = get_quantities_3D(Pos[p,:], Weights[p], HSML[p], Rho[p], M[p], len2pix)

        for k_periodic = k_start:8


            x, y, z, skip_k = get_xyz( _pos, HSML[p], k_periodic, 
                                   len2pix, param.Npixels[1], param.Npixels[2], param.Npixels[3], 
                                   param.boxsize, param.periodic,
                                   halfXsize, halfYsize, halfZsize, 3)

            if skip_k
                continue
            end
            
            # calculate relevant pixel range
            iMin, iMax = get_ijk_min_max( x, hsml, param.Npixels[1] )
            jMin, jMax = get_ijk_min_max( y, hsml, param.Npixels[2] )
            kMin, kMax = get_ijk_min_max( z, hsml, param.Npixels[3] )

            wk, n_distr_pix, distr_weight, distr_volume = calculate_weights_3D( wk, 
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

            @fastmath @inbounds for i = iMin:iMax, j = jMin:jMax, k = kMin:kMax
                idx = calculate_index_3D(i, j, k, 
                                         param.Npixels[1],
                                         param.Npixels[2])

                image[idx], w_image[idx] = update_image( image[idx], w_image[idx], 
                                                         wk[idx], bin_q, volume_norm)
                
            end # i, j, k

        end # k_periodic


         # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
            #unlock(P_lock)
        end
    end # p

    return reduce_image_3D( image, w_image, 
                            param.Npixels[1], param.Npixels[2], param.Npixels[3] )

end # function 