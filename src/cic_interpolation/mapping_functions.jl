"""
            Functions for sph mapping to grid.


            The property conservation is based on Smac by Dolag et. al. 2005:
            https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...29D/abstract

    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

using ProgressMeter
using SPHKernels

"""
    check_in_image(x::Real, y::Real, z::Real, hsml::Real,
                                halfXsize::Real, halfYsize::Real, halfZsize::Real)

Checks if a particle is in the image frame.
"""
@inline function check_in_image(x::Real, y::Real, z::Real, hsml::Real,
                                halfXsize::Real, halfYsize::Real, halfZsize::Real)

    if (( x + hsml ) > halfXsize ||  ( x - hsml ) < -halfXsize ||
        ( y + hsml ) > halfYsize ||  ( y - hsml ) < -halfYsize ||
        ( z + hsml ) > halfZsize ||  ( z - hsml ) < -halfZsize )

        return false
    else
        return true
    end
end

@inline function check_in_image(pos::Array{T}, hsml::T,
                                halfsize::Array{Float64}) where T

    if (( pos[1] + hsml ) > halfsize[1] ||  ( pos[1] - hsml ) < -halfsize[1] ||
        ( pos[2] + hsml ) > halfsize[2] ||  ( pos[2] - hsml ) < -halfsize[2] ||
        ( pos[3] + hsml ) > halfsize[3] ||  ( pos[3] - hsml ) < -halfsize[3] )

        return false
    else
        return true
    end
end


"""
        Positions and distances
"""

"""
    find_position_periodic( pos::Array{<:Real}, k::Integer, boxsize::Real)

Performs a periodic mapping of the particle position.
"""
@inline function find_position_periodic( pos::Array{T}, k::Integer, boxsize::Float64) where T
    
    # re-define here to avoid type-instability
    box   = eltype(pos[1])(boxsize)
    
    x = (k & 0x1) == 0 ? pos[1] : (pos[1] > 0 ? 
			pos[1] - box : pos[1] + box)
	
	y = (k & 0x2) == 0 ? pos[2] : (pos[2] > 0 ? 
			pos[2] - box : pos[2] + box)

	z = (k & 0x4) == 0 ? pos[3] : (pos[3] > 0 ? 
            pos[3] - box : pos[3] + box)
            
    return x, y, z
end


@inline function add_subtr_box(pos::T, boxsize::T) where T
    if pos > 0.0
        return pos - boxsize
    else
        return pos + boxsize
    end
end

@inline function find_position_periodic!( _pos::Array{T}, pos::Array{T}, k::Integer, boxsize::Float64) where T
    
    # re-define here to avoid type-instability
    box   = eltype(pos[1])(boxsize)

    _pos[1] = (k & 0x1) == 0 ? pos[1] : add_subtr_box(pos[1], box)
	
	_pos[2] = (k & 0x2) == 0 ? pos[2] : add_subtr_box(pos[2], box)

	_pos[3] = (k & 0x4) == 0 ? pos[3] : add_subtr_box(pos[3], box)

    return _pos
end

"""
    get_d_hsml_2D( dx::Real, dy::Real,
                   hsml_inv::Real )

Computes the distance in 2D to the pixel center in units of the kernel support.
"""
@inline function get_d_hsml_2D(dx::T, dy::T, hsml_inv::T) where T
    @fastmath sqrt( dx*dx + dy*dy ) * hsml_inv
end

"""
    get_d_hsml_3D( dx::Real, dy::Real, dz::Real,
                   hsml_inv::Real )

Computes the distance in 3D to the pixel center in units of the kernel support.
"""
@inline function get_d_hsml_3D(dx::T, dy::T, dz::T,
                               hsml_inv::T) where T
    @fastmath sqrt( dx*dx + dy*dy + dz*dz ) * hsml_inv
end

"""
    function get_dxyz(x::Real, hsml::Real, i::Integer)

Calculates the extent of the current particle size in units of pixels.
"""
@inline @fastmath function get_dxyz(x::Real, hsml::Real, i::Integer)
    return min(x + hsml, i + 1) - max(x - hsml, i)
end

"""
    function get_xyz( pos::Vector{<:Real}, hsml::Real, k::Integer,
                      len2pix::Real, x_pixels::Integer, y_pixels::Integer, z_pixels::Integer,
                      boxsize::Real, periodic::Bool,
                      halfXsize::Real, halfYsize::Real, halfZsize::Real,
                      Ndim::Integer)

Calculates `x, y, z` position in units of pixels and performs periodic mapping, if required.
"""
@inline @fastmath function get_xyz( pos::Vector{T}, hsml::T, k::Integer,
                                   len2pix::Float64, x_pixels::Int64, y_pixels::Int64, z_pixels::Int64,
                                   boxsize::Float64, periodic::Bool,
                                   halfXsize::Float64, halfYsize::Float64, halfZsize::Float64,
                                   Ndim::Integer) where T
    
    if periodic
        x, y, z = find_position_periodic(pos, k, boxsize)

        # check if the particle is still in the image
        if (( x - hsml ) > halfXsize ||  ( x + hsml ) < -halfXsize ||
            ( y - hsml ) > halfYsize ||  ( y + hsml ) < -halfYsize ||
            ( z - hsml ) > halfZsize ||  ( z + hsml ) < -halfZsize )

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

    x *= T(len2pix)
    y *= T(len2pix)
    x += T(0.5 * x_pixels)  
    y += T(0.5 * y_pixels)

    if Ndim == 2
        return x, y, false
    elseif Ndim == 3
        z *= T(len2pix)
        z += T(0.5 * z_pixels)
        return x, y, z, false
    end
end

@inline @fastmath function pos2pix( _pos::Vector{<:Real}, pos::Vector{<:Real}, hsml::Real, k::Integer,
                                   par::mappingParameters, Ndim::Integer)
    
    if par.periodic
        find_position_periodic!(_pos, pos, k, par.boxsize)

        # check if the particle is still in the image
        if !check_in_image(_pos, hsml, par.halfsize)
            return _pos, true
        end
    else
        _pos = pos
    end


    @inbounds for i = 1:Ndim
        _pos[i] *= par.len2pix
        _pos[i] += 0.5 * par.Npixels[i]
    end  

    return _pos, false
end


"""
            Indices
"""

"""
    function get_ijk_min_max( x::Real, hsml::Real,
                              x_pixels::Integer)

Calculates the minimum and maximum pixel to which a particle contributes.
"""
@inline @fastmath function get_ijk_min_max(x::T, hsml::T,
                                           x_pixels::Int64) where T

    iMin = max(floor(Integer, x - hsml), 0)
    iMax = min(floor(Integer, x + hsml), x_pixels-1)

    return iMin, iMax
end

"""
    function calculate_index_2D(i::Integer, j::Integer, x_pixels::Integer)

Calculates the index of a flattened 2D image array.
"""
@inline @fastmath function calculate_index_2D(i::T, j::T, x_pixels::T) where T
    return floor(T, i * x_pixels + j ) + 1
end

"""
    function calculate_index_3D(i::Integer, j::Integer, x_pixels::Integer)

Calculates the index of a flattened 3D image array.
"""
@inline @fastmath function calculate_index_3D(  i::T, j::T, k::T,
                                                x_pixels::T, y_pixels::T) where T
    return floor(T, i * x_pixels + j * y_pixels + k) + 1
end


"""
            Weights
"""


@inline @fastmath function calc_kernel_and_weights(  is_undersampled::Bool,
                    n_distr_pix::Integer, distr_weight::Float64, distr_area::Float64,
                    dx::T, dy::T, x_dist::T, y_dist::T, 
                    hsml_inv::T, kernel::SPHKernel) where T<:Real

    dxdy = dx * dy
    distr_area += dxdy

    if is_undersampled

        wk = Float64(dxdy)
        distr_weight += dxdy
        n_distr_pix  += 1

    else # is_undersampled

        u = get_d_hsml_2D(x_dist, y_dist, hsml_inv)

        if u <= 1.0

            wk = Float64(kernel_value_2D(kernel, u, hsml_inv))
            wk *= Float64(dxdy)
            distr_weight += wk
            n_distr_pix += 1
        else
            wk = Float64(0.0)
        end # u < 1.0

    end # is_undersampled

    return wk, distr_weight, distr_area, n_distr_pix
end

@inline @fastmath function get_x_dx(x::T, hsml::T, i::Integer) where T

    dx = get_dxyz(x, hsml, i)
    x_dist = T(x - i - 0.5)

    return x_dist, dx
end

"""
    function calculate_weights_2D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                    kernel::SPHKernel,
                                    x_pixels::Integer )

Calculates the kernel- and geometric weights of the pixels a particle contributes to.
"""
@inline @fastmath function calculate_weights_2D(  wk::Vector{Float64}, 
                                                iMin::Integer, iMax::Integer, 
                                                jMin::Integer, jMax::Integer,
                                                x::Real, y::Real, hsml::Real, hsml_inv::Real,
                                                kernel::SPHKernel,
                                                x_pixels::Integer )

    is_undersampled = false

    if hsml <= 1.0
        is_undersampled = true
    end

    distr_weight::Float64 = 0.0
    distr_area::Float64   = 0.0
    n_distr_pix           = 0


    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        @inbounds for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            idx = calculate_index_2D(i, j, x_pixels)

            wk[idx], distr_weight, distr_area, n_distr_pix = calc_kernel_and_weights(is_undersampled, n_distr_pix, 
                                                                                     distr_weight, distr_area, dx, dy, 
                                                                                     x_dist, y_dist, hsml_inv, kernel)
        end # j
    end # i

    return wk, n_distr_pix, distr_weight, distr_area
end

"""
    function calculate_weights_3D(  wk::Array{<:Real,1}, 
                                    iMin::Integer, iMax::Integer, 
                                    jMin::Integer, jMax::Integer,
                                    kMin::Integer, kMax::Integer,
                                    x::Real, y::Real, z::Real, 
                                    hsml::Real, hsml_inv::Real,
                                    kernel::SPHKernel,
                                    x_pixels::Integer, y_pixels::Integer )
                                                
Calculates the kernel- and geometric weights of the pixels a particle contributes to.
"""
@inline @fastmath function calculate_weights_3D(  wk::Vector{Float64}, 
                                                iMin::Integer, iMax::Integer, 
                                                jMin::Integer, jMax::Integer,
                                                kMin::Integer, kMax::Integer,
                                                x::T, y::T, z::T, 
                                                hsml::T, hsml_inv::T,
                                                kernel::SPHKernel,
                                                x_pixels::Integer, y_pixels::Integer ) where T

    is_undersampled = false

    if hsml <= 1.0
        is_undersampled = true
    end

    distr_weight::Float64 = 0.0
    distr_volume::Float64 = 0.0
    n_distr_pix  = 0


    @inbounds for i = iMin:iMax
        x_dist, dx = get_x_dx(x, hsml, i)

        @inbounds for j = jMin:jMax
            y_dist, dy = get_x_dx(y, hsml, j)

            @inbounds for k = kMin:kMax
                z_dist, dz = get_x_dx(z, hsml, k)

                idx = calculate_index_3D(i, j, k, x_pixels, y_pixels)

                dxdydz        = dx * dy * dz
                distr_volume += Float64(dxdydz)

                if is_undersampled

                    wk[idx] = Float64(dxdydz)
                    distr_weight += Float64(dxdydz)
                    n_distr_pix  += 1

                    continue

                else # is_undersampled

                    u = get_d_hsml_3D(x_dist, y_dist, z_dist, hsml_inv)

                    if u <= 1.0

                        wk[idx] =  Float64(kernel_value_3D(kernel, u, hsml_inv))
                        wk[idx] *= Float64(dxdydz)
                        distr_weight += wk[idx]
                        n_distr_pix += 1
                    
                    else
                        wk[idx] = Float64(0.0)
                        continue
                    end # u < 1.0

                end # is_undersampled
            end # k
        end # j
    end # i

    return wk, n_distr_pix, distr_weight, distr_volume
end


"""
    function update_image( image::Real, w_image::Real, 
                           wk::Real, bin_q::Real, 
                           geometry_norm::Real )

Applies the different contributions to the image and the weight image.
"""
@inline @fastmath function update_image( image::Float64, w_image::Float64, 
                                         wk::Float64, bin_q::T, 
                                         geometry_norm::Float64 ) where T<:Real

    if wk > 0.0
        pix_weight = geometry_norm * wk
        image   += bin_q * pix_weight
        w_image += pix_weight
    end

    return image, w_image
end

@inline @fastmath function update_image( image::Array{<:Real}, 
                                         wk::Real, bin_q::Real, 
                                         geometry_norm::Real )

    if wk > 0.0
        pix_weight = geometry_norm * wk
        image[1] += bin_q * pix_weight
        image[2] += pix_weight
    end

    return image[1], image[2]
end

@inline @fastmath function update_image!( idx::Integer, 
                                          image::Array{Float64}, 
                                          wk::Array{Float64}, bin_q::T, 
                                          geometry_norm::Float64 ) where T

    if wk[idx] > 0.0
        pix_weight = geometry_norm * wk[idx]
        image[idx,1] += bin_q * pix_weight
        image[idx,2] += pix_weight
    end

    image
end



"""
        Helper functions for variables
"""


@inline @fastmath function get_quantities_2D(pos::Array{T}, weight::T, hsml::T, 
                                          rho::T, m::T, len2pix::Float64) where T
        
    hsml *= T(len2pix)
    hsml_inv = T(1.0)/hsml
    area  = π * hsml^2

    rho *= T(1.0/len2pix*len2pix*len2pix)
    dz  = m / rho / area

    return pos, weight, hsml, hsml_inv, area, m, rho, dz
end

@inline @fastmath function get_quantities_3D(pos::Array{T}, weight::T, hsml::T, 
                                          rho::T, m::T, len2pix::Float64) where T
        hsml *= T(len2pix)
        hsml_inv = T(1.0)/hsml

        volume  = T(4π/3) * hsml^3
        rho /= T(len2pix)^3
        dz  = m / rho / volume

    return pos, weight, hsml, hsml_inv, volume, m, rho, dz
end






"""
   sphMapping_2D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::SPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 2D grid.
"""
function sphMapping_2D( Pos::Array{T}, HSML::Array{T}, 
                        M::Array{T}, Rho::Array{T}, 
                        Bin_Q::Array{T}, Weights::Array{T};
                        param::mappingParameters, kernel::SPHKernel,
                        show_progress::Bool=false ) where T

    N = size(M,1)  # number of particles
    
    N_distr = param.Npixels[1] * param.Npixels[2]

    image = zeros(Float64, N_distr, 2)

    # store this here for performance increase
    halfXsize = param.halfsize[1]
    halfYsize = param.halfsize[2]
    halfZsize = param.halfsize[3]

    if param.periodic
        k_start = 0
    else
        k_start = 7
    end

    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[2] * param.Npixels[1]

    # allocate arrays for weights
    wk = zeros(Float64, N_distr)

    # allocate array for positions
    #pos = Vector{eltype(Pos[1,1])}(undef, 3)

    if show_progress
        P = Progress(N)
        idx_p = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Bin_Q[p]

        if bin_q == 0.0
            continue
        end

        _pos, weight, hsml, hsml_inv, area, m, rho, dz = get_quantities_2D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        for k = k_start:7


            x, y, skip_k = get_xyz( _pos, HSML[p], k, 
                                   param.len2pix, param.Npixels[1], param.Npixels[2], param.Npixels[3], 
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

            weight_per_pix  = n_distr_pix / distr_weight
            kernel_norm     = area / n_distr_pix
            area_norm       = Float64(kernel_norm * weight_per_pix * weight * dz)

            @fastmath @inbounds for i = iMin:iMax, j = jMin:jMax
                idx = calculate_index_2D(i, j, param.Npixels[1])

                image[idx,1], image[idx,2] = update_image( image[idx,1], image[idx,2], 
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

    return image

end # function 


"""
   sphMapping_2D( Pos::Array{<:Real}, HSML::Array{<:Real}, 
                  M::Array{<:Real}, Rho::Array{<:Real}, 
                  Bin_Q::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));
                  param::mappingParameters, kernel::SPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 3D grid.
"""
function sphMapping_3D( Pos::Array{T}, HSML::Array{T}, 
                        M::Array{T}, Rho::Array{T}, 
                        Bin_Q::Array{T}, Weights::Array{T};
                        param::mappingParameters, kernel::SPHKernel,
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
        k_start = 8
    end

    # allocate arrays for weights
    wk = zeros(Float64, N_distr)

    # allocate array for positions
    #pos = Vector{eltype(Pos[1,1])}(undef, 3)

    if show_progress
        P = Progress(N)
        idx_p = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Bin_Q[p]

        if bin_q == 0.0
            continue
        end

        _pos, weight, hsml, hsml_inv, volume, m, rho, dz = get_quantities_3D(Pos[:,p], Weights[p], HSML[p], Rho[p], M[p], param.len2pix)

        for k_periodic = k_start:8


            x, y, z, skip_k = get_xyz( _pos, HSML[p], k_periodic, 
                                   param.len2pix, param.Npixels[1], param.Npixels[2], param.Npixels[3], 
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

                image[idx,1], image[idx,2] = update_image( image[idx,1], image[idx,2], 
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

    return image

end # function 