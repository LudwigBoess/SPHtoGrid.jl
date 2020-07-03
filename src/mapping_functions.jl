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

@inline function find_position_periodic(pos::Vector{Float64}, k::Int64, boxsize::Float64)

    x = !(k & 0x1) ? pos[1] : ( pos[1] > 0.0 ? pos[1] - boxsize : pos[1] + boxsize)
    y = !(k & 0x2) ? pos[2] : ( pos[2] > 0.0 ? pos[2] - boxsize : pos[2] + boxsize)
    z = !(k & 0x4) ? pos[3] : ( pos[3] > 0.0 ? pos[3] - boxsize : pos[3] + boxsize)

    return [x, y, z]
end


@inline function get_d_hsml_2D(dx::Float64, dy::Float64, hsml_inv::Float64)
    sqrt( dx*dx + dy*dy ) * hsml_inv
end

@inline function get_d_hsml_3D(dx::Float64, dy::Float64, dz::Float64,
                               hsml_inv::Float64)
    sqrt( dx*dx + dy*dy + dz*dz ) * hsml_inv
end

@inline function check_in_image(pos::Float64, hsml::Float64,
                                minCoords::Float64, maxCoords::Float64)

    if ( (minCoords - hsml) <= pos <= (maxCoords + hsml) )
        return true
    else
        return false
    end
end

@inline function find_min_pixel(pos::Float64, hsml::Float64,
                                minCoords::Float64,
                                pixsize_inv::Float64)

    pix = floor(Int64, (pos - hsml - minCoords) * pixsize_inv )

    return max(pix, 1)
end

@inline function find_max_pixel(pos::Float64, hsml::Float64,
                                minCoords::AbstractFloat, pixsize_inv::Float64,
                                max_pixel::Integer)

    pix = floor(Int64, (pos + hsml - minCoords) * pixsize_inv )

    return min(pix, max_pixel)
end


function sphMapping_2D(Pos, HSML, M, ρ, Bin_Quant;
                       param::mappingParameters, kernel,
                       conserve_quantities::Bool=false,
                       show_progress::Bool=false)

    # if this is not a float it has units, which need to be stripped
    if !(typeof(Bin_Quant[1,1]) <: AbstractFloat)

        if show_progress
            @info "Stripping units..."
        end

        Pos       = ustrip(Pos)
        HSML      = ustrip(HSML)
        M         = ustrip(M)
        ρ         = ustrip(ρ)
        Bin_Quant = ustrip(Bin_Quant)

    end

    N = length(M)  # number of particles

    image = zeros(param.Npixels[1], param.Npixels[2])

    minCoords = [param.x_lim[1], param.y_lim[1], param.z_lim[1]]
    maxCoords = [param.x_lim[2], param.y_lim[2], param.z_lim[2]]

    #max_pixel = [length(param.x), length(param.y)]

    cen = 0.5 .* [ param.x_size, param.y_size, param.z_size ]

    particles_in_image = 0

    if show_progress
        P = Progress(N)
        idx = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    @inbounds for p = 1:N

        # save stuff from array to single variables
        pos  = Float64.(Pos[p,:])
        hsml = Float64(HSML[p])

        in_image = false

        @inbounds for dim = 1:3

            in_image = check_in_image(pos[dim], hsml,
                                      minCoords[dim], maxCoords[dim])

            # exit the loop if the particle is not in the image frame
            if !in_image
                break
            end
        end

        # only calculate the properties if the particle is in the image
        if in_image

            # shift particles
            @inbounds for dim = 1:3
                pos[dim] += cen[dim]
            end
            
            # store this here for performance increase
            pixsize_inv = 1.0/param.pixelSideLength

            particles_in_image += 1

            # save rest of variables
            bin_q       = Float64(Bin_Quant[p])
            hsml_inv    = Float64(1.0/hsml)
            m           = Float64(M[p])
            rho_inv     = Float64(1.0/ρ[p])

            pixmin = Vector{Int64}(undef,2)
            pixmax = Vector{Int64}(undef,2)

            @inbounds for dim = 1:2

                pixmin[dim] = find_min_pixel(pos[dim], hsml, minCoords[dim],
                                             pixsize_inv)

                pixmax[dim] = find_max_pixel(pos[dim], hsml, minCoords[dim],
                                             pixsize_inv, param.Npixels[dim])

            end

            if conserve_quantities

                # calculate pixel area
                xp1 = (pos[1] - hsml) * pixsize_inv
                xp2 = (pos[1] + hsml) * pixsize_inv
                yp1 = (pos[2] - hsml) * pixsize_inv
                yp2 = (pos[2] + hsml) * pixsize_inv

                pix_area = (xp2 - xp1) * (yp2 - yp1) * param.pixelArea

                d3 = (m * rho_inv) / ( pix_area * param.pixelSideLength)

                # number of pixels over which the particle is distributed
                N_distr = Int64((pixmax[1] - pixmin[1] + 1 ) * (pixmax[2] - pixmin[2] + 1))

                # allocate arrays for weights
                kernel_tab = zeros(N_distr)
                d1_tab     = zeros(N_distr)
                d2_tab     = zeros(N_distr)
                dx_tab     = zeros(N_distr)

                # weights table
                wit1    = 0.
                witd    = 0.
                wit2    = 0
                wit2tot = 0

                # pixel count for attributing weights to pixels
                N_count = 1

                # first loop to calculate weights
                @inbounds for i = pixmin[1]:pixmax[1]

                    dx = pos[1] - ( i - 0.5 ) * param.pixelSideLength
                    dimin = max(xp1, i - 1.0)
                    dimax = min(xp2, i)

                    @inbounds for j = pixmin[2]:pixmax[2]

                        dy = pos[2] - ( j - 0.5 ) * param.pixelSideLength
                        djmin = max(xp1, j - 1.0)
                        djmax = min(xp2, j)

                        d1 = dimax - dimin
                        d2 = djmax - djmin

                        # compute distance to pixel center in units of hsml
                        dx_tab[N_count] = get_d_hsml_2D(dx, dy, hsml_inv)
                        # update pixel value
                        wi = kernel_value_2D(kernel, dx_tab[N_count], hsml_inv)

                        wit1 += wi * d1 * d2
                        witd += d1 * d2

                        if wi > 0
                            wit2 += 1
                        end

                        # store data in arrays
                        d1_tab[N_count]     = d1
                        d2_tab[N_count]     = d2
                        kernel_tab[N_count] = wi

                        wit2tot += 1

                        N_count += 1

                    end # end x-loop
                end # end y-loop

                #fak = 1.0/( (pixmax[2] - pixmin[2]) * (pixmax[1] - pixmin[1]) )
                if wit1 > 0
                    fak = wit2/wit1
                else
                    wit2 = wit2tot
                    fak = wit2tot / witd
                end

                fak_hsml = pix_area/(wit2 * param.pixelArea)

                # reset the counter
                N_count = 1

            end # conserve_quantities


            # to reduce number of calculations. Only needed if quantity conservation
            # is switched off. Otherwise it is incapsulated in d3
            if !conserve_quantities
                bin_prefac = bin_q * m * rho_inv
            end

            # second loop to calculate value
            @inbounds for i = pixmin[1]:pixmax[1]

                @inbounds for j = pixmin[2]:pixmax[2]


                    if !conserve_quantities

                        # calculate simple distance to pixel center
                        dx = pos[1] - ( i - 0.5 ) * param.pixelSideLength
                        dy = pos[2] - ( j - 0.5 ) * param.pixelSideLength
                        distance_hsml = get_d_hsml_2D(dx, dy, hsml_inv)

                        # update pixel value
                        image[i, j] += bin_prefac * kernel_value_2D(kernel, distance_hsml, hsml_inv)
                    else

                        if wit1 <= 1
                            wi = fak * fak_hsml
                        else
                            wi = kernel_tab[N_count] * fak * fak_hsml
                        end
                        # update pixel value with weights
                        image[i, j] += bin_q * wi *
                                       d1_tab[N_count] * d2_tab[N_count] * d3
                        N_count += 1
                    end

                end # end x-loop
            end # end y-loop

        end # end check if in image


        # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx += 1
            ProgressMeter.update!(P, idx)
            #unlock(P_lock)
        end
    end

    if show_progress
        @info "Mapped $particles_in_image / $N particles."
    end

   return image
end


function sphMapping_3D(Pos, HSML, M, ρ, Bin_Quant;
                       param::mappingParameters, kernel::SPHKernel,
                       conserve_quantities::Bool=false,
                       show_progress::Bool=false)

    # if this is not a float it has units, which need to be stripped
    if !(typeof(Bin_Quant[1,1]) <: AbstractFloat)

        if show_progress
            @info "Stripping units..."
        end

        Pos       = ustrip(Pos)
        HSML      = ustrip(HSML)
        M         = ustrip(M)
        ρ         = ustrip(ρ)
        Bin_Quant = ustrip(Bin_Quant)

    end

    N = length(M)  # number of particles

    image = zeros(param.Npixels[1], param.Npixels[2])

    minCoords = [param.x_lim[1], param.y_lim[1], param.z_lim[1]]
    maxCoords = [param.x_lim[2], param.y_lim[2], param.z_lim[2]]

    cen = 0.5 .* [ param.x_size, param.y_size, param.z_size ]

    particles_in_image = 0

    if show_progress
        P = Progress(N)
        idx = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    @inbounds for p = 1:N


        # save stuff from array to single variables
        bin_q = Float64(Bin_Quant[p])
        pos   = Float64.(Pos[p,:])
        hsml  = Float64(HSML[p])

        in_image = false

        @inbounds for dim = 1:3

            in_image = check_in_image(pos[dim], hsml,
                                      minCoords[dim], maxCoords[dim])

            # exit the loop if the particle is not in the image frame
            if !in_image
                break
            end
        end

        # only calculate the properties if the particle is in the image
        if in_image

            # shift particles
            @inbounds for dim = 1:3
                pos[dim] += cen[dim]
            end

            # store this here for performance increase
            pixsize_inv = 1.0/param.pixelSideLength

            particles_in_image += 1

            # save rest of variables
            hsml_inv    = Float64(1.0/hsml)
            m           = Float64(M[p])
            rho_inv     = Float64(1.0/ρ[p])

            pixmin = Vector{Int64}(undef,3)
            pixmax = Vector{Int64}(undef,3)

            @inbounds for dim = 1:3

                pixmin[dim] = find_min_pixel(pos[dim], hsml, minCoords[dim],
                                             pixsize_inv)

                pixmax[dim] = find_max_pixel(pos[dim], hsml, minCoords[dim],
                                             pixsize_inv, param.Npixels[dim])

            end


            if conserve_quantities

                # calculate pixel area
                xp1 = (pos[1] - hsml) * pixsize_inv
                xp2 = (pos[1] + hsml) * pixsize_inv
                yp1 = (pos[2] - hsml) * pixsize_inv
                yp2 = (pos[2] + hsml) * pixsize_inv
                zp1 = (pos[3] - hsml) * pixsize_inv
                zp2 = (pos[3] + hsml) * pixsize_inv

                pix_volume = (xp2 - xp1) * (yp2 - yp1) * (zp2 - zp1) * param.pixelArea * param.pixelSideLength

                d4 = (m * rho_inv) / pix_volume

                # number of pixels over which the particle is distributed
                N_distr = Int64((pixmax[1] - pixmin[1] + 1 ) * 
                                (pixmax[2] - pixmin[2] + 1 ) *
                                (pixmax[3] - pixmin[3] + 1 ) )

                # allocate arrays for weights
                kernel_tab = zeros(N_distr)
                d1_tab     = zeros(N_distr)
                d2_tab     = zeros(N_distr)
                d3_tab     = zeros(N_distr)
                dx_tab     = zeros(N_distr)

                # weights table
                wit1    = 0.
                witd    = 0.
                wit2    = 0
                wit2tot = 0

                # pixel count for attributing weights to pixels
                N_count = 1

                    # first loop to calculate weights
                @inbounds for i = pixmin[2]:pixmax[2]

                    dx = pos[1] - ( i - 0.5 ) * param.pixelSideLength
                    dimin = max(xp1, i - 1.0)
                    dimax = min(xp2, i)

                    @inbounds for j = pixmin[2]:pixmax[2]

                        dy = pos[2] - ( j - 0.5 ) * param.pixelSideLength
                        djmin = max(yp1, j - 1.0)
                        djmax = min(yp2, j)

                        @inbounds for k = pixmin[3]:pixmax[3]

                            dz = pos[3] - ( k - 0.5 ) * param.pixelSideLength
                            dkmin = max(xp1, k - 1.0)
                            dkmax = min(xp2, k)

                            d1 = dimax - dimin
                            d2 = djmax - djmin
                            d3 = dkmax - dkmin

                            # compute distance to pixel center in units of hsml
                            dx_tab[N_count] = get_d_hsml_3D(dx, dy, dz, hsml_inv)
                            # update pixel value
                            wi = kernel_value_3D(kernel, dx_tab[N_count], hsml_inv)

                            wit1 += wi * d1 * d2 * d3
                            witd += d1 * d2 * d3

                            if wi > 0
                                wit2 += 1
                            end

                            # store data in arrays
                            d1_tab[N_count]     = d1
                            d2_tab[N_count]     = d2
                            d3_tab[N_count]     = d3
                            kernel_tab[N_count] = wi

                            wit2tot += 1

                            N_count += 1

                        end # end x-loop
                    end # end y-loop
                end # end z-loop

                #fak = 1.0/( (pixmax[2] - pixmin[2]) * (pixmax[1] - pixmin[1]) )
                if wit1 > 0
                    fak = wit2/wit1
                else
                    wit2 = wit2tot
                    fak = wit2tot / witd
                end

                fak_hsml = pix_volume/(wit2 * param.pixelArea * param.pixelSideLength)

                # reset the counter
                N_count = 1

            end # conserve_quantities

            if !conserve_quantities
                bin_prefac = bin_q * m * rho_inv
            end

            # second loop to calculate value
            @inbounds for i = pixmin[1]:pixmax[1]
                if !conserve_quantities
                    dx = pos[1] - ( i - 0.5 ) * param.pixelSideLength
                end
                @inbounds for j = pixmin[2]:pixmax[2]
                    if !conserve_quantities
                        dy = pos[2] - ( j - 0.5 ) * param.pixelSideLength
                    end
                    @inbounds for k = pixmin[3]:pixmax[3]

                        if !conserve_quantities

                            dz = pos[3] - ( k - 0.5 ) * param.pixelSideLength
                            # calculate simple distance to pixel center
                            distance_hsml = get_d_hsml_3D(dx, dy, dz, hsml_inv)

                            # update pixel value
                            image[i, j, k] += bin_prefac * kernel_value_3D(kernel, distance_hsml, hsml_inv)
                        
                        else

                            if wit1 <= 1
                                wi = fak * fak_hsml
                            else
                                wi = kernel_tab[N_count] * fak * fak_hsml
                            end
                            # update pixel value with weights
                            image[i, j, k] += bin_q * wi * d4
                                              d1_tab[N_count] * d2_tab[N_count] * d3_tab[N_count] 
                            
                        end

                    end # end x-loop
                end # end y-loop
            end # end z-loop

        end # end check if in image

        # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx += 1
            ProgressMeter.update!(P, idx)
            #unlock(P_lock)
        end
    end

    if show_progress
        @info "Mapped $particles_in_image / $N particles."
    end

   return image
end

