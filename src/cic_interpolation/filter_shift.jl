"""
    check_center_and_move_particles(x, par::mappingParameters)

Mapping only works if all coordinates are positive. This function shifts the particles accordingly.
"""
function check_center_and_move_particles(x::Array{T}, par::mappingParameters) where T

    # explicitly copy to its own variable to avoid memory overwrite
    cen  = copy(par.center)
    xlim = copy(par.x_lim)
    ylim = copy(par.y_lim)
    zlim = copy(par.z_lim)
    shift = 0.0
    
    @inbounds for i = 1:size(x,2), dim = 1:3
        x[dim, i] -= cen[dim]
    end

    xlim .-= cen[1]
    ylim .-= cen[2]
    zlim .-= cen[3]

    return x, mappingParameters(center=[0.0, 0.0, 0.0], x_lim=xlim, y_lim=ylim, z_lim=zlim, #pixelSideLength=par.pixelSideLength, 
                                Npixels=maximum(par.Npixels), 
                                boxsize=par.boxsize)
end


"""
    filter_particles_in_image(x, hsml, param::mappingParameters)

Checks if a particle is contained in the image and returns an array of `Bool`.
"""
function filter_particles_in_image(pos::Array{T}, hsml::Array{T}, param::mappingParameters) where T

    N = size(hsml,1)

    p_in_image = falses(N)

    if param.periodic
        k_start = 0
    else
        k_start = 7
    end

    in_image = false

    @inbounds for p = 1:N

        @inbounds for k_periodic = k_start:7

            if param.periodic
                x, y, z = find_position_periodic(pos[:,p], k_periodic, param.boxsize)
            else
                x, y, z = pos[:,p]
            end

            in_image = particle_in_image(x, y, z, hsml[p], param.halfsize)
            
            if in_image
                break
            end

        end # k_periodic

        p_in_image[p] = in_image

    end # p

    return p_in_image
end