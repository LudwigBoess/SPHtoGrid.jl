"""
    check_center_and_move_particles(x, par::mappingParameters)

Mapping only works if all coordinates are positive. This function shifts the particles accordingly.
"""
function check_center_and_move_particles(x::Matrix{T}, par::mappingParameters) where T

    # explicitly copy to its own variable to avoid memory overwrite
    cen  = copy(par.center)
    xlim = copy(par.x_lim)
    ylim = copy(par.y_lim)
    zlim = copy(par.z_lim)
    
    @inbounds for i = 1:size(x,2), dim = 1:3
        x[dim, i] -= cen[dim]

        # do periodic mapping here
        if par.periodic
            if abs(x[dim,i]) > par.boxsize/2
                x[dim, i] = x[dim, i] > 0 ? x[dim,i] - par.boxsize/2 : x[dim,i] + par.boxsize/2 
            end
        end
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
function filter_particles_in_image(pos::Array{T}, par::mappingParameters) where T

    N = size(pos,2)

    corner_lower_right = par.center - par.halfsize
    corner_upper_right = par.center + par.halfsize
    
    # allocate array to store if particle is in image
    p_in_image = trues(N)

    @inbounds for i = 1:N
        in_image = true
        for dim = 1:3
            if !(corner_lower_right[dim] <= pos[dim, i] <= corner_upper_right[dim])
                in_image = false
            end
        end
        p_in_image[i] = in_image
    end

    return p_in_image
end