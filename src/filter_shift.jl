"""
    check_center_and_move_particles(x, par::mappingParameters)

Mapping only works if all coordinates are positive. This function shifts the particles accordingly.
"""
function check_center_and_move_particles(x::Array{<:Real}, par::mappingParameters)

    # explicitly copy to its own variable to avoid memory overwrite
    cen  = copy(par.center)
    xlim = copy(par.x_lim)
    ylim = copy(par.y_lim)
    zlim = copy(par.z_lim)
    shift = 0.0
    
    @inbounds for i = 1:length(x[:,1]), dim = 1:3
        x[i, dim] -= cen[dim]
    end

    xlim .-= cen[1]
    ylim .-= cen[2]
    zlim .-= cen[3]

    return x, mappingParameters(center=cen, x_lim=xlim, y_lim=ylim, z_lim=zlim, Npixels=maximum(par.Npixels), boxsize=par.boxsize)
end


"""
    filter_particles_in_image(x, hsml, param::mappingParameters)

Checks if a particle is contained in the image and returns an array of `Bool`.
"""
function filter_particles_in_image(pos::Array{<:Real}, hsml::Array{<:Real}, param::mappingParameters)

    N = length(hsml)

    p_in_image = falses(N)

    if param.periodic
        k_start = 0
    else
        k_start = 8
    end

    in_image   = false

    @inbounds for p = 1:N

        for k_periodic = k_start:8

            x, y, z = find_position_periodic(pos, k_periodic, param.boxsize)

            in_image = check_in_image(x, y, z, hsml[p], param.halfsize[1], param.halfsize[2], param.halfsize[3])
            
            if in_image
                break
            end

        end # k_periodic

        p_in_image[p] = in_image

    end # p

    return p_in_image
end