"""
    center_particles(x, par::mappingParameters)

Shifts all particles so that the image is centered on [0, 0, 0].
"""
function center_particles(x::Matrix{T}, par::mappingParameters) where T

    # explicitly copy to its own variable to avoid memory overwrite
    cen  = copy(par.center)
    xlim = copy(par.x_lim)
    ylim = copy(par.y_lim)
    zlim = copy(par.z_lim)
    
    @inbounds for i = 1:size(x,2), dim = 1:3
        x[dim, i] -= cen[dim]

        # do periodic mapping here (minimum-image convention:
        # shift by the *full* boxsize, not half the box)
        if par.periodic
            if abs(x[dim,i]) > par.boxsize/2
                x[dim, i] = x[dim, i] > 0 ? x[dim,i] - par.boxsize : x[dim,i] + par.boxsize
            end
        end
    end

    xlim .-= cen[1]
    ylim .-= cen[2]
    zlim .-= cen[3]

    return x, mappingParameters(center=[0.0, 0.0, 0.0], x_lim=xlim, y_lim=ylim, z_lim=zlim,
                                pixelSideLength=par.pixelSideLength,
                                boxsize=par.boxsize)
end


"""
    filter_particles_in_image(pos, par::mappingParameters, sort_z::Bool=false)

Checks if a particle is contained in the image and returns an array of `Bool`.
A particle is kept if its *center* lies inside `center ± halfsize`.
"""
function filter_particles_in_image(pos::Array{T}, par::mappingParameters, sort_z::Bool=false) where T

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

    if sort_z
        @info "Sorting particles..."
        sorted = reverse(sortperm(pos[3,:]))
        return sorted[p_in_image]
    else
        return p_in_image
    end
end

"""
    filter_particles_in_image(pos, hsml, par::mappingParameters, sort_z::Bool=false)

`hsml`-aware variant: a particle is kept if its smoothing kernel *overlaps*
the image, i.e. its center lies inside `center ± (halfsize + hsml)`.
This prevents edge particles whose kernel reaches into the image (but whose
center is just outside) from being dropped, which would lose flux/mass at
the image boundary (Correctness issue C11 / mass-conservation B2).
"""
function filter_particles_in_image(pos::Array{T}, hsml::AbstractVector,
                                    par::mappingParameters, sort_z::Bool=false) where T

    N = size(pos,2)

    corner_lower_right = par.center - par.halfsize
    corner_upper_right = par.center + par.halfsize

    # allocate array to store if particle is in image
    p_in_image = trues(N)

    @inbounds for i = 1:N
        h = hsml[i]
        in_image = true
        for dim = 1:3
            if !((corner_lower_right[dim] - h) <= pos[dim, i] <= (corner_upper_right[dim] + h))
                in_image = false
            end
        end
        p_in_image[i] = in_image
    end

    if sort_z
        @info "Sorting particles..."
        sorted = reverse(sortperm(pos[3,:]))
        return sorted[p_in_image]
    else
        return p_in_image
    end
end