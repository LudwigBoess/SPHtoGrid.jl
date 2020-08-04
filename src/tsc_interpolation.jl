using TriangularShapedCloudInterpolation


function filter_particles_tsc(pos::Array{<:Real}, par::mappingParameters)
    return findall( (par.x_lim[1] .<= pos[:,1] .< par.x_lim[2]) .&
                    (par.y_lim[1] .<= pos[:,2] .< par.y_lim[2]) .&
                    (par.z_lim[1] .<= pos[:,3] .< par.z_lim[2]) )
end

function construct_positions_tsc(pos::Array{<:Real}, par::mappingParameters)

    pos_tsc = zeros(length(pos[:,1]),3)

    dx   = -(par.x_lim[1] - par.x_lim[2]) / par.Npixels[1]
    pos_tsc[:,1] = ( pos[:,1] .- par.x_lim[1] ) ./ dx 

    dy   = -(par.y_lim[1] - par.y_lim[2]) / par.Npixels[2]
    pos_tsc[:,2] = ( pos[:,2] .- par.y_lim[1] ) ./ dy 

    dz   = -(par.z_lim[1] - par.z_lim[2]) / par.Npixels[3]
    pos_tsc[:,3] = ( pos[:,3] .- par.z_lim[1] ) ./ dz

    return pos_tsc
end

function reduce_image_2D_tsc(tsc::Array{<:Real}, par::mappingParameters)
    
    image = zeros(par.Npixels[1], par.Npixels[2])                    
    @inbounds for ix = 1:par.Npixels[1], iy = 1:par.Npixels[2], iz = 1:par.Npixels[3]

        image[ix, iy] += tsc[ix, iy, iz]

    end # ix, iy, iz

    return image
end

function reduce_image_3D_tsc(tsc::Array{<:Real}, par::mappingParameters)
    
    return tsc
end

"""
    sphMapping( Pos::Array{<:Real}, Bin_Quant::Array{<:Real};
                param::mappingParameters,
                show_progress::Bool=true,
                dimensions::Int=2)

Maps the data in `Bin_Quant` to a grid using triangualar shaped cloud (TSC) interpolation.

# Arguments
- `Pos`: Array with particle positions.
- `Bin_Quant`: Array with particle quantity to be mapped.
- `show_progress::Bool=true`: Show progress bar.
- `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
"""
function sphMapping(Pos::Array{<:Real}, Bin_Quant::Array{<:Real};
                    param::mappingParameters,
                    show_progress::Bool=true,
                    dimensions::Int=2)

    
    # store number of input particles
    N_in = length(Bin_Quant)

    if show_progress
        @info "Filtering particles..."
        t1 = time_ns()
    end

    # filter particles if they are contained in the image
    p_in_image = filter_particles_tsc(Pos, param)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    # if this is not a float it has units, which need to be stripped
    if !(typeof(Bin_Quant[1,1]) <: AbstractFloat)

        if show_progress
            @info "Stripping units..."
            t1 = time_ns()
        end

        # allocate reduced arrays
        x       = ustrip(Pos[p_in_image, :])
        bin_q   = ustrip(Bin_Quant[p_in_image])

    else
        if show_progress
            @info "Assigning arrays..."
            t1 = time_ns()
        end
        
        # allocate reduced arrays
        x       = Pos[p_in_image, :]
        bin_q   = Bin_Quant[p_in_image]
    end

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end
    
    N_map = length(bin_q)

    @info "Particles in image: $N_map / $N_in"

    # First check if all particles are centered around 0 and shift them if they are not
    if show_progress
        @info "constructing grid positions..."
        t1 = time_ns()
    end
    pos_tsc = construct_positions_tsc(x, param)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end
    

    if show_progress
        @info "Running TSC interpolation..."
        t1 = time_ns()
    end

    tsc = TSCInterpolation( bin_q, 
                            pos_tsc[:,1], param.Npixels[1], 
                            pos_tsc[:,2], param.Npixels[2], 
                            pos_tsc[:,3], param.Npixels[3], 
                            average=true)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    if show_progress
        @info "Constructing image..."
        t1 = time_ns()
    end

    if (dimensions == 2)
        image = reduce_image_2D_tsc(tsc, param)
    elseif (dimensions == 3 )
        image = reduce_image_3D_tsc(tsc, param)
    end

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    return image
end