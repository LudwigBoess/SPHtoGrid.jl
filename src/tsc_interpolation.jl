using TriangularShapedCloudInterpolation


function filter_particles_tsc(pos::Array{<:Real}, par::mappingParameters)
    return findall( (par.x_lim[1] .<= pos[:,1] .< par.x_lim[2]) .&
                    (par.y_lim[1] .<= pos[:,2] .< par.y_lim[2]) .&
                    (par.z_lim[1] .<= pos[:,3] .< par.z_lim[2]) )
end


function reduce_image_2D_tsc(tsc::Array{<:Real}, par::mappingParameters)
    
    image = zeros(par.Npixels[2], par.Npixels[1])                    
    @inbounds for ix = 1:par.Npixels[1], iy = 1:par.Npixels[2], iz = 1:par.Npixels[3]

        image[iy, ix] += tsc[ix, iy, iz]

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
- `filter_particles::Bool=true`: If particles should be filtered before mapping.
- `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
"""
function sphMapping(Pos::Array{<:Real}, Bin_Quant::Array{<:Real};
                    param::mappingParameters,
                    show_progress::Bool=true,
                    filter_particles::Bool=true,
                    parallel::Bool=false,
                    dimensions::Int=2)

    
    # store number of input particles
    N_in = length(Bin_Quant)

    if show_progress
        @info "Filtering particles..."
        t1 = time_ns()
    end

    if filter_particles
        # filter particles if they are contained in the image
        p_in_image = filter_particles_tsc(Pos, param)
    end

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

        if filter_particles
            # allocate reduced arrays
            x       = ustrip(Pos[p_in_image, :])
            bin_q   = ustrip(Bin_Quant[p_in_image])
        else
            x       = ustrip(Pos)
            bin_q   = ustrip(Bin_Quant)
        end

    else
        if show_progress
            @info "Assigning arrays..."
            t1 = time_ns()
        end
        
        if filter_particles
            # allocate reduced arrays
            x       = Pos[p_in_image, :]
            bin_q   = Bin_Quant[p_in_image]
        else
            x       = Pos
            bin_q   = Bin_Quant
        end
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
    pos_tsc = get_tsc_positions(x, param.Npixels)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end
    

    if show_progress
        @info "Running TSC interpolation..."
        t1 = time_ns()
    end

    if !parallel
        tsc = TSCInterpolation( bin_q, 
                                pos_tsc, param.Npixels, 
                                average=true)
    else
        if show_progress
            @info "Running on $(nworkers()) cores."
        end

        # Number of particles
        N = length(bin_q)

        # allocate an array of Future objects
        futures = Array{Future}(undef, nworkers())

        # 'Domain decomposition':
        # calculate array slices for each worker
        batch = domain_decomposition(N, nworkers())

        # start remote processes
        for (i, id) in enumerate(workers())
            futures[i] = @spawnat id TSCInterpolation(bin_q[batch[i],:], pos_tsc[batch[i]],
                                                      param.Npixels, 
                                                      average=true)
        end
        
        tsc = sum(fetch.(futures))
    end


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
        #image = copy(transpose(image))
    elseif (dimensions == 3 )
        image = reduce_image_3D_tsc(tsc, param)
    end

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    return image
end