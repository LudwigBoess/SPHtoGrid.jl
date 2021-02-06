"""
    function sphMapping(Pos, HSML, M, ρ, Bin_Quant,
                        Weights=ρ;
                        param::mappingParameters,
                        kernel::SPHKernel [,
                        show_progress::Bool=true,
                        parallel::Bool=false,
                        filter_particles::Bool=true,
                        dimensions::Int=2])

Maps the data in `Bin_Quant` to a grid. Parameters of mapping are supplied in
`param` and the kernel to be used in `kernel`.

# Arguments
- `Pos`: Matrix (3xNpart) with particle positions.
- `HSML`: Array with particle hsml.
- `M`: Array with particle masses.
- `ρ`: Array with particle densities.
- `Bin_Quant`: Array with particle quantity to be mapped.
- `Weights`: Array with weights. Defaults to density-weighted.
- `kernel::SPHKernel`: Kernel object to be used.
- `show_progress::Bool=true`: Show progress bar.
- `parallel::Bool=true`: Run on multiple processors.
- `filter_particles::Bool=true`: Find the particles that are actually contained in the image.
- `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
"""
function sphMapping(Pos, HSML, M, ρ, Bin_Quant, 
                    Weights=ρ;
                    param::mappingParameters,
                    kernel::SPHKernel,
                    show_progress::Bool=true,
                    parallel::Bool=false,
                    filter_particles::Bool=true,
                    dimensions::Int=2)

    
    # store number of input particles
    N_in = size(Bin_Quant,1)

    # check if weights need to be applied
    if Weights == part_weight_one(N_in) || Weights == part_weight_physical(N_in, param)
        reduce_image = false
    else
        reduce_image = true
    end

    # First check if all particles are centered around 0 and shift them if they are not
    if show_progress
        @info "Centering on [0.0, 0.0, 0.0]"
        t1 = time_ns()
    end
    Pos, par = check_center_and_move_particles(Pos, param)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    # filter particles if they are contained in the image
    if show_progress
        @info "Filtering particles..."
        t1 = time_ns()
    end

    if filter_particles
        p_in_image = filter_particles_in_image(Pos, HSML, param)
    else
        p_in_image = trues(N_in)
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

        # allocate reduced arrays
        x       = ustrip(Pos[:,p_in_image])
        hsml    = ustrip(HSML[p_in_image])
        m       = ustrip(M[p_in_image])
        rho     = ustrip(ρ[p_in_image])
        bin_q   = ustrip(Bin_Quant[p_in_image])
        weights = ustrip(Weights[p_in_image])

    else
        if show_progress
            @info "Assigning arrays..."
            t1 = time_ns()
        end
        
        # allocate reduced arrays
        x       = Pos[:,p_in_image]
        hsml    = HSML[p_in_image]
        m       = M[p_in_image]
        rho     = ρ[p_in_image]
        bin_q   = Bin_Quant[p_in_image]
        weights = Weights[p_in_image]
    end

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end
    
    N_map = length(m)

    @info "Particles in image: $N_map / $N_in"
    

    if show_progress
        @info "Mapping..."
        t1 = time_ns()
    end

    if (dimensions == 2)

        if !parallel

            image = sphMapping_2D(x, hsml, m, rho, bin_q, weights;
                                param=par, kernel=kernel,
                                show_progress=show_progress)

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if !reduce_image 
                image[:,2] .= 1.0
            end

            return reduce_image_2D( image,
                            param.Npixels[1], param.Npixels[2] )
        else
            @info "Running on $(nworkers()) cores."

            # Number of particles
            N = length(m)

            # allocate an array of Future objects
            futures = Array{Future}(undef, nworkers())

            # 'Domain decomposition':
            # calculate array slices for each worker
            batch = domain_decomposition(N, nworkers())

            # start remote processes
            for (i, id) in enumerate(workers())
                futures[i] = @spawnat id sphMapping_2D(x[:,batch[i]], hsml[batch[i]],
                                                        m[batch[i]], rho[batch[i]],
                                                        bin_q[batch[i]], weights[batch[i]];
                                                        param=par, kernel=kernel,
                                                        show_progress=false)
            end
            
            image = sum(fetch.(futures))

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if !reduce_image 
                image[:,2] .= 1.0
            end

            return reduce_image_2D( image, 
                        param.Npixels[1], param.Npixels[2] )
        end

    elseif (dimensions == 3 )
        if !parallel
            image = sphMapping_3D(x, hsml, m, rho, bin_q, weights;
                                param=par, kernel=kernel,
                                show_progress=show_progress)

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if !reduce_image 
                image[:,2] .= 1.0
            end
                                
            return reduce_image_3D( image, 
                        param.Npixels[1], param.Npixels[2], param.Npixels[3] )
        else
            @info "Running on $(nworkers()) cores."

            N = length(m)
            futures = Array{Future}(undef, nworkers())

            # 'Domain decomposition':
            # calculate array slices for each worker
            batch = domain_decomposition(N, nworkers())

            # start remote processes
            for (i, id) in enumerate(workers())
                futures[i] = @spawnat id sphMapping_3D(x[:,batch[i]], hsml[batch[i]],
                                                    m[batch[i]], rho[batch[i]],
                                                    bin_q[batch[i]], weights[batch[i]];
                                                    param=par, kernel=kernel,
                                                    show_progress=false)
            end

            # get and reduce results
            image = sum(fetch.(futures))

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if !reduce_image 
                image[:,2] .= 1.0
            end

            return reduce_image_3D( image,
                        param.Npixels[1], param.Npixels[2], param.Npixels[3] )
        end
    end



end