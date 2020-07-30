module SPHtoGrid

    using Distributed
    using SPHKernels

    include("parameters.jl")
    include("filter_shift.jl")
    include("domain_decomp.jl")
    include("reduce_image.jl")
    include("reconstruct_grid.jl")
    include("mapping_functions.jl")
    include("smac1_utility.jl")
    include("smac2_utility.jl")
    include("rotate_particles.jl")


    export mappingParameters,                         # parameters for SPH mapping
           sphMapping,                                # main function for mapping 
           filter_particles_in_image,                 # helper function to preselect particles
           get_map_grid_2D,
           get_map_grid_3D,
           read_smac1_binary_info,
           read_smac1_binary_image,
           write_smac1_par,
           write_smac2_par,
           rotate_3D, rotate_3D!,
           project_along_axis,
           part_weight_one,
           part_weight_physical, 
           part_weight_emission,
           part_weight_XrayBand

    

    """
        function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, 
                            M::Array{<:Real}, ρ::Array{<:Real}, 
                            Bin_Quant::Array{<:Real},
                            Weights::Array{<:Real}=ρ;
                            param::mappingParameters,
                            kernel::SPHKernel [,
                            show_progress::Bool=true,
                            parallel::Bool=false,
                            filter_particles::Bool=true,
                            dimensions::Int=2])

    Maps the data in `Bin_Quant` to a grid. Parameters of mapping are supplied in
    `param` and the kernel to be used in `kernel`.

    # Arguments
    - `Pos`: Array with particle positions.
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
    function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, 
                        M::Array{<:Real}, ρ::Array{<:Real}, 
                        Bin_Quant::Array{<:Real},
                        Weights::Array{<:Real}=ρ;
                        param::mappingParameters,
                        kernel::SPHKernel,
                        show_progress::Bool=true,
                        parallel::Bool=false,
                        filter_particles::Bool=true,
                        dimensions::Int=2)

        
        # store number of input particles
        N_in = length(Bin_Quant)

        # filter particles if they are contained in the image
        if filter_particles
            p_in_image = filter_particles_in_image(Pos, HSML, param)
        else
            p_in_image = trues(N_in)
        end

        # if this is not a float it has units, which need to be stripped
        if !(typeof(Bin_Quant[1,1]) <: AbstractFloat)

            if show_progress
                @info "Stripping units..."
            end

            # allocate reduced arrays
            x       = ustrip(Pos[p_in_image, :])
            hsml    = ustrip(HSML[p_in_image])
            m       = ustrip(M[p_in_image])
            rho     = ustrip(ρ[p_in_image])
            bin_q   = ustrip(Bin_Quant[p_in_image])
            weights = ustrip(Weights[p_in_image])

        else
            # allocate reduced arrays
            x       = Pos[p_in_image, :]
            hsml    = HSML[p_in_image]
            m       = M[p_in_image]
            rho     = ρ[p_in_image]
            bin_q   = Bin_Quant[p_in_image]
            weights = Weights[p_in_image]
        end
        
        N_map = length(m)

        @info "Particles in image: $N_map / $N_in"

        # First check if all particles are centered around 0 and shift them if they are not
        x, par = check_center_and_move_particles(x, param)

        if (dimensions == 2)

            if !parallel
                image, w_image = sphMapping_2D(x, hsml, m, rho, bin_q, weights;
                                    param=par, kernel=kernel,
                                    show_progress=show_progress)

                return reduce_image_2D( image, w_image, 
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
                    futures[i] = @spawnat id sphMapping_2D(x[batch[i],:], hsml[batch[i]],
                                                           m[batch[i]], rho[batch[i]],
                                                           bin_q[batch[i]], weights[batch[i]];
                                                           param=par, kernel=kernel,
                                                           show_progress=false)
                end
                
                image, w_image = reduce_futures(fetch.(futures))

                return reduce_image_2D( image, w_image, 
                         param.Npixels[1], param.Npixels[2] )
            end

        elseif (dimensions == 3 )
            if !parallel
                image, w_image = sphMapping_3D(x, hsml, m, rho, bin_q, weights;
                                    param=par, kernel=kernel,
                                    show_progress=show_progress)
                                    
                return reduce_image_3D( image, w_image, 
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
                    futures[i] = @spawnat id sphMapping_3D(x[batch[i],:], hsml[batch[i]],
                                                        m[batch[i]], rho[batch[i]],
                                                        bin_q[batch[i]], weights[batch[i]];
                                                        param=par, kernel=kernel,
                                                        show_progress=false)
                end

                # get and reduce results
                image, w_image = reduce_futures(fetch.(futures))

                return reduce_image_3D( image, w_image, 
                         param.Npixels[1], param.Npixels[2], param.Npixels[3] )
            end
        end



    end


end # module
