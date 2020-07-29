module SPHtoGrid

    using Distributed
    using SPHKernels

    include("parameters.jl")
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
           WendlandC6



    """ 
        domain_decomposition(N::Int64, N_workers::Int64)

    Calculate relevant array slices for each worker.
    Could be done better!
    """
    function domain_decomposition(N::Int64, N_workers::Int64)

        batch = Array{typeof(1:2)}(undef, N_workers)
        size = Int(floor(N/N_workers))
        @inbounds for i = 1:N_workers-1
            batch[i] = 1+(i-1)*size:i*size
        end
        batch[N_workers] = 1 + (N_workers-1)*size:N

        return batch
    end

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

        return x, mappingParameters(center=cen, x_lim=xlim, y_lim=ylim, z_lim=zlim, Npixels=maximum(par.Npixels))
    end

    """
        filter_particles_in_image(x, hsml, param::mappingParameters)

    Checks if a particle is contained in the image and returns an array of Bool.
    """
    function filter_particles_in_image(x::Array{<:Real}, hsml::Array{<:Real}, param::mappingParameters)

        p_in_image = falses(length(hsml))
        in_image   = false

        minCoords  = [param.x_lim[1], param.y_lim[1], param.z_lim[1]]
        maxCoords  = [param.x_lim[2], param.y_lim[2], param.z_lim[2]]

        @inbounds for p = 1:length(hsml)

            @inbounds for dim = 1:3

                in_image = check_in_image(Float64(x[p, dim]), Float64(hsml[p]),
                                        minCoords[dim], maxCoords[dim])

                # exit the loop if the particle is not in the image frame
                if !in_image
                    break
                end
            end

            p_in_image[p] = in_image

        end

        return p_in_image
    end

    function get_map_grid_2D(par::mappingParameters)

        x_grid = zeros(par.Npixels[1])
        y_grid = zeros(par.Npixels[2])

        for i = 1:par.Npixels[1]
            x_grid[i] = par.x_lim[1] + ( i - 0.5 ) * par.pixelSideLength
        end

        for i = 1:par.Npixels[2]
            y_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
        end

        return x_grid, y_grid
    end

    function get_map_grid_3D(par::mappingParameters)

        x_grid = zeros(par.Npixels[1])
        y_grid = zeros(par.Npixels[2])
        z_grid = zeros(par.Npixels[3])

        for i = 1:par.Npixels[1]
            x_grid[i] = par.x_lim[1] + ( i - 0.5 ) * par.pixelSideLength
        end

        for i = 1:par.Npixels[2]
            y_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
        end

        for i = 1:par.Npixels[3]
            z_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
        end

        return x_grid, y_grid, z_grid
    end

    """
        sphMapping(Pos, HSML, M, ρ, Bin_Quant, [Weights];
                   param::mappingParameters,
                   kernel::SPHKernel[,
                   show_progress::Bool=true,
                   parallel::Bool=true,
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
                d = sphMapping_2D(x, hsml, m, rho, bin_q, weights;
                                    param=par, kernel=kernel,
                                    show_progress=show_progress)

                return d
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

                # get and reduce results
                return sum(fetch.(futures))
            end

        elseif (dimensions == 3 )
            if !parallel
                d = sphMapping_3D(x, hsml, m, rho, bin_q, weights;
                                    param=par, kernel=kernel,
                                    show_progress=show_progress)
                return d
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
                                                        conserve_quantities=conserve_quantities,
                                                        show_progress=false)
                end

                # get and reduce results
                d = sum(fetch.(futures))
                par = nothing
                return d
            end
        end



    end


end # module
