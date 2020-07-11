module SPHtoGrid

    using Base.Threads
    using Distributed

    include("kernels.jl")
    include("sph_types.jl")
    include("mapping_functions.jl")
    include("smac1_utility.jl")
    include("smac2_utility.jl")


    export Cubic, Quintic, WendlandC4, WendlandC6,    # Kernels
           kernel_value_2D, kernel_value_3D,
           mappingParameters,                         # parameters for SPH mapping
           sphMapping,                                # main function for mapping 
           read_smac1_binary_info,
           read_smac1_binary_image,
           write_smac1_par,
           write_smac2_par


    """ 
    'Domain decomposition':
    Calculate array slices for each worker.
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
        sphMapping(Pos, HSML, M, ρ, Bin_Quant;
                param::mappingParameters,
                kernel::SPHKernel[,
                show_progress::Bool=true,
                conserve_quantities::Bool=true,
                parallel::Bool=true,
                dimensions::Int=2])

    Maps the data in `Bin_Quant` to a grid. Parameters of mapping are supplied in
    `param` and the kernel to be used in `kernel`.

    # Arguments
    - `Pos`: Array with particle positions.
    - `HSML`: Array with particle hsml.
    - `M`: Array with particle masses.
    - `ρ`: Array with particle densities.
    - `Bin_Quant`: Array with particle quantity to be mapped.
    - `kernel::SPHKernel`: Kernel object to be used.
    - `show_progress::Bool=true`: Show progress bar.
    - `parallel::Bool=true`: Run on multiple processors.
    - `conserve_quantities::Bool=true`: If quantities should be conserved while mapping, like in Smac (Dolag et. al. 2005).
    - `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
    """
    function sphMapping(Pos, HSML, M, ρ, Bin_Quant;
                        param::mappingParameters,
                        kernel::SPHKernel,
                        show_progress::Bool=true,
                        conserve_quantities::Bool=false,
                        parallel::Bool=true,
                        dimensions::Int=2)

        if (dimensions == 2)

            if !parallel
                return sphMapping_2D(Pos, HSML, M, ρ, Bin_Quant;
                                    param=param, kernel=kernel,
                                    conserve_quantities=conserve_quantities,
                                    show_progress=show_progress)
            else
                @info "Running on $(nworkers()) cores."

                N = length(M)
                futures = Array{Future}(undef, nworkers())

                # 'Domain decomposition':
                # calculate array slices for each worker
                batch = domain_decomposition(N, nworkers())

                # start remote processes
                for (i, id) in enumerate(workers())
                    futures[i] = @spawnat id sphMapping_2D(Pos[batch[i],:], HSML[batch[i]],
                                                        M[batch[i]], ρ[batch[i]],
                                                        Bin_Quant[batch[i]];
                                                        param=param, kernel=kernel,
                                                        conserve_quantities=conserve_quantities,
                                                        show_progress=false)
                end

                # get and reduce results
                return sum(fetch.(futures))
            end

        elseif (dimensions == 3 )
            if !parallel
                return sphMapping_3D(Pos, HSML, M, ρ, Bin_Quant;
                                    param=param, kernel=kernel,
                                    show_progress=show_progress)
            else
                @info "Running on $(nworkers()) cores."

                N = length(M)
                futures = Array{Future}(undef, nworkers())

                # 'Domain decomposition':
                # calculate array slices for each worker
                batch = domain_decomposition(N, nworkers())

                # start remote processes
                for (i, id) in enumerate(workers())
                    futures[i] = @spawnat id sphMapping_3D(Pos[batch[i],:], HSML[batch[i]],
                                                        M[batch[i]], ρ[batch[i]],
                                                        Bin_Quant[batch[i]];
                                                        param=param, kernel=kernel,
                                                        conserve_quantities=conserve_quantities,
                                                        show_progress=false)
                end

                # get and reduce results
                return sum(fetch.(futures))
            end
        end



    end


end # module
