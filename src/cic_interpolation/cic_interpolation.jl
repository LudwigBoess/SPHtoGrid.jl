"""
    sphMapping( Pos::Array{<:Real}, HSML::Array{<:Real}, M::Array{<:Real}, 
                Rho::Array{<:Real}, Bin_Quant::Array{<:Real}, 
                Weights::Array{<:Real}=Rho;
                param::mappingParameters,
                kernel::AbstractSPHKernel,
                show_progress::Bool=true,
                parallel::Bool=false,
                reduce_image::Bool=true,
                return_both_maps::Bool=false,
                filter_particles::Bool=true,
                dimensions::Int=2,
                calc_mean::Bool=false)

Maps the data in `Bin_Quant` to a grid. Parameters of mapping are supplied in
`param` and the kernel to be used in `kernel`.

# Arguments
- `Pos`: Matrix (3xNpart) with particle positions.
- `HSML`: Array with particle hsml.
- `M`: Array with particle masses.
- `Rho`: Array with particle densities.
- `Bin_Quant`: Array with particle quantity to be mapped.
- `Weights`: Array with weights. Defaults to density-weighted.
- `kernel::AbstractSPHKernel`: Kernel object to be used.
- `show_progress::Bool=true`: Show progress bar.
- `parallel::Bool=true`: Run on multiple processors.
- `reduce_image::Bool=true`: If weights need to be applied or not. Set to `false` for [`part_weight_physical`](@ref).
- `return_both_maps::Bool=true`: Returns the full image array. To be used with parallel mapping of subfiles.
- `filter_particles::Bool=true`: Find the particles that are actually contained in the image.
- `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
- `calc_mean::Bool=false`: Calculates the mean value along the line of sight. If set to `false` the particle only contributes if its `Bin_Quant` is larger than 0.
"""
function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, M::Array{<:Real}, 
                    Rho::Array{<:Real}, Bin_Quant::Array{<:Real}, 
                    Weights::Array{<:Real}=Rho;
                    param::mappingParameters,
                    kernel::AbstractSPHKernel,
                    show_progress::Bool=true,
                    parallel::Bool=false,
                    reduce_image::Bool=true,
                    return_both_maps::Bool=false,
                    dimensions::Int=2,
                    calc_mean::Bool=false,
                    sort_z::Bool=false)

    
    # store number of input particles
    N_in = size(Bin_Quant,1)

    # First check if all particles are centered around 0 and shift them if they are not
    if show_progress
        @info "Centering on [0.0, 0.0, 0.0]"
        t1 = time_ns()
    end

    Pos, par = center_particles(Pos, param)

    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    # filter particles if they are contained in the image
    if show_progress
        @info "Filtering particles..."
        t1 = time_ns()
    end

    p_in_image = filter_particles_in_image(Pos, par)


    if show_progress
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    # if this is not a float it has units, which need to be stripped
    if !(typeof(Bin_Quant[1]) <: AbstractFloat)

        if show_progress
            @info "Stripping units..."
            t1 = time_ns()
        end

        # allocate reduced arrays
        x       = ustrip.(Pos[:,p_in_image])
        hsml    = ustrip.(HSML[p_in_image])
        m       = ustrip.(M[p_in_image])
        rho     = ustrip.(Rho[p_in_image])
        bin_q   = ustrip.(Bin_Quant[p_in_image])
        weights = ustrip.(Weights[p_in_image])

    else
        if show_progress
            @info "Assigning arrays..."
            t1 = time_ns()
        end
        
        # allocate reduced arrays
        x       = Pos[:,p_in_image]
        hsml    = HSML[p_in_image]
        m       = M[p_in_image]
        rho     = Rho[p_in_image]
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

            image = cic_mapping_2D(x, hsml, m, rho, bin_q, weights;
                                param=par, kernel=kernel,
                                show_progress=show_progress,
                                calc_mean=calc_mean)

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if return_both_maps
                return image
            end

            return reduce_image_2D( image,
                            param.Npixels[1], param.Npixels[2],
                            reduce_image)
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
                futures[i] = @spawnat id cic_mapping_2D(x[:,batch[i]], hsml[batch[i]],
                                                        m[batch[i]], rho[batch[i]],
                                                        bin_q[batch[i]], weights[batch[i]];
                                                        param=par, kernel=kernel,
                                                        show_progress=false,
                                                        calc_mean=calc_mean)
            end
            
            image = sum(fetch.(futures))

            if show_progress
                t2 = time_ns()
                @info "  elapsed: $(output_time(t1,t2)) s"
            end

            if return_both_maps
                return image
            end

            return reduce_image_2D( image, 
                        param.Npixels[1], param.Npixels[2], 
                        reduce_image)
        end

    elseif (dimensions == 3 )
        if !parallel
            image = cic_mapping_3D(x, hsml, m, rho, bin_q, weights;
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
                futures[i] = @spawnat id cic_mapping_3D(x[:,batch[i]], hsml[batch[i]],
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


using Statistics
using Printf

"""
    map_it(pos_in, hsml, mass, rho, bin_q, weights, k, 
           snap, units, image_path, reduce_image, param, 
           filter_particles=true, parallel=true)

Small helper function to copy positions, map particles and save the fits file.
"""
function map_it(pos_in, hsml, mass, rho, bin_q, weights;
                param::mappingParameters,
                kernel::AbstractSPHKernel=WendlandC6(), 
                snap::Integer=0, 
                units::AbstractString="", 
                image_prefix::String="dummy",
                reduce_image::Bool=true, 
                parallel=true,
                calc_mean::Bool=true, show_progress::Bool=true,
                sort_z::Bool=false,
                projection="xy")

    # copy the positions to new array to be able to shift particles 
    # and also continue mapping with the same array
    pos = copy(pos_in)

    if projection == "xy"
        pos = pos
    elseif projection == "xz"
        pos = rotate_to_xz_plane!(pos)
    elseif projection == "yz"
        pos = rotate_to_yz_plane!(pos)
    elseif typeof(projection) <: AbstractVector
        pos = rotate_3D(pos, projection...)
        projection = "alpha=$(@sprintf("%0.2f", projection[1]))beta=$(@sprintf("%0.2f", projection[2]))gamma=$(@sprintf("%0.2f", projection[3]))"
    else
        error("projection must be either along in 'xy', 'xz', or 'yz' plane of defined by a vector of Euler angles!")
    end 

    # get cic map
    quantitiy_map = sphMapping( pos, hsml, mass, rho,
                                bin_q, weights; show_progress,
                                param, kernel, parallel,
                                reduce_image,
                                calc_mean,
                                sort_z)


    fo_image = image_prefix * ".fits"

    @info "Map Properties:"
    @info "  Max:  $(maximum(quantitiy_map))"
    @info "  Min:  $(minimum(quantitiy_map))"
    @info "  Mean: $(mean(quantitiy_map))"
    @info "  Sum:  $(sum(quantitiy_map))"
    println()
    Npixels = size(quantitiy_map,1) * size(quantitiy_map,2)
    @info "  Nr. of pixels:     $Npixels"
    @info "  Nr. of 0 pixels:   $(length(findall(iszero.(quantitiy_map))))"
    @info "  Nr. of NaN pixels: $(length(findall(isnan.(quantitiy_map))))"
    @info "  Nr. of Inf pixels: $(length(findall(isinf.(quantitiy_map))))"

    write_fits_image(fo_image, quantitiy_map, param, snap = snap, units = units)
end