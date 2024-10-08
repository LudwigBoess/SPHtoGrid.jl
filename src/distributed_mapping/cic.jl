"""
    distributed_cic_map(cic_filename::String, Nsubfiles::Integer,
                        mapping_function::Function,
                        Nimages::Integer, param::mappingParameters;
                        reduce_image::Bool=true,
                        Ndim::Integer=2,
                        snap::Integer=0,
                        units::String="",
                        vtk::Bool=true)

Dynamically dispatches workers to compute one cic map per subfile, sum up the results and save to a fits file.

# Arguments
- `cic_filename::String`: Name of the file under which the image should be saved
- `Nsubfiles::Integer`: Number of subfiles the snapshot is distributed over.
- `mapping_function::Function`: The function to be executed per subfile. Must return a `Tuple` of the cic image(s) and the weight_image.
- `Nimages::Integer=1`: Number of cic image that are constructed simultaniously. Only useful if you map multiple quantities at once, so it defaults to `1`.
- `param::mappingParameters`: Parameters for the mapping, needed to save the image in the usual way.
- `reduce_image`: If the final image should be divided by the weight image set to `true`.
- `Ndim::Integer=2`: Dimensions of image, must be either `2` or `3`.
- `snap::Integer=0`: Number of the input snapshot, used for saving the image.
- `units::String=""`: Units of the image, used for saving the image.
"""
function distributed_cic_map(cic_filename::String, Nsubfiles::Integer,
                             mapping_function::Function,
                             param::mappingParameters,
                             Nimages::Integer=1;
                             reduce_image::Bool=true,
                             Ndim::Integer=2,
                             snap::Integer=0,
                             units::String="",
                             vtk::Bool=true)

    println("starting workers")

    for p ∈ workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results, mapping_function)
    end
    
    println("allocating image")
    N_distr = 1
    for i = 1:Ndim
        N_distr *= param.Npixels[i]
    end
    sum_map = zeros(Float64, N_distr, Nimages)
    sum_weights = zeros(Float64, N_distr)

    # feed the jobs channel with "Nsubfiles" jobs
    errormonitor(@async make_jobs(Nsubfiles)); 
    
    println("running")
    flush(stdout); flush(stderr)

    # loop over all subfiles
    @time while Nsubfiles > 0

        # collect result for subfile
        local_map, local_weight_map = take!(results)

        if !isnothing(local_map)
            # sum up contribution
            @inbounds for i ∈ eachindex(local_map)
                if !isnan(local_map[i]) && !isinf(local_map[i])
                    sum_map[i]  += local_map[i]
                end

                if !isnan(local_weight_map[i]) && !isinf(local_weight_map[i])
                    sum_weights[i] += local_weight_map[i]
                end
            end
        end

        local_map = local_weight_map = nothing 
        GC.gc()

        # count down
        Nsubfiles -= 1
    end

    println("reducing image")
    flush(stdout); flush(stderr)

    sum_map = [sum_map sum_weights]

    if Ndim == 2
        image = reduce_image_2D(sum_map, param.Npixels[1], param.Npixels[2], reduce_image)
    elseif Ndim == 3
        image = reduce_image_3D(sum_map, param.Npixels[1], param.Npixels[2], param.Npixels[3])
    else
        error("Only 2D and 3D are possible")
    end

    println("saving")
    flush(stdout); flush(stderr)

    if Ndim == 3 && vtk
        # save 3D grid in vtk file format
        write_vtk_image(cic_filename, image, "map", param, snap = snap, units = units)
    else
        # defaults to fits file
        write_fits_image(cic_filename, image, param, snap = snap, units = units)
    end

    println("Image values")
    println("    Min:    $(minimum(image[.!isnan.(image)]))")
    println("    Max:    $(maximum(image[.!isnan.(image)]))")
    println("    Mean:   $(mean(image[.!isnan.(image)]))")
    println("    Median: $(median(image[.!isnan.(image)]))")
    flush(stdout); flush(stderr)
end
