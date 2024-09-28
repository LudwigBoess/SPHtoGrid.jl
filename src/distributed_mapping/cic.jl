# the job channels take their ID as a 32-bit integer
const cic_jobs = RemoteChannel(()->Channel{Int}(32))

# the maps are pointers with 80 bits each
const cic_results = RemoteChannel(()->Channel{Matrix}(160))


"""
    distributed_allsky_map( allsky_filename::String, 
                            Nside::Integer, Nsubfiles::Integer, 
                            mapping_function::Function;
                            reduce_image::Bool=true)

Dynamically dispatches workers to compute one allsky map per subfile, sum up the cic_results and save to a fits file.

# Arguments
- `allsky_filename::String`: Name of the file under which the image should be saved
- `Nside::Integer`: `Nside` for healpix map, must be a power of 2! `Nside = 2^N`.
- `Nsubfiles::Integer`: Number of subfiles the snapshot is distributed over.
- `mapping_function::Function`: The function to be executed per subfile. Must have a call to [`allsky_map`](@ref) as return value.
- `reduce_image`: If the final image should be divided by the weight image set to `true`
"""
function distributed_cic_map(cic_filename::String, Nsubfiles::Integer,
                             mapping_function::Function,
                             Nimages::Integer, param::mappingParameters;
                             reduce_image::Bool=true,
                             Ndim::Integer=2,
                             snap::Integer=0,
                             units::String="")

    println("starting workers")

    for p ∈ workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, cic_jobs, cic_results, mapping_function)
    end
    
    println("allocating image")
    N_distr = 1
    for i = 1:Ndim
        N_distr *= param.Npixels[i]
    end
    sum_map = zeros(Float64, N_distr, Nimages+1)

    # feed the cic_jobs channel with "Nsubfiles" cic_jobs
    errormonitor(@async make_jobs(Nsubfiles)); 
    
    println("running")
    flush(stdout); flush(stderr)

    # loop over all subfiles
    @time while Nsubfiles > 0

        # collect result for subfile
        map = take!(cic_results)

        # sum up contribution
        @inbounds for i ∈ eachindex(allsky_map)
            if !isnan(allsky_map[i]) && !isinf(allsky_map[i])
                sum_map[i]  += map[i]
            end
        end

        map = nothing 
        GC.gc()

        # count down
        Nsubfiles -= 1
    end

    println("reducing image")
    flush(stdout); flush(stderr)

    if Ndim == 2
        image = reduce_image_2D(sum_map, Npixels, Npixels, reduce_image)
    elseif Ndim == 3
        image = reduce_image_3D(sum_map, Npixels, Npixels, Npixels)
    else
        error("Only 2D and 3D are possible")
    end

    println("saving")
    flush(stdout); flush(stderr)
    write_fits_image(cic_filename, image, param, snap = snap, units = units)

    println("Image values")
    println("    Min:    $(minimum(image[.!isnan.(image)]))")
    println("    Max:    $(maximum(image[.!isnan.(image)]))")
    println("    Mean:   $(mean(image[.!isnan.(image)]))")
    println("    Median: $(median(image[.!isnan.(image)]))")
    flush(stdout); flush(stderr)
end
