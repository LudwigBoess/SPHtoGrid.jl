using Statistics
using Distributed

# the job channels take their ID as a 32-bit integer
const jobs = RemoteChannel(()->Channel{Int}(32))

# the maps are pointers with 80 bits each
const results = RemoteChannel(()->Channel{Tuple}(160))


"""
    do_work(jobs, results, f)

Defines the work function `f` on all processes. 
"""
function do_work(jobs, results, f) 
    while true
        job_id = take!(jobs)
        put!(results, f(job_id-1))
    end
end

"""
    make_jobs(n)

Starts up `n` number of jobs.
"""
function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end

"""
    distributed_allsky_map( allsky_filename::String, 
                            Nside::Integer, Nsubfiles::Integer, 
                            mapping_function::Function)

Dynamically dispatches workers to compute one allsky map per subfile, sum up the results and save to a fits file.

# Arguments
- `allsky_filename::String`: Name of the file under which the image should be saved
- `Nside::Integer`: `Nside` for healpix map, must be a multiple of 2! `Nside = 2^N`.
- `Nsubfiles::Integer`: Number of subfiles the snapshot is distributed over.
- `mapping_function::Function`: The function to be executed per subfile. Must have a call to [`allsky_map`](@ref) as return value.
"""
function distributed_allsky_map(allsky_filename::String, 
                                Nside::Integer, Nsubfiles::Integer, 
                                mapping_function::Function)

    println("starting workers")

    for p ∈ workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results, mapping_function)
    end
    
    println("allocating images")
    sum_allsky  = HealpixMap{Float64,RingOrder}(Nside)
    sum_weights = HealpixMap{Float64,RingOrder}(Nside)

    # feed the jobs channel with "Nsubfiles" jobs
    errormonitor(@async make_jobs(Nsubfiles)); 
    
    println("running")
    flush(stdout); flush(stderr)

    # loop over all subfiles
    @time while Nsubfiles > 0

        # collect result for subfile
        allsky_map, weight_map = take!(results)

        # sum up contribution
        @inbounds for i ∈ eachindex(allsky_map)
            if !isnan(allsky_map[i]) && !isinf(allsky_map[i])
                sum_allsky[i]  += allsky_map[i]
            end

            if !isnan(weight_map[i]) && !isinf(weight_map[i])
                sum_weights[i] += weight_map[i]
            end
        end

        allsky_map = nothing 
        weight_map = nothing
        GC.gc()

        # count down
        Nsubfiles -= 1
    end

    println("reducing image")
    flush(stdout); flush(stderr)
    @inbounds for i ∈ eachindex(sum_allsky)
        if !isnan(sum_weights[i]) && !iszero(sum_weights[i]) && !isinf(sum_weights[i])
            sum_allsky[i]  /= sum_weights[i]
        end
    end

    println("saving")
    flush(stdout); flush(stderr)
    if isfile(allsky_filename)
        rm(allsky_filename)
    end
    saveToFITS(sum_allsky, allsky_filename)

    println("Image values")
    println("    Min:    $(minimum(sum_allsky[.!isnan.(sum_allsky)]))")
    println("    Max:    $(maximum(sum_allsky[.!isnan.(sum_allsky)]))")
    println("    Mean:   $(mean(sum_allsky[.!isnan.(sum_allsky)]))")
    println("    Median: $(median(sum_allsky[.!isnan.(sum_allsky)]))")

    flush(stdout); flush(stderr)
end
