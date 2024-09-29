# the job channels take their ID as a 32-bit integer
const jobs = RemoteChannel(()->Channel{Int}(32))
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