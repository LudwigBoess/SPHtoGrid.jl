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