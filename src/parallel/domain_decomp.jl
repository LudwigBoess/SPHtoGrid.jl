""" 
    domain_decomposition(N::Int64, N_workers::Int64)

Calculate relevant array slices for each worker.
Could be done better!
"""
function domain_decomposition(N::Integer, N_workers::Integer)

    batch = Array{typeof(1:2)}(undef, N_workers)
    size = Int(floor(N/N_workers))
    @inbounds for i = 1:N_workers-1
        batch[i] = 1+(i-1)*size:i*size
    end
    batch[N_workers] = 1 + (N_workers-1)*size:N

    return batch
end