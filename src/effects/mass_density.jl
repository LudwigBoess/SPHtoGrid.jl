"""
    mass_density(pos::Matrix{<:Real}, mass::Vector{<:Real}; 
                 kernel::AbstractSPHKernel=Cubic(Float64, 3),
                 Nneighbors::Integer=32,
                 verbose::Bool=false)

Compute the density of an arbitrary particle distribution using the SPH method. 
The density is computed in input units, so additional unit conversion to cgs units is required, if input units are not cgs.

## Arguments:
- `pos`: particle position in physical code units.
- `mass`: particle mass in physical code units.

## Keyword Arguments
- `kernel::AbstractSPHKernel=Cubic(Float64, 3)`: SPH kernel to use for the density estimate.
- `Nneighbors::Integer=32`: Number of neighboring particles to use for the density estimate.
- `verbose::Bool=false`: If set to true gives progress reports and progress bar.

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function mass_density(pos::Matrix{<:Real}, mass::Vector{<:Real}; 
                      kernel::AbstractSPHKernel=Cubic(Float64, 3),
                      Nneighbors::Integer=32,
                      boxsize::Vector{<:Real}=zeros(3),
                      verbose::Bool=false)

    # initialize arrays for density and smoothing length
    rho = Vector{Float64}(undef, length(mass))
    hsml = Vector{Float64}(undef, length(mass))

    # define metric 
    if boxsize == zeros(3)
        println("Non-periodic box")
        metric = Euclidean()
    else
        println("Periodic box")
        metric = PeriodicEuclidean(boxsize)
    end

    # build tree
    if verbose
        println("Building BallTree")
    end
    tree = BallTree(pos, metric)
    if verbose
        println("BallTree built")
    end

    if verbose
        println("Running SPH density loop on $(nthreads()) threads")
        p = Progress(length(hsml))
    end

    # parallel loop over all particles
    @threads for i ∈ eachindex(mass)

        # get the Nneighbors nearest neighbors and their distances
        idxs, rᵢⱼ = knn(tree, pos[:, i], Nneighbors+1, true)

        # we take the distance to the Nneighbors-th neighbor as the smoothing length
        hsml[i] = maximum(rᵢⱼ)

        # precompute inverse here to save time
        h_inv = 1 / hsml[i]

        # initialize density with zero
        rho[i] = 0.0

        # loop over neighbors to compute density
        @inbounds for j = 1:length(rᵢⱼ)
            # SPH density estimate
            rho[i] += mass[idxs[j]] * 𝒲(kernel, rᵢⱼ[j] * h_inv, h_inv)
        end

        # Corrects the density estimate for the kernel bias
        # See Dehnen&Aly 2012, eq. 18+19
        rho[i] = bias_correction(kernel, rho[i], mass[i], h_inv, Nneighbors)
        
        # update progress
        if verbose
            next!(p)
            flush(stdout)
            flush(stderr)
        end
    end

    # return density and smoothing length
    return rho, hsml
end

