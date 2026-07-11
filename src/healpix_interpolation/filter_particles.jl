"""
    find_in_shell(Δx, radius_limits)

Checks if a particle is contained in a shell around 
"""
function find_in_shell(Δx::Vector{T}, radius_limits::Vector{T}) where T
    @. ( radius_limits[1] <= Δx <= radius_limits[2] )
end


"""
    filter_sort_particles(Pos, Hsml, Bin_q, Weights, center, radius_limits)

Filters the particles that are in the shell that should be mapped and sorts them by according to their position along the line of sight.
Returns arrays with only relevant particles
"""
function filter_sort_particles(Pos, Hsml, M, Rho, Bin_q, Weights, center, radius_limits, calc_mean)

    # subtract center
    Pos .-= center

    # calculate radii of all particles
    Δx = @. √(Pos[1, :]^2 + Pos[2, :]^2 + Pos[3, :]^2)

    # select contributing particles (boolean mask over the *original* order)
    sel = find_in_shell(Δx, radius_limits)

    if !calc_mean
        # additionally require a non-zero bin quantity.
        # NOTE: this must be an element-wise AND over the full-length mask.
        # The previous `sel = sel[Bin_q[sel] .> 0.0]` indexed the mask with a
        # shorter logical vector, which throws a BoundsError whenever the shell
        # cut removed any particle.
        sel = sel .& (Bin_q .> 0.0)
    end

    # sort particles by radial distance (descending: far -> near)
    sorted = reverse(sortperm(Δx))

    # reorder the mask into sorted order, then keep only the selected
    # particles. This yields the indices of the in-shell particles in
    # far-to-near order.
    # NOTE: `sorted[sel]` (used previously) mixes index spaces -- it selects
    # permutation entries at the *original* positions of in-shell particles,
    # which silently includes out-of-shell particles and drops in-shell ones
    # whenever `radius_limits` is not the default [0, Inf].
    idx = sorted[sel[sorted]]

    # allocate arrays only for relevant particles
    pos = Pos[:, idx]
    hsml = Hsml[idx]
    bin_q = Bin_q[idx]
    weights = Weights[idx]
    rho = Rho[idx]
    m = M[idx]

    # free memory allocated for full arrays
    Pos = nothing
    Hsml = nothing
    Bin_q = nothing
    Weights = nothing
    M = nothing 
    Rho = nothing
    GC.gc()

    # return relevant particle arrays
    return pos, hsml, m, rho, bin_q, weights
end