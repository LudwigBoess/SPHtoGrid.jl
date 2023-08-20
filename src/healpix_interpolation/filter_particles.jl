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

    # select contributing particles
    sel = find_in_shell(Δx, radius_limits)

    if !calc_mean
        sel = sel[Bin_q[sel] .> 0.0]
    end

    # sort particles by radial distance
    sorted = reverse(sortperm(Δx))

    # allocate arrays only for relevant particles
    pos = Pos[:, sorted[sel]]
    hsml = Hsml[sorted[sel]]
    bin_q = Bin_q[sorted[sel]]
    weights = Weights[sorted[sel]]
    rho = Rho[sorted[sel]]
    m = M[sorted[sel]]

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