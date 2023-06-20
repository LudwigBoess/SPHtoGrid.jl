"""
    rotation_measure(n_cm3::Real, B_par::Real; ν_obs=nothing)

Computes the rotation measure of the parallel magnetic field along the LOS.
"""
function rotation_measure(n_cm3::Real, B_par::Real; ν_obs=nothing)
    RM = 812 * n_cm3 * 1.e3 * B_par * 1.e6
    if !isnothing(ν_obs)
        RM *= c_light^2 * 1.e-4 / ν_obs
    end
    return mod(RM, π)
end