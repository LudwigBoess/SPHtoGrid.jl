# """
#     rotation_measure(n_cm3::Real, B_par::Real; ν_obs=nothing)

# Computes the rotation measure of the parallel magnetic field along the LOS.
# """
# function rotation_measure(n_cm3::Real, B_par::Real; ν_obs=nothing)
#     RM = 812 * n_cm3 * 1.e3 * B_par * 1.e6
#     if !isnothing(ν_obs)
#         RM *= c_light^2 * 1.e-4 / ν_obs
#     end
#     return mod(RM, π)
# end

"""
    get_dz(M, rho, hsml)

Helper function to get the particle depth along the LOS.
"""
get_dz(M, rho, hsml) = M / rho / (4 * hsml^2)

"""
    rotation_measure(n_cm3::Real, B_los::Real, dz::Real; ν_obs=nothing)

Computes the rotation measure of the parallel magnetic field along the LOS.
"""
function rotation_measure(n_cm3::Real, B_los::Real, dz::Real; ν_obs=nothing)
    RM = faraday_prefac * n_cm3 * B_los
    if !isnothing(ν_obs)
        RM *= c_light^2 * dz / ν_obs
    end
    return mod(RM, π)
end