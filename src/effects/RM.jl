"""
    get_dz(M, rho, hsml)

Helper function to get the particle depth along the LOS in the given units.
All additional unit conversion must be performed by hand.
"""
get_dz(M, rho, hsml) = M / rho / (4 * hsml^2)

"""
    rotation_measure(n_cm3::Real, B_los::Real, dz::Real; ν_obs=nothing)

Computes the rotation measure of the parallel magnetic field along the LOS.

## Arguments
- `n_cm3::Real`: Electron number density in [cm^-3].
- `B_los::Real`: Magnetic field strength along the LOS in [G].

## Returns
- `RM::Real`: Rotation measure in [rad/m^2].

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function rotation_measure(n_cm3::Real, B_los::Real)
    # RM in rad/m^2
    return faraday_prefac * n_cm3 * B_los * 100.0^2
end


"""
    rotation_measure(n_cm3::Real, B_los::Real, dz::Real, ν_obs::Real)

Computes the rotation measure of the parallel magnetic field along the LOS at a given frequency.
To be used for continuous rotation of polarized emission along the LOS.

## Arguments
- `n_cm3::Real`: Electron number density in [cm^-3].
- `B_los::Real`: Magnetic field strength along the LOS in [G].
- `dz::Real`: Depth along the LOS in [cm]. See [`get_dz`](@ref) for a convenient helper function.
- `ν_obs::Real`: Observing frequency in [Hz].

## Returns
- `RM::Real`: Rotation measure in [rad/cm^2].

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
- stokes: `true`
- sort_z: `true`
- `RM`: use output of this function.
"""
function rotation_measure(n_cm3::Real, B_los::Real, dz::Real, ν_obs::Real)
    # RM in rad/cm^2
    RM = faraday_prefac * n_cm3 * B_los

    # for Stokes parameters we need rotation at given frequency
    RM *= c_light^2 * dz / ν_obs^2

    # store sign of rotation 
    _sign = sign(RM)

    # reduce to minimum angle
    return _sign * mod(RM, π)
end