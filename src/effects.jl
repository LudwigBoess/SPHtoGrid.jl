"""
    density_2D( rho::Real, pixelSideLength::Real, 
                Mass::Real=1.989e43, Length::Real=3.085678e21)

Computes the 2D density.
"""
function density_2D(rho::Real, pixelSideLength::Real, 
                    Mass::Real=1.989e43, Length::Real=3.085678e21)
    return rho * ( (Mass / Length^2 ) * pixelSideLength )
end

"""
    Tcmb(z::Real)

Computes the temperature of the CMB at redshift `z`.
"""
Tcmb(z::Real) = ( 2.728 * ( 1.0 + z ) )

k_B = 1.38066e-16
c   = 2.9979e10
σ_T = 6.65245e-25
m_e = 9.10953e-28

"""
    kSzPrefac  = -1.0 * σ_T / c

Prefactor for the kinetic Sunyaev-Zel'dovich effect.
"""
global const kSzPrefac  = -1.0 * σ_T / c

"""
    yPrefac    = σ_T * k_B / (m_e * c * c)

Precator for the Compton-Y parameter.
"""
global const yPrefac    = σ_T * k_B / (m_e * c * c)

"""
    kinetic_SZ(n_cm3::Real, vel_y_cgs::Real)

Computes the kinetic Sunyaev-Zel'dovich effect from electron density `n_cm3` and velocity in y-direction to the projection plane in cgs units `vel_y_cgs`.
"""
function kinetic_SZ(n_cm3::Real, vel_y_cgs::Real)
    return kSzPrefac * n_cm3 * vel_y_cgs
end

"""
    comptonY(n_cm3::Real, T::Real, z::Real)

Computes the Compton-Y parameter from electron density `n_cm3` and temperature `T` at redshift `z`.
"""
function comptonY(n_cm3::Real, T::Real, z::Real)
    return yPrefac * n_cm3 * ( T - Tcmb(z) )
end

"""
    tSzPrefac(ν::Real, z::Real)

Computes the prefactor for the thermal Sunyaev-Zel'dovich effect.
"""
function tSzPrefac(ν::Real, z::Real)
    h   = 6.6261e-27
    k_B = 1.38066e-16
    x   = h * ν / ( k_B * Tcmb(z) )
    return (x * (exp(x) + 1.0) / (exp(x) - 1.0) - 4.0)
end

"""
    thermal_SZ(n_cm3::Real, T::Real, z::Real=0.0, ν::Real=1.44e9)

Computes the thermal Sunyaev-Zel'dovich effect for electron density `n_cm3` and temperature `T` at redshift `z` and observer frequency `ν`.
"""
function thermal_SZ(n_cm3::Real, T::Real, z::Real=0.0, ν::Real=1.44e9)
    return tSzPrefac(ν, z) * comptonY(n_cm3, T, z)
end