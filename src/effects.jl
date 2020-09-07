"""
    density_2D( rho::Real, pixelSideLength::Real, 
                Mass::Real=1.989e43, Length::Real=3.085678e21)

Computes the 2D density.

# Returns
2D surface density in untis of [g/cm^2].
"""
function density_2D(rho::Real, pixelSideLength::Real, 
                    Mass::Real=1.989e43, Length::Real=3.085678e21)
    return rho * ( (Mass / Length^2 ) * pixelSideLength )
end


"""
    x_ray_emission(n_cm3::Real, T_eV::Real; 
                   Emin::Real=5.e4, Emax::Real=1.e10, 
                   xH::Real=0.76)

X-Ray emission of a particle with number density `n_cm3` in ``1/cm^3`` and temperature `T` in ``eV``.
`Emin` and `Emax` give the minimum and maximum energy of the oberservation.
`xH` gives the hydrogen fraction used in the simulation.

# Returns
X-Ray surface brightness contribution in units of [erg/cm^2/s/Hz].
"""
function x_ray_emission(n_cm3::Real, T_eV::Real; 
                        Emin::Real=5.e4, Emax::Real=1.e10, 
                        xH::Real=0.76)

    prefac = 4.0 * 2.42e-24 / (1 + xH)
    k_B = 1.38066e-16

    return prefac * sqrt( k_B * T_eV * 1.e-3 ) * n_cm3^2 *
            ( exp( -Emin / (k_B * T_eV) ) - exp( -Emax / (k_B * T_eV) ) )
end


"""
    Tcmb(z::Real)

Computes the temperature of the CMB at redshift `z`.
"""
Tcmb(z::Real) = ( 2.728 * ( 1.0 + z ) )


"""
    kSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

Prefactor for the kinetic Sunyaev-Zel'dovich effect.
"""
function kSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

    c   = 2.9979e10
    σ_T = 6.65245e-25
    
    kSzPrefac  = -1.0 * σ_T / c

    if DI_over_I
        h   = 6.6261e-27
        k_B = 1.38066e-16
        x   = h * ν / ( k_B * Tcmb(z) )

        kSzPrefac *= exp(x) - 1 / (x * exp(x))
    end

    return kSzPrefac
end


"""
    yPrefac()

Prefactor for the Compton-Y parameter.
"""
function yPrefac()   

    k_B = 1.38066e-16
    c   = 2.9979e10
    σ_T = 6.65245e-25
    m_e = 9.10953e-28

    return σ_T * k_B / (m_e * c * c)
end

"""
    kinetic_SZ(n_cm3::Real, vel_y_cgs::Real, 
                    ν::Real=1.e9, z::Real=0.0; 
                    DI_over_I::Bool=false)

Computes the kinetic Sunyaev-Zel'dovich effect from electron density `n_cm3` and velocity in y-direction to the projection plane in cgs units `vel_y_cgs`.
If `DI_over_I` is set to `true` you also need to provide an observation frequency `ν` and redshift `z`.
"""
function kinetic_SZ(n_cm3::Real, vel_y_cgs::Real, 
                    ν::Real=1.e9, z::Real=0.0; 
                    DI_over_I::Bool=false)
    return kSzPrefac(ν, z, DI_over_I) * n_cm3 * vel_y_cgs
end

"""
    comptonY(n_cm3::Real, T_K::Real, z::Real)

Computes the Compton-Y parameter from electron density `n_cm3` and temperature `T` in Kelvin at redshift `z`.
"""
function comptonY(n_cm3::Real, T_K::Real, z::Real)
    return yPrefac() * n_cm3 * ( T_K - Tcmb(z) )
end

"""
    tSzPrefac(ν::Real, z::Real)

Computes the prefactor for the thermal Sunyaev-Zel'dovich effect.
"""
function tSzPrefac(ν::Real, z::Real, DI_over_I::Bool)
    h   = 6.6261e-27
    k_B = 1.38066e-16
    x   = h * ν / ( k_B * Tcmb(z) )
    tSzPrefac = (x * (exp(x) + 1.0) / (exp(x) - 1.0) - 4.0)

    if DI_over_I
        tSzPrefac *= exp(x) - 1 / (x * exp(x))
    end

    return tSzPrefac

end

"""
    thermal_SZ(n_cm3::Real, T::Real, z::Real=0.0, ν::Real=1.44e9)

Computes the thermal Sunyaev-Zel'dovich effect for electron density `n_cm3` and temperature `T` in Kelvin at redshift `z` and observer frequency `ν`.
"""
function thermal_SZ(n_cm3::Real, T::Real, 
                    z::Real=0.0, ν::Real=1.44e9; 
                    DI_over_I::Bool=false)
    return tSzPrefac(ν, z, DI_over_I) * comptonY(n_cm3, T, z)
end