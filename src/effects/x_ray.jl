"""
    x_ray_emission(n_cm3::Real, T_K::Real; 
                   Emin::Real=5.e4, Emax::Real=1.e10, 
                   xH::Real=0.76)

X-Ray emission of a particle with number density `n_cm3` in ``1/cm^3`` and temperature `T` in ``K``.
`Emin` and `Emax` give the minimum and maximum energy of the oberservation.
`xH` gives the hydrogen fraction used in the simulation.

# Returns
X-Ray surface brightness contribution in units of [erg/cm^2/s/Hz].
"""
function x_ray_emission(n_cm3::Real, T_K::Real; 
                        Emin::Real=5.e4, Emax::Real=1.e10, 
                        xH::Real=0.76)

    prefac = 4.0 * 2.42e-24 / (1 + xH)
    T_eV   = T_K * cgs2eV

    return prefac * sqrt( k_B * T_eV * 1.e-3 ) * n_cm3^2 *
            ( exp( -Emin / (k_B * T_eV) ) - exp( -Emax / (k_B * T_eV) ) )
end
