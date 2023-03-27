"""
    get_T_keV(U::Vector{<:Real}, mass::Vector{<:Real}, T_eV::Real)

Helper function to compute energy in keV used for Xray emission.
Takes `U` and `mass` in code units and converts it with `T_eV` as the temperature->eV factor from `GadgetUnits`.
"""
function get_T_keV(U::Vector{<:Real}, mass::Vector{<:Real}, T_eV::Real)
    @. U * T_eV * 1.e-3 * mass
end

"""
    x_ray_emission(T_keV::Vector{<:Real}, 
                    m_cgs::Vector{<:Real}, 
                    rho_cgs::Vector{<:Real}; 
                    E0::Real=0.1, E1::Real=2.4, 
                    xH::Real = 0.76)

X-Ray emissivity for particles with temperature `T_keV` in ``keV``, mass `m_cgs` in ``g`` and density `rho_cgs` in ``g/cm^3``.
`Emin` and `Emax` give the minimum and maximum energy of the oberservation.
`xH` gives the hydrogen fraction used in the simulation.

# Returns
X-Ray emissivity in units of [erg/s].

## Arguments:
- `T_keV`: SPH particle temperature [keV]
- `m_cgs`: SPH particle mass in [g]
- `rho_cgs`: SPH particle density in [g/cm^3]
- `E0`: Minimum photon energy for Xray spectrum [keV]
- `E1`: Maximum photon energy for Xray spectrum [keV]
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings
- weight function: [`part_weight_one`](@ref)
- reduce image: `true`
"""
function x_ray_emission(T_keV::Vector{<:Real}, 
                        m_cgs::Vector{<:Real}, 
                        rho_cgs::Vector{<:Real}; 
                        E0::Real=0.1, E1::Real=2.4, 
                        xH::Real = 0.76)

    mol  = 4 / (5 * xH + 3);
    n2ne = (xH + 0.5 * (1 - xH)) /
           (2xH + 0.75 * (1 - xH));

    cutoff = @. exp(-E0 / T_keV) - exp(-E1 / T_keV)

    """
        Steinmetz & Bartelmann, based on Spizer 1968, gg = 1.0 (!?)
    Beside the fact, that it is not clear which value they used for the Gaunt factor
    it is the best formulation, as composition H/He (fr) and conversion from particle
    number to electron number (n2ne) is explicite formulated.
    """
    Lxbol = @. m_cgs * rho_cgs *
            √(T_keV) * 
            4C_j * gg / (1 + xH) * (n2ne / (mol * m_p))^2

    @. Lxbol * cutoff
end