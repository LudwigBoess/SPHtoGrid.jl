using SpecialFunctions
using QuadGK

"""
    α_γ(slope::T) where T

Slope of γ-ray spectrum as a function of the slope of a proton spectrum in energy space. 
"""
α_γ(slope::T) where {T} = 4 / 3 * (slope - 1 / 2)

"""
    δ_γ(α::T) where T

Shape parameter
"""
δ_γ(α::T) where {T} = 0.14 * α^(-1.6) + 0.44


"""
    σ_pp(α::T) where T

Scattering cross-section of proton-proton scatterin in ``cm^2``
"""
σ_pp(α::T) where {T} = 32 * (0.96 + exp(4.4 - 2.4α)) * 1.e-27

"""
    E_γ_powδ(Eγ::T, δ::T) where T

Helper function for the sum in the Integrand.
"""
E_γ_powδ(Eγ::T, δ::T) where {T} = (2Eγ / E_π0)^δ


"""
    E_integrant(Eγ::T, δ::T, α_δ::T) where T

Helper function for the Integrand.
"""
E_integrant(Eγ::T, δ::T, α_δ::T) where {T} = (E_γ_powδ(Eγ, δ) + E_γ_powδ(Eγ, -δ))^α_δ


"""
    gamma_source_PE04(Eγ::T, f_p::T, q_p::T, nH::T) where {T}

Source function of gamma-ray photons at energy `Eγ` in units of `N_photons erg^-1 s^-1 cm^-3` as given in Pfrommer&Enßlin (2004), Eq.19.

"""
function gamma_source_PE04(Eγ::Real, n_crp::Real, α_p::Real, nH::Real)

    # slope of γ-spectrum
    α = α_γ(α_p)

    # compute expensive stuff here once
    δ = δ_γ(α)

    # α/δ
    α_δ = -α / δ

    # Werhahn+21, B1
    σ_pp(α_γ(α_p)) * c_light * nH * 2^(2 - α) * n_crp * 4 / 3α * (1 / E_π0)^α * E_integrant(Eγ, δ, α_δ)
end



"""
    ñ_crp(rho_cgs::Real, T_K::Real, α_p::Real, Xcr::Real=0.01; 
            xH::Real=0.76)

CR proton normalisation in [1/cm^3]
"""
function ñ_crp(rho_cgs::Real, T_K::Real, α_p::Real, Xcr::Real=0.01;
            xH::Real=0.76)

    # mean molucular weight in hydr. mass
    umu = 4 / (5xH + 3)

    # energy density in gas [g/cm^3/g*K*erg/K] -> PE04, Eq. 10
    ϵNtherm = rho_cgs / (umu * m_p) * T_K * k_B
    α_m1 = α_p - 1

    # PE04, Eq. 8 solved for n_CRp
    return Xcr * ϵNtherm / (m_p * c_light^2) * 2 * α_m1 * (E_p0)^α_m1 /
           beta(0.5 * α_p - 1, (3 - α_p) / 2)

end

"""
    get_number_densities(rho_cgs::Real, T_K::Real, α_p::Real, Xcr::Real, xH::Real)

Thermal proton number density and CR proton normalisation.
"""
function get_number_densities(rho_cgs::Real, T_K::Real, α_p::Real, Xcr::Real, xH::Real)

    # number density of thermal protons [1/cm^3]
    nH = xH * rho_cgs / m_p

    # CR proton normalisation [1/cm^3]
    ñ = ñ_crp(rho_cgs, T_K, α_p, Xcr; xH)

    return nH, ñ
end


"""
    qγ_PE04(rho_cgs, T_K, α_p, Eγ; 
            Xcr::Real=0.5, 
            xH=0.752)

Gamma-ray sorce function at photon energy `Eγ` for thermal gas with properties `rho_cgs` [g/cm^3] and `T_K` [K].
Sets up a CR proton spectrum with energy slope `α_p` as a fraction `Xcr` of the thermal energy density.
Returns source function in units [γ cm^-3 s^-1 GeV^-1].
See Pfrommer&Enßlin (2004), Eq. 19.
"""
function qγ_PE04(rho_cgs::Real, T_K::Real, α_p::Real, Eγ::Real; 
                Xcr::Real=0.5, xH=0.752)

    # number densities of thermal protons and CR protons in [1/cm^3]
    nH, ñ = get_number_densities(rho_cgs, T_K, α_p, Xcr, xH)

    return gamma_source_PE04(Eγ, ñ, α_p, nH)
end

"""
    jγ_PE04(rho_cgs::Real, T_K::Real, α_p::Real, Eγ::Real; 
            Xcr::Real=0.5, xH::Real=0.752)

Gamma-ray emissivity at photon energy `Eγ` [GeV] for thermal gas with properties `rho_cgs` [g/cm^3] and `T_K` [K].
Sets up a CR proton spectrum with energy slope `α_p` as a fraction `Xcr` of the thermal energy density.
Returns emissivity in units [GeV cm^-3 s^-1 ].
See Pfrommer&Enßlin (2004), Eq. 19.

## Function Arguments:
- `rho_cgs`: SPH particle density in [g/cm^3]
- `T_K`: SPH particle temperature [K]
- `α_p`: Slope of proton energy spectrum `S ~ 2.0 - 2.5`
- `Eγ`: Photon energy [GeV]
- `Xcr`: CR proton to thermal pressure ratio.
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings

For mean value along line-of-sight:
- `weights`: `rho` (weight with density)
- `reduce_image`: `true`

For integral along line-of-sight, aka surface brightness:
- `weights`: [`part_weight_physical`](@ref)
- `reduce_image`: `false`
"""
function jγ_PE04(rho_cgs::Real, T_K::Real, α_p::Real, Eγ::Real; 
                Xcr::Real=0.5, xH::Real=0.752)

    Eγ * qγ_PE04(rho_cgs, T_K, α_p, Eγ; Xcr, xH)
end


"""
    λγ_PE04(rho_cgs::Real, T_K::Real, α_p::Real; 
            Xcr::Real=0.5,
            Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
            xH::Real=0.752)

Number of γ-ray photons produced per time and volume from a proton spectrum given as a fraction `Xcr` of the energy density defined by `rho_cgs` [g/cm^3] and `T_K` [K], with a powerlaw slope in energy `α_p`.
Integrated between photon energies `Eγ_π0_min` and `Eγ_π0_max` [GeV].
Returns number of photons in energy band in units of [γ cm^-3 s^-1].
See Pfrommer&Enßlin (2004), Eq. 25.

## Function Arguments:
- `rho_cgs`: SPH particle density in [g/cm^3]
- `T_K`: SPH particle temperature [K]
- `α_p`: Slope of proton energy spectrum `S ~ 2.0 - 2.5`
- `Xcr`: CR proton to thermal pressure ratio.
- `Eγ_π0_min`: Minimum photon energy for γ-ray spectrum [GeV]
- `Eγ_π0_max`: Maximum photon energy for γ-ray spectrum [GeV]
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings

For mean value along line-of-sight:
- `weights`: `rho` (weight with density)
- `reduce_image`: `true`

For integral along line-of-sight, aka surface brightness:
- `weights`: [`part_weight_physical`](@ref)
- `reduce_image`: `false`
"""
function λγ_PE04(rho_cgs::Real, T_K::Real, α_p::Real; 
                Xcr::Real=0.5,
                Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
                xH::Real=0.752)

    # integrand for Eq. 25
    qγ(Eγ) = qγ_PE04(rho_cgs, T_K, α_p, Eγ; Xcr, xH)

    # integrate over energy spectrum
    λγ, int_error = quadgk(qγ, Eγ_π0_min, Eγ_π0_max)

    return λγ
end


"""
    gamma_luminosity_pions_PE04(rho_cgs::Real, m_cgs::Real, T_K::Real, α_p::Real;
                                Xcr::Real=0.5,
                                Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
                                xH::Real=0.752)


γ-ray luminosity produced from a proton spectrum given as a fraction `Xcr` of the energy density defined by `rho_cgs` [g/cm^3] and `T_K` [K], 
with a powerlaw slope in energy `α_p`. Integrated over SPH particle volume for particle of mass `m_cgs` [g].
Integrated between photon energies `Eγ_π0_min` and `Eγ_π0_max` [GeV].
Returns total luminosity integrated over energy band in units of [GeV s^-1].
See Pfrommer&Enßlin (2004), Eq. 25.

## Arguments:
- `rho_cgs`: SPH particle density in [g/cm^3]
- `m_cgs`: SPH particle mass in [g]
- `T_K`: SPH particle temperature [K]
- `α_p`: Slope of proton energy spectrum `S ~ 2.0 - 2.5`
- `Xcr`: CR proton to thermal pressure ratio.
- `Eγ_π0_min`: Minimum photon energy for γ-ray spectrum [GeV]
- `Eγ_π0_max`: Maximum photon energy for γ-ray spectrum [GeV]
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings
- weight function: [`part_weight_one`](@ref)
- reduce image: `true`
"""
function gamma_luminosity_pions_PE04(rho_cgs::Real, m_cgs::Real, T_K::Real, α_p::Real;
                                    Xcr::Real=0.5,
                                    Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
                                    xH::Real=0.752)

    # integrand for Eq. 25
    jγ(Eγ) = jγ_PE04(rho_cgs, T_K, α_p, Eγ; Xcr, xH)

    # integrate over energy spectrum
    Iγ, int_error = quadgk(jγ, Eγ_π0_min, Eγ_π0_max)
    
    # Lγ = Iγ * Vol
    Iγ * m_cgs / rho_cgs
end


"""
    gamma_flux_pions_PE04(rho_cgs::Real, m_cgs::Real, T_K::Real, α_p::Real, d::Real;
                          Xcr::Real=0.5,
                          Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
                          xH::Real=0.752)


Flux of γ-ray photons produced from a proton spectrum given as a fraction `Xcr` of the energy density defined by `rho_cgs` [g/cm^3] and `T_K` [K], 
with a powerlaw slope in energy `α_p`. 
Integrated over SPH particle volume for particle of mass `m_cgs` [g].
Flux from a distance `d` [cm].
Integrated between photon energies `Eγ_π0_min` and `Eγ_π0_max` [GeV].
Returns total number of photons in energy band in untis of [γ cm^-2 s^-1].
See Pfrommer&Enßlin (2004), Eq. 25.

## Arguments:
- `rho_cgs`: SPH particle density in [g/cm^3]
- `m_cgs`: SPH particle mass in [g]
- `T_K`: SPH particle temperature [K]
- `α_p`: Slope of proton energy spectrum `S ~ 2.0 - 2.5`
- `d`: Distance to SPH particle or halo [cm].
- `Xcr`: CR proton to thermal pressure ratio.
- `Eγ_π0_min`: Minimum photon energy for γ-ray spectrum [GeV]
- `Eγ_π0_max`: Maximum photon energy for γ-ray spectrum [GeV]
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings
- weight function: [`part_weight_one`](@ref)
- reduce image: `true`
"""
function gamma_flux_pions_PE04(rho_cgs::Real, m_cgs::Real, T_K::Real, α_p::Real, d::Real;
                                Xcr::Real=0.5,
                                Eγ_π0_min::Real=0.1, Eγ_π0_max::Real=200.0,
                                xH::Real=0.752)
    # volume of SPH particle
    V = m_cgs / rho_cgs

    # flux is total luminosity divided by surface of sphere with radius d
    V / (4π * d^2) * λγ_PE04(rho_cgs, T_K, α_p; Xcr, Eγ_π0_min, Eγ_π0_max, xH)
end