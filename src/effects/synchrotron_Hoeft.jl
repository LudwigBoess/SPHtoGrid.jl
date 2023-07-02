"""
    B_cmb(z::Real)

Magnetic field equivalent of the CMB for a given redshift `z` in `[G]`.
"""
B_cmb(z::Real) = 3.2*(1 + z)^2 * 1.e-6

"""
    analytic_synchrotron_HB07( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                   T_K::Array{<:Real}, Mach::Array{<:Real};
                                   xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                   integrate_pitch_angle::Bool=true )

Computes the analytic synchrotron emission with the simplified approach described in Hoeft&Brüggen (2007), following approach by Wittor et. al. (2017).
Returns J_ν in units [erg/cm^3/Hz/s].

# Arguments
- `rho_cgs::Array{<:Real}`: Density in ``g/cm^3``.
- `B_cgs::Array{<:Real}`:   Magnetic field in Gauss.
- `T_K::Array{<:Real}`:     Temperature in Kelvin.
- `Mach::Array{<:Real}`:    Mach number.

# Keyword Arguments
- `xH::Float64 = 0.76`:               Hydrogen fraction of the simulation, if run without chemical model.
- `dsa_model::Integer=1`:             Diffuse-Shock-Acceleration model. Takes values `0...4`, see next section.
- `ν0::Real=1.44e9`:                  Observation frequency in ``Hz``.
- `K_ep::Real=0.01`:                  Ratio of CR proton to electron energy density.
- `integrate_pitch_angle::Bool=true`: Integrates over the pitch angle as in Longair Eq. 8.87.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
"""
function analytic_synchrotron_HB07(rho_cgs::Array{<:Real}, m_cgs::Array{<:Real}, hsml_cgs::Array{<:Real},
                                    B_cgs::Array{<:Real}, T_keV::Array{<:Real}, Mach::Array{<:Real};
                                    xH::Real = 0.752, ν0::Real = 1.4e9,
                                    ξe::Real = 0.01, z::Real=0.0,
                                    dsa_model::Integer=1)

    # select DSA model
    if dsa_model == 0
        η_model = Kang07()
    elseif dsa_model == 1
        η_model = KR13()
    elseif dsa_model == 2
        η_model = Ryu19()
    elseif dsa_model == 3
        η_model = CS14()
    elseif dsa_model == 4
        η_model = P16()
    else
        error("Invalid DSA model selection!")
    end

    # prefactor to Eq. 32
    j_ν_prefactor = 6.4e34 # erg/s/Hz

    # electron energy density ratio 
    ξe_factor = ξe/0.05

    # factor to convert from [g/cm^3] to [1/cm^3]
    n2ne = (xH + 0.5 * (1 - xH)) /
           (2xH + 0.75 * (1 - xH))
    umu = 4 / (5xH + 3)
    rho2ne = n2ne / (umu * m_p)

    # storage array for synchrotron emissivity
    j_nu = Vector{Float64}(undef, length(rho_cgs))

    @threads for i = 1:length(rho_cgs)

        # particle depth in [cm]
        dz = get_dz(m_cgs[i], rho_cgs[i], hsml_cgs[i])

        # particle volume in [cm^3]
        V = m_cgs[i] / rho_cgs[i]

        # particle area in [Mpc^2]
        A = V / dz / (1.e3 * kpc)^2

        # spectral index of electrons
        s = dsa_spectral_index(Mach[i])

        # electron density
        ne = rho2ne * rho_cgs[i] * 1.e4

        # scaling with magnetic field
        Bfactor = (B_cgs[i] * 1.e6)^(1 + 0.5 * s) / ((B_cmb(z) * 1.e6)^2 + (B_cgs[i] * 1.e6)^2)
        println(B_cmb(z))

        # Wittor+17, Eq. 16, but in untis erg/s/Hz/cm^3
        j_nu[i] = j_ν_prefactor * A * ne * ξe_factor * (ν0/1.4e9)^(-0.5s) * 
                  √(T_keV[i] / 7.0)^3 * Bfactor * η_Ms_acc(η_model, Mach[i]) / V
    end

    return j_nu
end

