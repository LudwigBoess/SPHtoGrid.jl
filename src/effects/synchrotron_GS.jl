"""
    cre_spec_norm_particle(M::Real, η_model::AbstractShockAccelerationEfficiency)

Computes the CR electron norm of the particles. 
This depends on the Mach number `M` and the acceleration efficiency given by `η_model`.
"""
function cre_spec_norm_particle(η_model::AbstractShockAccelerationEfficiency, M::T) where {T}
    η_Ms_acc(η_model, M) / 7.6e14 # Donnert et al 2016, eq. 40, p_0 = 0.1 me c 
end

"""
    dsa_spectral_index(M::Real)

Spectral index given by standard Diffuse-Shock-Acceleration (DSA).
"""
function dsa_spectral_index(M::Real)

    # diverges for Mach smaller 1.01
    if M <= 1.01
        return maxintfloat(Float64)
    else
        return 2 * (M^2 + 1) / (M^2 - 1) 
    end
end


"""
    EpsNtherm(rho_cgs::Real, T_K::Real)

Thermal energy density in cgs.
"""
function EpsNtherm(rho_cgs::Real, T_K::Real; xH::Real = 0.76)
    u_mol = (4 / (5xH + 3))
    return (rho_cgs / (m_p * u_mol) * k_B * T_K)
end



"""
    analytic_synchrotron_GS( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                             T_K::Array{<:Real}, Mach::Array{<:Real};
                             xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                             integrate_pitch_angle::Bool=true )

Computes the analytic synchrotron emission with the simplified approach described in Ginzburg & Syrovatskii 1965, "Cosmic Magnetobremsstrahlung".
Uses the implementaion from Donnert J. M. F., Stroe A., Brunetti G., Hoang D., Roettgering H., 2016, MNRAS, 462, 2014.
Returns J_ν in units [erg/s/Hzcm^3].

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

# DSA models imported from DSAModels.jl
- `0`: Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1.
- `1`: Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95 .
- `2`: Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476 .
- `3`: Efficiency model from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91 .
- `4`: Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 .
"""
function analytic_synchrotron_GS(rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                 T_K::Array{<:Real}, Mach::Array{<:Real};
                                 xH::Real = 0.76, dsa_model::Integer = 1, ν0::Real = 1.4e9,
                                 K_ep::Real = 0.01)

    # default ratio for electron to protons from Donnert+16, already accounted for in spectral norm
    K_ep_default = 0.01

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

    # allocate storage array
    J_ν = Vector{Float64}(undef, length(T_K))


    @threads for i = 1:length(T_K)

        n0 = K_ep / K_ep_default * cre_spec_norm_particle(η_model, Mach[i]) * EpsNtherm(rho_cgs[i], T_K[i], xH = xH)

        if (n0 > 0.0) && (B_cgs[i] > 0.0)

            # spectral index of the electron energy spectrum
            s  = dsa_spectral_index(Mach[i])

            agam = 2^(s / 2 - 1 / 2) * √(3) * gamma(s / 4 - 1 / 12) * gamma(s / 4 + 19 / 12) /
                    (8√(π) * (s + 1)) * gamma((s + 5) / 4) / gamma((s + 7) / 4) # GS eq 3.32

            # GS, eq. 3.31 [erg/s/Hz/cm^3]
            J_ν[i] = agam * q_e^3 / (m_e * c_light^2) * (3 * q_e / (m_e^3 * c_light^3 * 4 * π))^((s - 1) / 2) *
                n0 * ν0^(-(s - 1) / 2) * B_cgs[i]^((s + 1) / 2)
        else
            J_ν[i] = 0.0
        end
    end

    return J_ν
end