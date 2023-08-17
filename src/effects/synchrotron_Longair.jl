
"""
    analytic_synchrotron_emission( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                   T_K::Array{<:Real}, Mach::Array{<:Real};
                                   xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                   integrate_pitch_angle::Bool=true )
Computes the analytic synchrotron emission with the simplified approach described in Longair Eq. 8.128.
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
# DSA models imported from DiffusiveShockAccelerationModels.jl
- `0`: Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1.
- `1`: Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95 .
- `2`: Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476 .
- `3`: Efficiency model from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91 .
- `4`: Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 .
"""
function analytic_synchrotron_Longair(rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                        T_K::Array{<:Real}, Mach::Array{<:Real};
                                        xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                        K_ep::Real=0.01,
                                        integrate_pitch_angle::Bool=true,
                                        convert_to_mJy::Bool=false)

    Npart = length(T_K)

    # default ratio for electron to protons from Donnert+16.
    K_ep_default = 0.01

    # allocate storage array
    J_ν = Vector{Float64}(undef, Npart)

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

    # bracket of Longair eq. 8.128
    nufac = (3q_e) / (m_e^3 * c_light * c_light * c_light * c_light * c_light * 2π * ν0)


    @inbounds for i = 1:Npart

        B = B_cgs[i]
        n0 = K_ep / K_ep_default * cre_spec_norm_particle(η_model, Mach[i]) * EpsNtherm(rho_cgs[i], T_K[i], xH=xH)
        s = dsa_spectral_index(Mach[i])

        if n0 > 0.0

            # Longair eq. 8.129
            a_p = gamma(s / 4 + 19 / 12) * gamma(s / 4 - 1 / 12)

            if integrate_pitch_angle
                # Longair eq 8.87
                a_p *= 0.5 * √(π) * gamma(s / 4 + 5 / 4) /
                       gamma(s / 4 + 7 / 4)
            end

            # Longair eq 8.128 prefactor
            prefac = √(3) * q_e^3 / (m_e * c_light^2 * (s + 1))

            # Longair eq 8.128
            J_ν[i] = prefac * n0 *
                     nufac^(0.5 * (s - 1)) * B^(0.5 * (s + 1)) *
                     a_p
        else
            J_ν[i] = 0.0
        end
    end

    if convert_to_mJy
        J_ν .*= mJy_factor
    end

    return J_ν
end