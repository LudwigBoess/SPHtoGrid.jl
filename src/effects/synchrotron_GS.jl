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
Uses the implementaion from [Donnert et. al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.2014D/abstract).

Returns synchrotron emissivity `j_nu` in units [erg/s/Hzcm^3].

# Arguments
- `rho_cgs::Array{<:Real}`: Density in ``g/cm^3``.
- `B_cgs::Array{<:Real}`:   Magnetic field in Gauss.
- `T_K::Array{<:Real}`:     Temperature in Kelvin.
- `Mach::Array{<:Real}`:    Mach number.
- `θ_B::Union{Nothing,Array{<:Real}}=nothing`: Shock obliquity (optional).

## Keyword Arguments
- `xH::Float64 = 0.76`:        Hydrogen fraction of the simulation, if run without chemical model.
- `ν0::Real=1.44e9`:           Observation frequency in ``Hz``.
- `dsa_model::Integer=1`:      Diffusive Shock Acceleration model. Takes values `0...4`, see next section.
- `K_ep::Real=0.01`:           Ratio of CR proton to electron energy acceleration.
- `show_progress::Bool=false`: Enables a progress bar if set to true

## DSA Models
See [DSAModels.jl](https://github.com/LudwigBoess/DSAModels.jl) for details!
- `0`: [Kang et. al. (2007)](https://ui.adsabs.harvard.edu/abs/2007ApJ...669..729K/abstract)
- `1`: [Kang & Ryu (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...95K/abstract)
- `2`: [Ryu et. al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...883...60R/abstract)
- `3`: [Caprioli & Spitkovsky (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...783...91C/abstract)
- `4`: [Pfrommer et. al. (2006)](https://ui.adsabs.harvard.edu/abs/2006MNRAS.367..113P/abstract)

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function analytic_synchrotron_GS(rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                 T_K::Array{<:Real}, Mach::Array{<:Real},
                                 θ_B::Union{Nothing,Array{<:Real}}=nothing;
                                 xH::Real = 0.76,  ν0::Real = 1.4e9,
                                 dsa_model::Integer=1, K_ep::Real=0.01,
                                 show_progress::Bool=false)

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
    j_nu = Vector{Float64}(undef, length(T_K))

    # enable progress meter
    if show_progress
        p = Progress(length(j_nu))
    end

    @threads for i = 1:length(T_K)

        if !isnothing(θ_B)
            ηB = ηB_acc_e(θ_B[i])
        else
            ηB = 1.0
        end

        n0 = K_ep / K_ep_default * ηB * cre_spec_norm_particle(η_model, Mach[i]) * EpsNtherm(rho_cgs[i], T_K[i], xH=xH)

        if (n0 > 0.0) && (B_cgs[i] > 0.0)

            # spectral index of the electron energy spectrum
            s  = dsa_spectral_index(Mach[i])

            agam = 2^(s / 2 - 1 / 2) * √(3) * gamma(s / 4 - 1 / 12) * gamma(s / 4 + 19 / 12) /
                    (8√(π) * (s + 1)) * gamma((s + 5) / 4) / gamma((s + 7) / 4) # GS eq 3.32

            # GS, eq. 3.31 [erg/s/Hz/cm^3]
            j_nu[i] = agam * q_e^3 / (m_e * c_light^2) * (3 * q_e / (m_e^3 * c_light^3 * 4 * π))^((s - 1) / 2) *
                        n0 * ν0^(-(s - 1) / 2) * B_cgs[i]^((s + 1) / 2)
        else
            j_nu[i] = 0.0
        end

        # update progress meter
        if show_progress
            next!(p)
            flush(stdout)
            flush(stderr)
        end
    end

    return j_nu
end