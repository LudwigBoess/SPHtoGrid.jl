using SpecialFunctions
using LinearAlgebra
using ProgressMeter

"""
    Helper functions
"""
# power functions
@inline p2(x) = x*x
@inline p3(x) = x*x*x

"""
    kr_fitting_function(x::Real, 
                             a0::Real, a1::Real, a2::Real, a3::Real, a4::Real)

Helper function to use the fitting function from KR07.
"""
@inline function kr_fitting_function(x::Real, 
                             a0::Real, a1::Real, a2::Real, a3::Real, a4::Real)
    mm = x - 1.0

	return ( a0 + a1*mm + a2*p2(mm) + a3*p3(mm) + a4*p2(p2(mm)) ) / p2(p2(x))
end

"""
    Efficiency functions
"""

"""
    KR07_acc(M::Float64)

Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1
"""
function KR07_acc(M::Real)
    if M <= 2.0
        return 1.96e-3*(M^2 - 1.)             # eq. A3
    else
        return kr_fitting_function(M, 5.46, -9.78, 4.17, -0.337, 0.57)
    end
end

"""
    KR13_acc(M::Real)

Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95
"""
function KR13_acc(M::Real)

    if M < 2.0
        return 0.0
    elseif 2.0 <= M <= 5.0
        return -0.0005950569221922047 + 1.880258286365841e-5 * M^5.334076006529829 
    elseif 5.0 < M <= 15.0
        return kr_fitting_function(M, -2.8696966498579606, 9.667563166507879,
                                  -8.877138312318019, 1.938386688261113, 0.1806112438315771)
    else
        return 0.21152
    end
end

"""
    Ryu19_acc(M::Real)

Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476
Values for 2.25 < M <= 5.0 extrapolated to entire range
"""
function Ryu19_acc(M::Real)

    if M < 2.25
        return 0.0
    elseif M <= 34.0
        return kr_fitting_function(M, -1.5255114554627316, 2.4026049650156693,
                 -1.2534251472776456, 0.22152323784680614, 0.0335800899612107)
    else
        return 0.0348
    end
end

"""
    CS14_acc(M::Real)

Efficiency from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91
Same simplified approach as Vazza+12 -> is roughly half the efficiency of Kang&Ryu 2013.
"""
function CS14_acc(M::Real)
    vazza_factor = 0.5
    return vazza_factor * KR13_acc(M)
end

"""
    P16_acc(M::Real)

Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 
"""
function P16_acc(M::Real)
    return 0.5
end


"""
    cre_spec_norm_particle(M::Real, eff_function::Function)

Computes the CR electron norm of the particles. 
This depends on the Mach number `M` and the acceleration efficiency given by `eff_function`.
"""
function cre_spec_norm_particle(M::Real, eff_function::Function)
    norm = eff_function(M)
    return norm / 7.6e14 # Donnert et al 2016, eq. 40, p_0 = 0.1 me c 
end

"""
    dsa_spectral_index(M::Real)

Spectral index given by standard Diffuse-Shock-Acceleration (DSA).
"""
function dsa_spectral_index(M::Real)
    return 2.0 * ( M^2 + 1.0 ) / ( M^2 - 1 )
end

"""
    EpsNtherm(rho_cgs::Real, T_K::Real)

Thermal energy density in cgs.
"""
function EpsNtherm(rho_cgs::Real, T_K::Real; xH::Real=0.76)
    u_mol = (4.0/(5.0*xH+3.0))
    return (rho_cgs / (m_p * u_mol) * k_B * T_K)
end


"""
    analytic_synchrotron_emission( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                   T_K::Array{<:Real}, Mach::Array{<:Real};
                                   xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                   integrate_pitch_angle::Bool=true )

Computes the analytic synchrotron emission with the simplified approach described in Longair Eq. 8.128.

# Arguments
- `rho_cgs::Array{<:Real}`: Density in ``g/cm^3``.
- `B_cgs::Array{<:Real}`:   Magnetic field in Gauss.
- `T_K::Array{<:Real}`:     Temperature in Kelvin.
- `Mach::Array{<:Real}`:    Mach number.

# Keyword Arguments
- `xH::Float64 = 0.76`:               Hydrogen fraction of the simulation, if run without chemical model.
- `dsa_model::Integer=1`:             Diffuse-Shock-Acceleration model. Takes values `0...4`, see next section.
- `ν0::Real=1.44e9`:                   Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=true`: Integrates over the pitch angle as in Longair Eq. 8.87.

# DSA models
- `0`: [`KR07_acc`](@ref). Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1.
- `1`: [`KR13_acc`](@ref). Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95 .
- `2`: [`Ryu19_acc`](@ref). Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476 .
- `3`: [`CS14_acc`](@ref). Efficiency model from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91 .
- `4`: [`P16_acc`](@ref). Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 .

"""
function analytic_synchrotron_emission(rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                       T_K::Array{<:Real}, Mach::Array{<:Real};
                                       xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                       integrate_pitch_angle::Bool=true)

    Npart = length(T_K)

    S_ν = zeros(Npart)

    if length(B_cgs[1,:]) == 1
        B_1dim = true
    else
        B_1dim = false
    end

    if dsa_model == 0
        acc_function = KR07_acc
    elseif dsa_model == 1
        acc_function = KR13_acc
    elseif dsa_model == 2
        acc_function = Ryu19_acc
    elseif dsa_model == 3
        acc_function = CS14_acc
    elseif dsa_model == 4
        acc_function = P16_acc
    else
        error("Invalid DSA model selection!")
    end

    nufac = (3q_e)/(p3(m_e) * c_light * c_light * c_light * c_light * c_light * 2π * ν0)

    @showprogress for i = 1:Npart
        if B_1dim
            B = B_cgs[i]
        else
            B  = sqrt( B_cgs[i,1]^2 + B_cgs[i,2]^2 + B_cgs[i,3]^2)
        end
        n0 = cre_spec_norm_particle(Mach[i], acc_function) * EpsNtherm(rho_cgs[i], T_K[i], xH=xH)
        s  = dsa_spectral_index(Mach[i])

        if n0 > 0.0
            # Longair eq 8.128
            prefac = sqrt( 3.0 )*p3(q_e)/(m_e * p2(c_light)*(s + 1.0)) * 
                    gamma(s / 4.0 + 19.0 / 12.0) * gamma(s / 4.0 - 1.0 / 12.0) *
                    nufac^( 0.5 * ( s - 1.0 ) )
    

            if integrate_pitch_angle
                # Longair eq 8.87
                prefac *= 0.5 * sqrt(π) * gamma(s / 4.0 + 5.0 / 4.0 ) / 
                                        gamma(s / 4.0 + 7.0 / 4.0 )
            end

            S_ν[i] = prefac * n0 * B^( 0.5 * ( s + 1.0 ) )
        else
            S_ν[i] = 0.0
        end
    end

    return S_ν
end
