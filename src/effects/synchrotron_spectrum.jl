using SynchrotronKernel
#using QuadGK
global const X_MAX = 70.0  # -> F(x) < 2e-6
global const X_MIN = 1e-20 # -> F(x) < 5e-6

"""
    min_synch_loud_momentum(ν0::Real, B::Real)

Computes the minimum synchrotron loud momentum as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 22.
"""
min_synch_loud_momentum(ν0::Real, B::Real) = √( ν0 / (C_crit * B)  )


"""
    ν_over_ν_crit(p::Real, B::Real, ν::Real)

Computes the fraction ``x = \frac{ν}{ν_c}`` needed for the synchrotron Kernel.
See Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 19.
"""
ν_over_ν_crit(p::Real, B::Real, ν::Real, sinθ::Real=1.0) = ν / ( C_crit * B * sinθ * p^2 )

"""
    dsa_spectral_index_momentum(M::Real)

Compute the spectral index from DSA for a distribution in momentum space.
"""
function dsa_spectral_index_momentum(M::Real) 
    xs = ( γ + 1 ) * M^2 / ( ( γ - 1 ) * M^2 + 2 )
    return 3.0 * xs / (xs - 1.0 )
end


"""
    spectrum(p::Real, p_min::Real, n0::Real, s::Real)

Energy density `n` for a given momentum `p` interpolated from `n0` at minimum momentum `p_min`.
"""
spectrum(p::Real, p_min::Real, n0::Real, s::Real) = p < p_min ? 0.0 : n0 * (p/p_min)^(-s)


"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17.
"""
function integrate_θ(x_in::Real, θ_steps::Integer=128)

    dθ = 0.5π / θ_steps
    K = 0.0
    @inbounds for θ ∈ LinRange(0.0, 0.5π, θ_steps)
        sinθ = sin(θ)
        x  = x_in / sinθ
        K += sinθ^2 * synchrotron_kernel(x) * dθ
    end
    return K
end

"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17 using Simpson's rule.
"""
function integrate_θ_simpson(x_in::Real, θ_steps::Integer=100)

    dθ = 0.5π / θ_steps
    #θ = LinRange(0.0, 0.5π, θ_steps)

    F     = Vector{Float64}(undef, θ_steps)
    F_mid = Vector{Float64}(undef, θ_steps)

    # sin(0.0) = 0.0
    F[1]     = 0.0
    F_mid[1] = sin(0.5dθ)^2 * synchrotron_kernel(x_in/sin(0.5dθ))
    K = 0.0
    @inbounds for i = 2:θ_steps
        # boundary
        sinθ = sin((i-1)*dθ)
        x  = x_in / sinθ
        F[i] = sinθ^2 * synchrotron_kernel(x)

        # mid point
        sinθ = sin((i-0.5)*dθ)
        x    = x_in / sinθ
        F_mid[i] = sinθ^2 * synchrotron_kernel(x)

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        K += dθ / 6.0 * ( F[i] + F[i-1] + 4F_mid[i] )
    end
    return K
end


"""
    analytic_synchrotron_emission( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                                   T_K::Array{<:Real}, Mach::Array{<:Real};
                                   xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                                   integrate_pitch_angle::Bool=true )

Computes the synchrotron emission for a powerlaw spectrum as described in Donnert+16, MNRAS 462, 2014–2032 (2016).

# Arguments
- `rho_cgs::Array{<:Real}`: Density in ``g/cm^3``.
- `B_cgs::Array{<:Real}`:   Magnetic field in Gauss.
- `T_K::Array{<:Real}`:     Temperature in Kelvin.
- `Mach::Array{<:Real}`:    Mach number.

# Keyword Arguments
- `xH::Float64 = 0.76`:               Hydrogen fraction of the simulation, if run without chemical model.
- `dsa_model::Integer=1`:             Diffuse-Shock-Acceleration model. Takes values `0...4`, see next section.
- `ν0::Real=1.44e9`:                  Observation frequency in ``Hz``.
- `Xcre::Real=0.01`:                  Ratio of CR proton to electron energy density.
- `integrate_pitch_angle::Bool=true`: Integrates over the pitch angle as in Longair Eq. 8.87.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz*s]`` to ``mJy``.

# DSA models
- `0`: [`KR07_acc`](@ref). Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1.
- `1`: [`KR13_acc`](@ref). Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95 .
- `2`: [`Ryu19_acc`](@ref). Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476 .
- `3`: [`CS14_acc`](@ref). Efficiency model from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91 .
- `4`: [`P16_acc`](@ref). Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 .

"""
function spectral_synchrotron_emission(rho_cgs::Real, B_cgs::Real,
                                       T_K::Real, Mach::Real, t::Real=0.0;
                                       xH::Real=0.76, dsa_model::Integer=1, 
                                       ν0::Real=1.44e9,
                                       Xcre::Real=0.01,
                                       Emin::Real=5.0e4,
                                       Emax::Real=1.e10,
                                       p_inj::Real=0.1, # in [me*c]
                                       convert_to_mJy::Bool=false,
                                       N_sample_bins::Integer=128)

    if ( B_cgs == 0 || Mach < 1 )
        return 0
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

    # # in dimensionless momentum!
    p_inj = 0.1 
    p_min = Emin * eV2cgs / ( m_e * c_light^2 )
    p_max = Emax * eV2cgs / ( m_e * c_light^2 )

    # prefactor to Eq. 17
    j_ν_prefac = q_e * √(3) / (m_e * c_light^2)

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # get spectral index 
    #s  = dsa_spectral_index(Mach)
    s  = dsa_spectral_index_momentum(Mach)

    # compute norm CR energy spectrum (in [erg/cm^3])
    # uses p_inj = 0.1 * m_e * c_light
    n0 = Xcre * cre_spec_norm_particle(Mach, acc_function) * EpsNtherm(rho_cgs, T_K, xH=xH)

    println("ϵ_th_dw = $(EpsNtherm(rho_cgs, T_K, xH=xH))")
    println("")
    println("n0 = $n0")
    
    
    # Smac2:
    #n0 = Xcre * ( s - 2 ) * (Emin * eV2cgs)^(s - 2) * EpsNtherm(rho_cgs, T_K, xH=xH)

    # width of momentum bins
    #di    = log(p_max/p_min) / (N_sample_bins-1)
    di    = log10(p_max/p_min) / (N_sample_bins-1)

    # momentum bins
    #p     = @. p_min * exp( di * 0:N_sample_bins-1)
    p     = @. p_min * 10.0^( di * 0:N_sample_bins-1)

    # width of momentum bins
    dp = zeros(N_sample_bins)
    #dp[1] = p[1] - p_min * exp(-1di)
    dp[1] = p[1] - p_min * 10.0^(-1di)
    @inbounds for i = 2:N_sample_bins
        dp[i] = p[i] - p[i-1]
    end

    # storage array for integrand of Eq. 17
    F        = Vector{Float64}(undef, N_sample_bins)
    F_mid    = Vector{Float64}(undef, N_sample_bins)
    F[1]     = synchrotron_kernel(ν_over_ν_crit(p[1], B_cgs, ν0))
    F_mid[1] = synchrotron_kernel(0.5 * ( ν_over_ν_crit(p[1], B_cgs, ν0) + 
                                          ν_over_ν_crit(p[2], B_cgs, ν0) ))

    # store total synchrotron emissivity
    jν = 0.0

    # Simpson rule
    @inbounds for i = 2:N_sample_bins

        # beginning of bin
        # energy density at momentum p
        N_E = spectrum(p[i], p_inj, n0, s)
        # x from Eq. 19 (without sinθ -> for integration)
        x   = ν_over_ν_crit(p[i], B_cgs, ν0)
        K   = integrate_θ_simpson(x)

        # p-Integrand of Eq. 17
        F[i] = N_E * K 

        # middle of bin 
        x_mid = 0.5 * ( x[i-1] + x[i] )
        p_mid = 0.5 * ( p[i-1] + p[i] )

        N_E_mid  = spectrum(p_mid, p_inj, n0, s)
        K_mid    = integrate_θ_simpson(x_mid)

        F_mid[i] = N_E_mid * K_mid

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        jν += dp[i] / 6.0 * ( F[i] + F[i-1] + 4F_mid[i] )

    end

    return j_ν_prefac * jν

end




"""
    analytic_synchrotron_emission( n_p::Vector{<:Real}, 
                                   p::Vector{<:Real},
                                   B_cgs::Real;
                                   ν0::Real=1.44e9,
                                   convert_to_mJy::Bool=false)

Computes the synchrotron emission (in ``[erg/cm^3/Hz*s]``) for a CR spectrum as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

# Arguments
- `n_p::Vector{<:Real}`: Energy density in ``erg/cm^3`` for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for energy densities.
- `B_cgs::Real`:         Magnetic field strength (absolute value).

# Keyword Arguments
- `ν0::Real=1.44e9`:                  Observation frequency in ``Hz``.
- `Xcre::Real=0.01`:                  Ratio of CR proton to electron energy density.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz*s]`` to ``mJy/cm``.

"""
function spectral_synchrotron_emission(n_p::Vector{<:Real}, 
                                       p::Vector{<:Real},
                                       B_cgs::Real;
                                       ν0::Real=1.44e9,
                                       convert_to_mJy::Bool=false)

    if B_cgs == 0
        return 0
    end

    # get the number of momentums for which the energy density is defined
    Nbins = length(p)

    # prefactor to Eq. 17
    j_ν_prefac = q_e * √(3) / (m_e * c_light^2)

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # width of momentum bins
    dp = zeros(N_sample_bins)
    dp[1] = p[1] - p[1] * exp(-1di)
    #dp[1] = p[1] - p_min * 10.0^(-1di)
    @inbounds for i = 2:N_sample_bins
        dp[i] = p[i] - p[i-1]
    end

    # storage array for integrand of Eq. 17
    F        = Vector{Float64}(undef, N_sample_bins)
    F_mid    = Vector{Float64}(undef, N_sample_bins)
    F[1]     = synchrotron_kernel(ν_over_ν_crit(p[1], B_cgs, ν0))
    F_mid[1] = synchrotron_kernel(0.5 * ( ν_over_ν_crit(p[1], B_cgs, ν0) + 
                                          ν_over_ν_crit(p[2], B_cgs, ν0) ))

    # store total synchrotron emissivity
    jν = 0.0

    # Simpson rule
    @inbounds for i = 2:N_bins

        # beginning of bin
        
        # x from Eq. 19
        x   = ν_over_ν_crit(p[i], B_cgs, ν0)
        # pitch-angle integral
        K   = integrate_θ_simpson(x)

        # energy density at momentum p * integrated synchrotron kernel
        F[i] = n_p[i] * K 

        # middle of bin 
        N_E_mid  = 0.5 * ( n_p[i-1] + n_p[i]  )

        p_mid = 0.5 * ( p[i-1] + p[i] )
        x     = ν_over_ν_crit(p_mid, B_cgs, ν0)
        K_mid = integrate_θ_simpson(x)

        F_mid[i] = N_E_mid * K_mid

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        jν += dp[i] / 6.0 * ( F[i] + F[i-1] + 4F_mid[i] )

    end

    return jν * j_ν_prefac

end 