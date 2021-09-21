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

Number density `n` for a given momentum `p` interpolated from `n0` at minimum momentum `p_min`.
"""
spectrum(p::Real, p_min::Real, n0::Real, s::Real) = p < p_min ? 0.0 : n0 * (p/p_min)^(-s)


"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17.
"""
function integrate_θ(x_in::Real, θ_steps::Integer=64)

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
function integrate_θ_simpson(x_in::Real, θ_steps::Integer=50)

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
    get_F( p::Real, p_inj::Real, ϵ_cr0::Real, s::Real, B_cgs::Real, ν0::Real)

Helper function to get synchrotron emission at momentum `p`.
"""
function get_F( p::Real, p_inj::Real, ϵ_cr0::Real, s::Real, B_cgs::Real, ν0::Real, integrate_pitch_angle::Bool)

    # energy density at momentum p
    N_E = spectrum(p, p_inj, ϵ_cr0, s)
    # x from Eq. 19 (without sinθ -> for integration)
    x   = ν_over_ν_crit(p, B_cgs, ν0)

    if integrate_pitch_angle
        K = integrate_θ_simpson(x)
    else
        K = synchrotron_kernel(x)
    end

    # integral Donnert+16, Eq. 17 evaluated at momenum p
    return N_E * K 

end

"""
    get_log_mid(p_low::Real, p_up::Real)

Helper function to compute the momentum in the middle of the bin.
"""
get_log_mid(q_low::Real, q_up::Real) = 10.0^( 0.5 * ( log10(q_low) + log10(q_up) ) )

"""
    spectral_synchrotron_emission(rho_cgs::Real, B_cgs::Real,
                                  T_K::Real, Mach::Real;
                                  xH::Real=0.76, dsa_model::Integer=1, 
                                  ν0::Real=1.44e9,
                                  K_ep::Real=0.01,
                                  Emin::Real=5.0e4,
                                  Emax::Real=1.e10,
                                  p_inj::Real=0.1, # in [me*c]
                                  integrate_pitch_angle::Bool=true,
                                  convert_to_mJy::Bool=false,
                                  N_sample_bins::Integer=128)

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
- `K_ep::Real=0.01`:                  Ratio of CR proton to electron energy density.
- `integrate_pitch_angle::Bool=true`: Explicitly integrates over the pitch angle. IF `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.

# DSA models
- `0`: [`KR07_acc`](@ref). Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1.
- `1`: [`KR13_acc`](@ref). Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95 .
- `2`: [`Ryu19_acc`](@ref). Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476 .
- `3`: [`CS14_acc`](@ref). Efficiency model from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91 .
- `4`: [`P16_acc`](@ref). Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 .

"""
function spectral_synchrotron_emission(rho_cgs::Real, B_cgs::Real,
                                       T_K::Real, Mach::Real;
                                       xH::Real=0.76, dsa_model::Integer=1, 
                                       ν0::Real=1.44e9,
                                       K_ep::Real=0.01,
                                       Emin::Real=5.0e4,
                                       Emax::Real=1.e10,
                                       p_inj::Real=0.1, # in [me*c]
                                       integrate_pitch_angle::Bool=true,
                                       convert_to_mJy::Bool=false,
                                       N_sample_bins::Integer=128)

    if ( B_cgs == 0 || Mach < 2 )
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

    # in dimensionless momentum!
    p_inj = 0.1 
    p_min = Emin * eV2cgs / ( m_e * c_light^2 )
    p_max = Emax * eV2cgs / ( m_e * c_light^2 )

    # prefactor to Eq. 17
    j_ν_prefac = q_e * √(3) / (m_e * c_light^2)

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # get spectral index s ∈ ]-∞, 4.0[
    s  = dsa_spectral_index_momentum(Mach)

    # energy density of thermal gas [erg/cm^3]
    ϵ_th = EpsNtherm(rho_cgs, T_K, xH=xH)

    # ! this gives the same result as the analytic solution by Longair
    # Normalisation of powerlaw from thermal energy density and CR acceleration efficiency
    # as in Donnert+16, Eq. 60. Uses p_inj = 0.1 me c
    ϵ_cr0 = cre_spec_norm_particle(Mach, acc_function) * ϵ_th

    # ! This gives the same result as in Smac2
    # # CR Energy density of whole powerlaw
    # ϵ_cr0 = K_ep * get_rel_energy_density(Mach, acc_function) * ϵ_th
    # # compensate for injection momentum
    # ϵ_cr0 *= (s - 2) / (p_inj^( 2 - s ) )
    # ! This gives the same result as in Smac2


    # width of momentum bins
    di    = log10(p_max/p_min) / (N_sample_bins-1)

    # momentum bins
    p     = @. p_min * 10.0^( di * 0:N_sample_bins-1)

    # width of momentum bins
    dp = zeros(N_sample_bins)
    dp[1] = p[1] - p_min * 10.0^(-1di)
    @inbounds for i = 2:N_sample_bins
        dp[i] = p[i] - p[i-1]
    end

    # storage array for integrand of Eq. 17
    F        = Vector{Float64}(undef, N_sample_bins)
    F_mid    = Vector{Float64}(undef, N_sample_bins)

    
    F[1]     = get_F(p[1], p_inj, ϵ_cr0, s, B_cgs, ν0, integrate_pitch_angle)
    # not needed
    F_mid[1] = 0.0

    # store total synchrotron emissivity
    jν = 0.0

    # Simpson rule
    @inbounds for i = 2:N_sample_bins

        # end of bin
        F[i]     = get_F(p[i], p_inj, ϵ_cr0, s, B_cgs, ν0, integrate_pitch_angle)

        # middle of bin
        p_mid    = get_log_mid(p[i-1], p[i])
        F_mid[i] = get_F(p_mid, p_inj, ϵ_cr0, s, B_cgs, ν0, integrate_pitch_angle)

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        jν += dp[i] / 6.0 * ( F[i] + F[i-1] + 4F_mid[i] )

    end

    if convert_to_mJy
        j_ν_prefac *= mJy_factor
    end

    return j_ν_prefac * jν 

end




"""
    spectral_synchrotron_emission( n_p::Vector{<:Real}, 
                                   p::Vector{<:Real},
                                   B_cgs::Real;
                                   ν0::Real=1.44e9
                                   integrate_pitch_angle::Bool=false,
                                   convert_to_mJy::Bool=false)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a CR spectrum as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

# Arguments
- `n_p::Vector{<:Real}`: Number density in ``1/cm^3`` for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for number densities.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.

# Keyword Arguments
- `ν0::Real=1.44e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.

"""
function spectral_synchrotron_emission(n_p::Vector{<:Real}, 
                                       p::Vector{<:Real},
                                       B_cgs::Real;
                                       ν0::Real=1.44e9,
                                       integrate_pitch_angle::Bool=false,
                                       convert_to_mJy::Bool=false)

    if ( B_cgs == 0 ) || sum(n_p) == 0
        return 0
    end

    # avoid error for n_p = -0.0
    n_p[ n_p .== -0.0 ] .= 0.0
    
    # if there are Infs or NaNs present we get a wrong result, so skip this particle
    # ToDo: Check if there is a better way to fix this! 
    if (length(findall( isnan.(n_p) )) > 0) ||
       (length(findall( isinf.(n_p) )) > 0) ||
       (length(findall(   n_p .< 0  )) > 0)
        return 0
    end

    # get the number of momentums for which the energy density is defined
    Nbins = length(n_p)

    # prefactor to Eq. 17
    j_ν_prefac = q_e * √(3) / (m_e * c_light^2)

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # width of momentum bins
    dp = Vector{Float64}(undef, Nbins)
    @inbounds for i = 1:Nbins
        dp[i] = p[i+1] - p[i]
    end

    # storage array for integrand of Eq. 17
    F        = Vector{Float64}(undef, Nbins)
    F_mid    = Vector{Float64}(undef, Nbins-1)

    @inbounds for i = 1:Nbins-1

        # beginning of bin
        
        # x from Eq. 19
        x   = ν_over_ν_crit(p[i], B_cgs, ν0)
        # pitch-angle integral
        if integrate_pitch_angle
            K = integrate_θ_simpson(x)
        else
            K = synchrotron_kernel(x)
        end

        # energy density at momentum p * integrated synchrotron kernel
        F[i] = n_p[i] * K 

        # middle of bin 
        N_E_mid  = get_log_mid( n_p[i],  n_p[i+1] )

        p_mid = get_log_mid( p[i], p[i+1] )
        x     = ν_over_ν_crit(p_mid, B_cgs, ν0)

        if integrate_pitch_angle
            K_mid = integrate_θ_simpson(x)
        else
            K_mid = synchrotron_kernel(x)
        end

        F_mid[i] = N_E_mid * K_mid

    end

    # last bin seperate
    x   = ν_over_ν_crit(p[end], B_cgs, ν0)
    # pitch-angle integral
    if integrate_pitch_angle
        K = integrate_θ_simpson(x)
    else
        K = synchrotron_kernel(x)
    end

    # energy density at momentum p * integrated synchrotron kernel
    F[end] = n_p[end] * K 

    # store total synchrotron emissivity
    jν = 0.0
    @inbounds for i = 1:Nbins-1
        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        jν += dp[i] / 6.0 * ( F[i] + F[i+1] + 4F_mid[i] )
    end

    if convert_to_mJy
        j_ν_prefac *= mJy_factor
    end

    return j_ν_prefac * jν

end 
