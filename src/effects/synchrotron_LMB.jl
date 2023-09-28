using SynchrotronKernel


"""
    powerlaw_spectrum(E::Real, s::Real, CR_Emin::Real)

Interpolate in powerlaw spectrum.
"""
function powerlaw_spectrum(E::Real, s::Real, CR_Emin::Real)
    Nprefac = (s - 2) * CR_Emin^(s - 2)
    return Nprefac / E^s
end

"""
    ys_f(xs::Real, ξ::Real)

Helper function to get `ys` (see Pfrommer+2017)
"""
function ys_f(xs::Real, ξ::Real)
    return (ξ * (xs * (γ_cr - 1.0) - (γ_cr + 1.0)) * xs^γ_th - xs * (γ_th + 1.0) + (γ_th - 1.0)) /
           (ξ * (xs * (γ_cr - 1.0) - (γ_cr + 1.0)) + xs * (γ_th - 1.0) - (γ_th + 1.0))
end

"""
    𝒳_f(xs::Real, ξ::Real)

Helper function.
"""
function 𝒳_f(xs::Real, ξ::Real)
    ys = ys_f(xs, ξ)
    return (γ_cr - 1.0) / (γ_th - 1.0) * ξ * (ys - xs^γ_th)
end

"""
    get_xs(M)

Shock compression ratio in the ideal hydrodynamical case.
"""
function get_xs(M)
    (γ_th + 1) * M^2 / (2 + (γ_th - 1) * M^2)
end

"""
    X_cr(Mach::Real, ξ::Real)

Ratio between CR and thermal pressure for an injected CR spectrum following Pfrommer et. al. (2017).
Uses the assumption that the CRs do not significantly alter the shock properties, which is not valid in the case of efficient CR proton injection!
"""
function X_cr(Mach::Real, ξ::Real)
    # shock compression ratio
    xs = get_xs(Mach)
    # helper functions
    𝒳 = 𝒳_f(xs, ξ)
    ys = ys_f(xs, ξ)
    return 𝒳 / (ys + 𝒳)
end


"""
    integrate_θ(x_in::Real, synch_F::Function, θ_steps::Integer=50)

Pitch angle integration in Donnert+16, Eq. 17.
`synch_F` can be either the first or second synchrotron function.
"""
function integrate_θ(x_in::Real, synch_F::Function, θ_steps::Integer=50)

    dθ = 0.5π / θ_steps
    θ = 0.0

    # first half step: sin(0) = 0
    K = 0.0

    # actual integration
    @inbounds for i ∈ 1:θ_steps-1
        θ += dθ
        sinθ = sin(θ)
        x = x_in / sinθ
        K += sinθ^2 * synch_F(x)
    end

    # last step: sin(0.5π) = 1
    K += 0.5 * synch_F(x_in)

    # multiply by step length
    K *= dθ

    return K
end


"""
    solve_synchrotron_function( x::Real, synch_F::Function, 
                                integrate_pitch_angle::Bool)

Helper function for solving the pitch angle integration.
"""
function solve_synchrotron_function(x::Real, synch_F::Function,
    integrate_pitch_angle::Bool)
    if integrate_pitch_angle
        return integrate_θ(x, synch_F)
    else
        return synch_F(x)
    end
end

"""

"""
function get_synch_emissivity_integral(ϵ_th::Real, Mach::Real, B::Real, etaB::Real; 
                                        spectrum::Function, K_ep::Real, ν0::Real, 
                                        η_model::AbstractShockAccelerationEfficiency,
                                        synch_F::Function,
                                        integrate_pitch_angle::Bool,
                                        N_steps=128)

    # combined acceleration efficiency coefficient
    ξ = etaB * η_Ms(η_model, Mach, 0.0) / (1.0 - etaB * η_Ms(η_model, Mach, 0.0))

    # no need to calculate emission if efficiency or magnetic field is zero!
    if iszero(ξ) || iszero(B)
        return 0.0
    end

    # energy density of CR electrons in [erg/cm^3]
    ϵ_cr = K_ep * X_cr(Mach, ξ) * ϵ_th

    # allocate storage arrays 
    x = Vector{Float64}(undef, N_steps)
    E = Vector{Float64}(undef, N_steps)
    dE = Vector{Float64}(undef, N_steps)
    F = Vector{Float64}(undef, N_steps)
    Fmid = Vector{Float64}(undef, N_steps)

    # store emissivity
    j_nu = 0.0

    # limits for energy integration
    E_min = sqrt(ν0 / (C_crit * B * 70.0))
    E_max = sqrt(ν0 / (C_crit * B * 1e-20))

    # step size in log space
    di = log(E_max / E_min) / (N_steps - 1)

    # initial step by itself
    x[1] = 70.0
    E[1] = E_min
    dE[1] = E[1] - E_min * exp(-1 * di)

    @inbounds for i = 2:N_steps
        E[i] = E_min * exp(di * i)
        dE[i] = E[i] - E[i-1]
        x[i] = ν0 / (C_crit * E[i]^2 * B)
    end

    # set initial value to 0
    F[1] = 0.0

    # integrate over energy space
    @inbounds for i = 2:N_steps
        # beginning of bin
        N_E = ϵ_cr * spectrum(E[i])
        F[i] = N_E * solve_synchrotron_function(x[i], synch_F, integrate_pitch_angle)

        # middle of bin 
        x_mid = 0.5 * (x[i-1] + x[i])
        e_mid = 0.5 * (E[i-1] + E[i])

        N_E_mid = ϵ_cr * spectrum(e_mid)
        Fmid[i] = N_E_mid * solve_synchrotron_function(x_mid, synch_F, integrate_pitch_angle)

        # simpson rule 
        j_nu += dE[i] / 6 * (F[i] + F[i-1] + 4 * Fmid[i])
    end

    return j_nu_prefac * B * j_nu
end

"""
    analytic_synchrotron(P_cgs::Array{<:Real}, B_cgs::Array{<:Real}, 
                         Mach::Array{<:Real}, θ_B::Union{Nothing,Array{<:Real}}=nothing;
                         dsa_model::Union{Integer,AbstractShockAccelerationEfficiency}=1, 
                         ν0::Real=1.4e9,
                         K_ep::Real=0.01, CR_Emin::Real=1.0,
                         spectrum::Union{Nothing,Function}=nothing,
                         integrate_pitch_angle::Bool=true,
                         polarisation::Bool=false,
                         show_progress::Bool=false)

Computes the analytic synchrotron emission from a spectrum of electrons by explicitly integrating over the distribution function.
The integral over the spectrum must be normalized to 1.
The total energy density of the relativistic electrons is given by the CR to thermal pressure ratio obtained by employing a DSA model and computing Xcr as in [Pfrommer et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.4500P/abstract).

Returns synchrotron emissivity `j_nu` in units `[erg/s/Hz/cm^3]`.

# Arguments
- `P_cgs::Array{<:Real}`:   Thermal energy density in `erg/cm^3`.
- `B_cgs::Array{<:Real}`:   Magnetic field in Gauss.
- `Mach::Array{<:Real}`:    Mach number.
- `θ_B::Union{Nothing,Array{<:Real}}=nothing`: Shock obliquity (optional).

## Keyword Arguments
- `ν0::Real=1.4e9`:           Observation frequency in `Hz`.
- `dsa_model`:      Diffusive Shock Acceleration model. Takes values `0...4`, or custom model. See next section.
- `K_ep::Real=0.01`:           Ratio of CR proton to electron energy acceleration.
- `CR_Emin::Real=1`:           Injection energy of CR electron population in `GeV`.
- `spectrum::Union{Nothing,Function}=nothing`: Spectrum function. Must be normalized so that the integral over it is 1.
- `integrate_pitch_angle::Bool=true`: Optional avoid pitch angle integration to reduce computational cost.
- `polarisation::Bool=false`: Set to `true` if you want to compute the polarized emission instead of the total intensity.
- `show_progress::Bool=false`: Enables a progress bar if set to true.

## DSA Models
Takes either your self-defined `AbstractShockAccelerationEfficiency` (see [DiffusiveShockAccelerationModels.jl](https://github.com/LudwigBoess/DiffusiveShockAccelerationModels.jl) for details!)
or a numerical value as input.
Numerical values correspond to:
- `0`: [Kang et. al. (2007)](https://ui.adsabs.harvard.edu/abs/2007ApJ...669..729K/abstract)
- `1`: [Kang & Ryu (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...95K/abstract)
- `2`: [Ryu et. al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...883...60R/abstract)
- `3`: [Caprioli & Spitkovsky (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...783...91C/abstract)
- `4`: [Pfrommer et. al. (2006)](https://ui.adsabs.harvard.edu/abs/2006MNRAS.367..113P/abstract)

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function analytic_synchrotron(P_cgs::Array{<:Real}, B_cgs::Array{<:Real}, 
                              Mach::Array{<:Real}, θ_B::Union{Nothing,Array{<:Real}}=nothing;
                              ν0::Real=1.4e9,
                              dsa_model::Union{Integer,AbstractShockAccelerationEfficiency}=1, 
                              K_ep::Real=0.01, CR_Emin::Real=1.0,
                              spectrum::Union{Nothing,Function}=nothing,
                              integrate_pitch_angle::Bool=true,
                              polarisation::Bool=false,
                              show_progress::Bool=false)

    # check if intensity or polarisation should be computed
    if !polarisation
        synch_F = ℱ
    else
        synch_F = 𝒢
    end

    # assign requested DSA model
    η_model = select_dsa_model(dsa_model)

    # allocate storage array
    j_nu = Vector{Float64}(undef, length(P_cgs))

    # enable progress meter
    if show_progress
        p = Progress(length(P_cgs))
    end

    # loop over all particles
    @threads for i = 1:length(P_cgs)

        # compute shock obliquity dependent acceleration efficiency if requested
        if !isnothing(θ_B)
            ηB = ηB_acc_e(θ_B[i])
        else
            ηB = 1.0
        end

        η_tot = η_Ms_acc(η_model, Mach[i]) * ηB 

        if η_tot > 0.0

            # default to powerlaw spectrum if no user input is set
            if isnothing(spectrum)
                # slope of injected spectrum
                s = dsa_spectral_index(Mach[i])

                # default to simple powerlaw spectrum
                spectrum(E::Real) = powerlaw_spectrum(E, s, CR_Emin * eV2cgs*1.e9)
            end

            # get emissivity in [erg/s/Hz/cm^3]
            j_nu[i] = get_synch_emissivity_integral(P_cgs[i], Mach[i], B_cgs[i], ηB; 
                                                    spectrum, K_ep, η_model, ν0, 
                                                    synch_F, integrate_pitch_angle)

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

