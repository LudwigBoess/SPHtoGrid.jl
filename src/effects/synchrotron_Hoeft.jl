"""
    B_cmb(z::Real)

Magnetic field equivalent of the CMB for a given redshift `z` in `[G]`.
"""
B_cmb(z::Real) = 3.2e-6*(1 + z)^2

"""
    analytic_synchrotron_HB07(rho_cgs::Array{<:Real}, m_cgs::Array{<:Real}, hsml_cgs::Array{<:Real},
                              B_cgs::Array{<:Real}, T_keV::Array{<:Real}, Mach::Array{<:Real};
                              xH::Real=0.752, ν0::Real=1.4e9, z::Real=0.0,
                              dsa_model::Integer=1, ξe::Real = 0.01 )

Computes the analytic synchrotron emission with the simplified approach described in [Hoeft&Brüggen (2007)](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375...77H/abstract), following approach by [Wittor et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.4448W/abstract).

Returns synchrotron emissivity `j_nu` in units [erg/s/Hzcm^3].

## Arguments
- `rho_cgs::Array{<:Real}`:  Density in ``g/cm^3``.
- `m_cgs::Array{<:Real}`:    Mass in ``g``.
- `hsml_cgs::Array{<:Real}`: `HSML` block in ``cm``.
- `B_cgs::Array{<:Real}`:    Magnetic field in ``G``.
- `T_keV::Array{<:Real}`:    Temperature in ``keV``.
- `Mach::Array{<:Real}`:     Sonic Mach number.
- `θ_B::Union{Nothing,Array{<:Real}}=nothing`: Shock obliquity (optional).

## Keyword Arguments
- `xH::Float64 = 0.76`:        Hydrogen fraction of the simulation, if run without chemical model.
- `ν0::Real=1.44e9`:           Observation frequency in ``Hz``.
- `z::Real=0.0`:               Redshift of the simulation.
- `dsa_model::Integer=1`:      Diffusive Shock Acceleration model. Takes values `0...4`, see next section.
- `ξe::Real=0.01`:             Ratio of CR proton to electron energy acceleration.
- `show_progress::Bool=false`: Enables a progress bar if set to true


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
function analytic_synchrotron_HB07(rho_cgs::Array{<:Real}, m_cgs::Array{<:Real}, hsml_cgs::Array{<:Real},
                                    B_cgs::Array{<:Real}, T_keV::Array{<:Real}, Mach::Array{<:Real},
                                    θ_B::Union{Nothing, Array{<:Real}}=nothing;
                                    xH::Real=0.752, ν0::Real=1.4e9, z::Real=0.0,
                                    dsa_model::Union{Integer,AbstractShockAccelerationEfficiency}=1,
                                    ξe::Real=0.01,
                                    show_progress::Bool=false )

    # assign requested DSA model
    η_model = select_dsa_model(dsa_model)

    # prefactor to HB07, Eq. 32
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

    # enable progress meter
    if show_progress
        p = Progress(length(j_nu))
    end

    @threads for i = 1:length(rho_cgs)

        if !isnothing(θ_B)
            ηB = ηB_acc_e(θ_B[i])
        else
            ηB = 1.0
        end

        η_tot = η_Ms_acc(η_model, Mach[i]) * ηB 

        if η_tot > 0

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

            # Wittor+17, Eq. 16, but in untis erg/s/Hz/cm^3
            j_nu[i] = j_ν_prefactor * A * ne * ξe_factor * (ν0/1.4e9)^(-s/2) * 
                    √(T_keV[i] / 7)^3 * Bfactor / V *
                    η_tot

        else
            j_nu[i] = 0.0
        end

        # update progress meter
        if show_progress
            next!(p)
            flush(stdout)
            flush(stderr)
        end

        if isnan(j_nu[i])
            error("NaN!\nA = $A 
                    ne = $ne 
                    s = $s
                    ξe_factor = $ξe_factor
                    (ν0/1.4e9)^(-s/2) = $((ν0/1.4e9)^(-s/2))
                    T = $(T_keV[i]) 
                    √(T_keV[i] / 7)^3 = $(√(T_keV[i] / 7)^3)
                    Bfactor = $Bfactor 
                    V = $V
                    etaM = $(η_Ms_acc(η_model, Mach[i])) 
                    etaB = $ηB
                    j_ν_prefactor = $j_ν_prefactor
                    ")
        end
    end

    return j_nu
end

