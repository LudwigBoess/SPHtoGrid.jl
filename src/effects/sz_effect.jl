"""
    Tcmb(z::Real)

Computes the temperature of the CMB at redshift `z`.
"""
Tcmb(z::Real) = (2.728 * (1 + z))

"""
    kSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

Prefactor for the kinetic Sunyaev-Zel'dovich effect.
"""
function kSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

    kSzPrefac = -1 * σ_T / c_light

    if DI_over_I
        x = h_planck * ν / (k_B * Tcmb(z))

        kSzPrefac *= exp(x) - 1 / (x * exp(x))
    end

    return kSzPrefac
end

"""
    kinetic_SZ(n_cm3::Real, vel_y_cgs::Real, 
               z::Real=0.0, ν::Real=1.e9;
               DI_over_I::Bool=false)

Computes the kinetic Sunyaev-Zel'dovich effect from electron density `n_cm3` and velocity in y-direction to the projection plane in cgs units `vel_y_cgs`.
If `DI_over_I` is set to `true` you also need to provide an observation frequency `ν` and redshift `z`.

## Arguments:
- `n_cm3`: SPH particle density in [1/cm^3]
- `vel_y_cgs`: SPH particle velocity in y-direction in [cm/s]
- `z`: Redshift
- `ν`: Observing frequency

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function kinetic_SZ(n_cm3::Vector{<:Real}, vel_y_cgs::Vector{<:Real},
                    z::Real=0.0, ν::Real=1.e9;
                    DI_over_I::Bool = false)

    # calculate prefator once
    prefac = kSzPrefac(ν, z, DI_over_I)

    # allocate storage vector
    kin_SZ = Vector{eltype(n_cm3[1])}(undef, length(n_cm3))

    @threads for i = 1:length(kin_SZ)
        kin_SZ[i] = prefac * n_cm3[i] * vel_y_cgs[i]
    end

    return kin_SZ
end

"""
    comptonY(n_cm3::Real, T_K::Real, z::Real)

Computes the Compton-Y parameter from electron density `n_cm3` and temperature `T` in Kelvin at redshift `z`.

## Arguments:
- `n_cm3`: SPH particle density in [1/cm^3]
- `T_K`: SPH particle temperature [K]
- `z`: Redshift.

## Mapping settings
- weight function: [`part_weight_one`](@ref)
- reduce image: `false`
"""
function comptonY(n_cm3::Real, T_K::Real, z::Real)
    return yPrefac * n_cm3 * (T_K - Tcmb(z))
end

"""
    comptonY(n_cm3::Vector{<:Real}, T_K::Vector{<:Real}, z::Real)

Computes the Compton-Y parameter from electron density `n_cm3` and temperature `T` in Kelvin at redshift `z`.

## Arguments:
- `n_cm3`: SPH particle density in [1/cm^3]
- `T_K`: SPH particle temperature [K]
- `z`: Redshift

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function comptonY(n_cm3::Vector{<:Real}, T_K::Vector{<:Real}, z::Real)

    T_cmb = Tcmb(z)

    compton_y = Vector{Float64}(undef, length(T_K))

    @threads for i ∈ eachindex(T_K)
        compton_y[i] = yPrefac * n_cm3[i] * (T_K[i] - T_cmb)
    end

    return compton_y
end

"""
    tSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

Computes the prefactor for the thermal Sunyaev-Zel'dovich effect.
"""
function tSzPrefac(ν::Real, z::Real, DI_over_I::Bool)

    x = h_planck * ν / (k_B * Tcmb(z))
    tSzPrefac = (x * (exp(x) + 1) / (exp(x) - 1) - 4)

    if DI_over_I
        tSzPrefac *= exp(x) - 1 / (x * exp(x))
    end

    return tSzPrefac
end

"""
    thermal_SZ( n_cm3::Vector{<:Real}, T_K::Vector{<:Real},
                z::Real=0.0, ν::Real=1.44e9; 
                DI_over_I::Bool=false )

Computes the thermal Sunyaev-Zel'dovich effect for electron density `n_cm3` and temperature `T_K` in Kelvin at redshift `z` and observer frequency `ν`.
`DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

## Arguments:
- `n_cm3`: SPH particle density in [1/cm^3]
- `T_K`: SPH particle temperature [K]
- `z`: Redshift
- `ν`: Observing frequency

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function thermal_SZ(n_cm3::Vector{<:Real}, T_K::Vector{<:Real},
                    z::Real = 0.0, ν::Real = 1.4e9;
                    DI_over_I::Bool = false)

    # calculate prefator once
    prefac = tSzPrefac(ν, z, DI_over_I)

    # allocate storage vector
    th_SZ = Vector{eltype(T_K[1])}(undef, length(T_K))

    @threads for i = 1:length(T_K)
        th_SZ[i] = prefac * comptonY(n_cm3[i], T_K[i], z)
    end

    return th_SZ
end