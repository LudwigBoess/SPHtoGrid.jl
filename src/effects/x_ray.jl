"""
    interpolate_table(Qi, Q_table)

Finds the first position in the provided table where `Qi` is larger than a table entry.
"""
function interpolate_table(Qi, Q_table)

    sel = findfirst(Q_table .>= Qi)
    if isnothing(sel)
        sel = length(Q_table)
    end

    return sel
end


"""
    get_cooling_emissivity(T_keV, rho_cgs, metalicity; E0, E1)

Computes the Xray emissivity from interpolating a cooling table.
"""
function get_cooling_emissivity(T_keV, rho_cgs, metalicity; E0, E1)

    # read all cooling tables
    Temp, Zmetal, ne_nH, mue, dLcool, E_low, E_high = read_cooling_tables()

    # get energy band indices
    iE_low  = interpolate_table(E0, E_low)
    iE_high = interpolate_table(E1, E_high)

    # store number of particles
    Npart = length(T_keV)

    # allocate storage arrays
    q_x = Vector{Float64}(undef, Npart)

    # fill in factors
    @threads for i = 1:Npart

        # interpolate in the tables 
        if !isnothing(metalicity)
            iMetal = interpolate_table(metalicity[i], Zmetal)
        else 
            iMetal = 1
        end

        iTemp = interpolate_table(T_keV[i], Temp)

        Lumfact = ne_nH[iMetal] * (mue[iMetal] * m_p)^2

        Lcool_band = 0.0
        for iE = iE_low:iE_high
            Lcool_band += dLcool[iMetal, iTemp, iE]
        end

        q_x[i] = rho_cgs[i]^2 * 1.e-23 * Lcool_band / Lumfact
    end

    return q_x
end

"""
    x_ray_emissivity(T_keV::Vector{<:Real}, 
                     rho_cgs::Vector{<:Real},
                     metalicity::Union{Vector{Float64, Nothing}}=nothing; 
                     E0::Real=0.1, E1::Real=2.4, 
                     xH::Real=0.752,
                     cooling_function::Bool=false,
                     z::Real=0.0)

X-Ray emissivity for particles with temperature `T_keV` in ``keV``, and density `rho_cgs` in ``g/cm^3``.
If available you can also add the `metalicity` in the gas.
`Emin` and `Emax` give the minimum and maximum energy of the oberservation.
`xH` gives the hydrogen fraction used in the simulation.

# Returns
X-Ray emissivity in units of [erg/s/cm^3].

## Arguments:
- `T_keV`: SPH particle temperature [keV]
- `m_cgs`: SPH particle mass in [g]
- `rho_cgs`: SPH particle density in [g/cm^3]
- `E0`: Minimum photon energy for Xray spectrum [keV]
- `E1`: Maximum photon energy for Xray spectrum [keV]
- `xH`: Hydrogen mass fraction in the simulation

## Mapping settings
- weight function: [`part_weight_physical`](@ref)
- reduce image: `false`
"""
function x_ray_emissivity(T_keV::Vector{<:Real}, 
                        rho_cgs::Vector{<:Real},
                        metalicity::Union{Vector{Float64}, Nothing}=nothing; 
                        E0::Real=0.1, E1::Real=2.4, 
                        xH::Real=0.752,
                        cooling_function::Bool=false,
                        z::Real=0.0,
                        smac_filter::Bool=false)

    # shift from observer to rest-frame
    E0 *= (1 + z)
    E1 *= (1 + z)

    # apply cooling function if requested
    if cooling_function
        q_x = get_cooling_emissivity(T_keV, rho_cgs, metalicity; E0, E1)
    else

        mol  = 4 / (5 * xH + 3);
        n2ne = (xH + 0.5 * (1 - xH)) /
            (2xH + 0.75 * (1 - xH));
        
        xray_prefactor = 4C_j * gg * (1 + xH)^2 * (n2ne / (mol * m_p))^2
        
        """
            Steinmetz & Bartelmann, based on Spizer 1968, gg = 1.2 (!?)
        Beside the fact, that it is not clear which value they used for the Gaunt factor
        it is the best formulation, as composition H/He (fr) and conversion from particle
        number to electron number (n2ne) is explicite formulated.
        """
        q_x = Vector{Float64}(undef, length(rho_cgs))

        @threads for i ∈ eachindex(q_x)
            cutoff = exp(-E0 / T_keV[i]) - exp(-E1 / T_keV[i])
            q_x[i] = xray_prefactor * cutoff * rho_cgs[i]^2 * √(T_keV[i])
        end
    end

    # apply smac filter if requested
    if smac_filter
        apply_smac_T_rho_filter!(q_x, T_keV .* keV2K, rho_cgs)
    end

    return q_x
end
