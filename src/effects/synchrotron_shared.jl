"""
    select_dsa_model(dsa_model::Integer)

Helper function to select the DSA model based on integer values.
"""
function select_dsa_model(dsa_model::Integer)

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

    return η_model
end

"""
    select_dsa_model(dsa_model::AbstractShockAccelerationEfficiency)

Helper function to assign the self-defined DSA model.
"""
select_dsa_model(dsa_model::AbstractShockAccelerationEfficiency) = dsa_model


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
function EpsNtherm(rho_cgs::Real, T_K::Real; xH::Real=0.76)
    u_mol = (4 / (5xH + 3))
    return (rho_cgs / (m_p * u_mol) * k_B * T_K)
end

