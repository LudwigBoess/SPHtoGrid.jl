global const ev2K = 1.160451812e4
global const K2eV = 1.0/1.160451812e4

"""
    part_weight_one(N::Integer)

Equivalent to no weighting. Returns an Array of ones.
"""
function part_weight_one(N::Integer)
    return ones(N)
end

"""
    part_weight_physical(N::Integer, par::mappingParameters)

Physical weighting function.
"""
function part_weight_physical(N::Integer, par::mappingParameters)
    return ones(N) .* par.pixelSideLength
end

"""
    part_weight_emission(rho::Array{<:Real}, T::Array{<:Real})

Emission weighted mapping. Takes density and temperature and computes weights.
"""
function part_weight_emission(rho::Array{<:Real}, T::Array{<:Real})
    return rho.^2 .* √.(T)
end


"""
    part_weight_spectroscopic(rho::Array{<:Real}, T::Array{<:Real})

Spectroscopic weighted mapping from Mazotta+ 04. Takes density and temperature and computes weights.
"""
function part_weight_spectroscopic(rho::Array{<:Real}, T::Array{<:Real})
    return rho.^2 .* T.^(0.75 - 1.5)
end

"""
    part_weight_XrayBand(T::Array{<:Real}, Emin::Real, Emax::Real)

Computes Xray weighted emission of a defined energy band. Emin and Emax are energies in eV.
"""
function part_weight_XrayBand(T::Array{<:Real}, Emin::Real, Emax::Real)
    # convert Kelvin to eV
    T_eV = T .* k2eV

    @. exp( -Emin / ( kB * T_eV )) - exp( -Emax / ( kB * T_eV ))
end