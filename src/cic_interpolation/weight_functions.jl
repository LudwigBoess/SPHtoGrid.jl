global const ev2K = 1.160451812e4
global const K2eV = 1.0/1.160451812e4
global const kB = 1.38066e-16

"""
    part_weight_one(N::Integer)

Equivalent to no weighting. Returns an Array of ones.
"""
function part_weight_one(N::Integer)
    return ones(N)
end

"""
    part_weight_physical(N::Integer, par::mappingParameters)

Physical weighting function in units of [cm/pix].
"""
function part_weight_physical(N::Integer, par::mappingParameters, x_cgs::Real=3.085678e21)
    return ones(N) .* par.pixelSideLength .* x_cgs
end

"""
    part_weight_emission(rho::Array{<:Real}, T_K::Array{<:Real})

Emission weighted mapping. Takes density in internal untis and temperature in K and computes weights.
"""
function part_weight_emission(rho::Array{<:Real}, T_K::Array{<:Real})
    return @. rho^2 * √(T_K)
end


"""
    part_weight_spectroscopic(rho::Array{<:Real}, T_K::Array{<:Real})

Spectroscopic weighted mapping from Mazotta+ 04. Takes density and temperature and computes weights.
"""
function part_weight_spectroscopic(rho::Array{<:Real}, T_K::Array{<:Real})
    return @. rho^2 / √(√(T_K))^3
end

"""
    part_weight_XrayBand(T_K::Array{<:Real}, Emin::Real, Emax::Real)

Computes Xray weighted emission of a defined energy band. Emin and Emax are energies in eV.
"""
function part_weight_XrayBand(T_K::Array{<:Real}, Emin::Real=5.0e4, Emax::Real=1.0e10)
    # convert Kelvin to eV
    T_eV = T_K .* cgs2eV

    @. exp( -Emin / ( kB * T_eV )) - exp( -Emax / ( kB * T_eV ))
end