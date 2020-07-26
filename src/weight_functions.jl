global const ev2K = 1.160451812e4
global const K2eV = 1.0/1.160451812e4

function part_weight_one(N::Integer)
    return ones(N)
end

function part_weight_physical(N::Integer, par::mappingParameters)
    return ones(N) .* par.pixelSideLength
end

function part_weight_emission(rho::Array{<:Real}, T::Array{<:Real})
    return rho.^2 .* √.(T)
end

function part_weight_XrayBand(T::Array{<:Real}, Emin::Real, Emax::Real)
    # convert Kelvin to eV
    T_eV = T .* k2eV

    @. exp( -Emin / ( kB * T_eV )) - exp( -Emax / ( kB * T_eV ))
end