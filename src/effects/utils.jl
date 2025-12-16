function apply_smac_T_rho_filter!(Q::Vector{<:Real},
                                 T_K::Vector{<:Real}, 
                                 rho_cgs::Vector{<:Real})

    @inbounds for i = 1:length(Q)
        if T_K[i] < 3.0e4 && rho_cgs[i] > 1.9e-28
            Q[i] = 0.0
        end
    end

    nothing
end

function get_smac_T_rho_filter(T_K::Vector{<:Real}, 
                               rho_cgs::Vector{<:Real})

    mask = falses(length(T_K))
    @inbounds for i = 1:length(T_K)
        if T_K[i] < 3.0e4 && rho_cgs[i] > 1.9e-28
            mask[i] = true
        end
    end
    return mask
end