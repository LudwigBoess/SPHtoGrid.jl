function get_norm(pos::Vector{T}) where T

    Δx = zero(T)
    
    @inbounds for dim  = 1:3
        Δx += pos[dim]^2
    end

    return √(Δx)
end
