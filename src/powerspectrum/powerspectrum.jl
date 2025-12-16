"""
    flatten_powerspec(dim::Int, data; verbose::Bool=false)

Computes a 1D powerspectrum for a grid of arbitrary dimensions.
"""
flatten_powerspec(dim::Int, data; verbose::Bool=false) = flatten_powerspec(Val(dim), data; verbose)

"""
    flatten_powerspec(::Val{2}, fft_field; verbose::Bool=false)

Computes a 1D powerspectrum from a 1D fourier grid by computing the average in shells.
"""
flatten_powerspec(::Val{1}, spectrum; verbose::Bool=false) = ifftshift(spectrum)

"""
    flatten_powerspec(::Val{2}, fft_field; verbose::Bool=false)

Computes a 1D powerspectrum from a 2D fourier grid by computing the average in shells.
"""
function flatten_powerspec(::Val{2}, fft_field; verbose::Bool=false)
    N = size(fft_field, 1)
    center = N ÷ 2 + 1
    
    # Maximum possible radius is corner of the box, but usually we just go to Nyquist (N/2)
    max_r = floor(Int, N / 2)
    
    # Arrays to store sum of magnitudes and count of pixels per radius bin
    radial_sum = zeros(Float64, max_r + 1)
    radial_count = zeros(Int, max_r + 1)
    
    if verbose
        println("Computing radial average (1D spectrum)...")
        p = Progress(N^2)
    end
    
    @inbounds for j in 1:N, i in 1:N
        # Calculate distance from center (DC component)
        dx = i - center
        dy = j - center
        dist = sqrt(dx^2 + dy^2)
        
        # Determine integer bin index (0 to max_r)
        bin_idx = round(Int, dist)
        
        # Only accumulate if within our interest range
        if bin_idx <= max_r
            # Julia arrays are 1-based, so map bin 0 -> index 1
            radial_sum[bin_idx + 1] += fft_field[i, j]
            radial_count[bin_idx + 1] += 1
        end

        if verbose
            # count up progress meter
            next!(p)
        end
    end
    
    # avoid division by zero
    radial_avg = zeros(Float64, length(radial_sum))
    @inbounds for i in 1:length(radial_sum)
        if radial_count[i] > 0
            radial_avg[i] = radial_sum[i] / radial_count[i]
        end
    end
    
    return radial_avg
end

"""
    flatten_powerspec(::Val{3}, fft_field; verbose::Bool=false)

Computes a 1D powerspectrum from a 3D fourier grid by computing the average in shells.
"""
function flatten_powerspec(::Val{3}, fft_field; verbose::Bool=false)
    N = size(fft_field, 1)
    center = N ÷ 2 + 1
    
    # Maximum possible radius is corner of the box, but usually we just go to Nyquist (N/2)
    max_r = floor(Int, N / 2)
    
    # Arrays to store sum of magnitudes and count of pixels per radius bin
    radial_sum = zeros(Float64, max_r + 1)
    radial_count = zeros(Int, max_r + 1)
    
    if verbose
        println("Computing radial average (1D spectrum)...")
        p = Progress(N)
    end
    
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        # Calculate distance from center (DC component)
        dx = i - center
        dy = j - center
        dz = k - center
        dist = sqrt(dx^2 + dy^2 + dz^2)
        
        # Determine integer bin index (0 to max_r)
        bin_idx = round(Int, dist)
        
        # Only accumulate if within our interest range
        if bin_idx <= max_r
            # Julia arrays are 1-based, so map bin 0 -> index 1
            radial_sum[bin_idx + 1] += fft_field[i, j, k]
            radial_count[bin_idx + 1] += 1
        end

        if verbose
            # count up progress meter
            next!(p)
        end
    end
    
    # avoid division by zero
    radial_avg = zeros(Float64, length(radial_sum))
    for i in 1:length(radial_sum)
        if radial_count[i] > 0
            radial_avg[i] = radial_sum[i] / radial_count[i]
        end
    end
    
    return radial_avg
end

"""
    power_spectrum(field::Array{<:Real}, box_width::Real; verbose::Bool=false)

Computes the power spectrum of a given real-valued field defined on a grid.
"""
function power_spectrum(field::Array{<:Real}, box_width::Real; verbose::Bool=false)

    if verbose
        @info "Fourier transform of grid"
    end
    # fourier transform of grid
    Fk = fftshift(fft(field))
    if verbose
        @info "Fourier transform done"
        @info "Computing powerspectrum"
    end
    # reduce from N-dim to 1-dim
    Pk = flatten_powerspec(ndims(Fk), abs2.(Fk); verbose)
    
    if verbose
        @info "Powerspectrum done"
    end
    # construct modes on grid
    N = size(Fk, 1)
    i_max = N ÷ 2
    k = [i/box_width for i = 0:i_max]

    return k, Pk
end