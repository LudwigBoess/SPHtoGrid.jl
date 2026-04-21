
# =====================================================================
# 2D Power Spectrum Implementation
# =====================================================================

"""
    power_spectrum_2D(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; cells_per_Mpc = 1.0)

Computes the 2D power spectrum of the given fields. The fields are expected to be 3D arrays with dimensions (nx, ny, nz) or 2D arrays with dimensions (nx, ny). 
The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum_2D(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; 
                            cells_per_Mpc::Real = 1.0, verbose::Bool = false)
    nx, ny = size(field_x, 1), size(field_x, 2)

    verbose && @info "Calculating 2D power spectrum with dimensions: ($nx, $ny)"

    verbose && @info "Calculating frequencies for 2D grid..."
    # Calculate frequencies
    kx = fftfreq(nx, cells_per_Mpc)
    ky = rfftfreq(ny, cells_per_Mpc)
    kmin = ky[2] 
    

    verbose && @info "Performing 2D RFFT on fields..."
    # Perform 2D real fourier transform 
    # (halves the second dimension as specified by the tuple (2,1))
    field_x_fft = rfft(field_x, (2, 1))
    field_y_fft = rfft(field_y, (2, 1))
    field_z_fft = rfft(field_z, (2, 1))
    
    nx_fft, ny_fft = size(field_x_fft, 1), size(field_x_fft, 2)
    
    # Determine maximum k to setup bins correctly
    kmax = sqrt(maximum(abs, kx)^2 + ky[end]^2)
    k_bins = collect(0:kmin:(kmax - kmin/10))
    nbins = length(k_bins)
    
    amp = zeros(nbins)

    P = Progress(ny_fft, desc="Binning Fourier modes...")
    
    # Loop over Fourier modes and bin by kperp
    for j in 1:ny_fft
        # Multiply energy by 2 for internal RFFT frequencies
        nf = (1 < j < ny_fft) ? 2.0 : 1.0
        ky_val = ky[j]
        
        for i in 1:nx_fft
            kx_val = kx[i]
            kperp = sqrt(kx_val^2 + ky_val^2)
            
            # Compute bin index
            bin_idx = floor(Int, kperp / kmin) + 1
            
            if 1 <= bin_idx <= nbins
                # Calculate energy on the fly
                energy = (abs2(field_x_fft[i, j]) + 
                          abs2(field_y_fft[i, j]) + 
                          abs2(field_z_fft[i, j])) * nf
                
                amp[bin_idx] += energy
            end
        end

        verbose && next!(P)
    end

    return k_bins, amp
end

"""
    power_spectrum_2D(field::AbstractArray; cells_per_Mpc = 1.0)
Computes the 2D power spectrum of a single scalar field. The field is expected to be a 2D array with dimensions (nx, ny). The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum_2D(field::AbstractArray; cells_per_Mpc::Real = 1.0, verbose::Bool = false)
    nx, ny = size(field)
    
    verbose && @info "Calculating 2D power spectrum with dimensions: ($nx, $ny)"
        
    verbose && @info "Calculating frequencies for 2D grid..."
    # Frequencies: rfft halves the dimension listed first in the region tuple.
    # Below, rfft(field, (2, 1)) will halve dimension 2 (ny).
    kx = fftfreq(nx, cells_per_Mpc)
    ky = rfftfreq(ny, cells_per_Mpc)
    kmin = ky[2] 
    
    verbose && @info "Performing 2D RFFT on field..."
    # 2D RFFT 
    field_fft = rfft(field, (2, 1))
    nx_fft, ny_fft = size(field_fft)
    
    # Setup bins
    kmax = sqrt(maximum(abs, kx)^2 + ky[end]^2)
    k_bins = collect(0:kmin:(kmax - kmin/10))
    nbins = length(k_bins)
    
    amp = zeros(nbins)

    P = Progress(ny_fft, desc="Binning Fourier modes...")

    # Single-pass loop fusion (Column-major iteration)
    @inbounds for j in 1:ny_fft
        # Reality condition: double energy for internal frequencies
        nf = (1 < j < ny_fft) ? 2.0 : 1.0
        ky_val = ky[j]
        
        for i in 1:nx_fft
            kx_val = kx[i]
            kperp = sqrt(kx_val^2 + ky_val^2)
            
            bin_idx = floor(Int, kperp / kmin) + 1
            
            if 1 <= bin_idx <= nbins
                energy = abs2(field_fft[i, j]) * nf
                amp[bin_idx] += energy
            end
        end

        verbose && next!(P)
    end
    
    return k_bins, amp
end

# =====================================================================
# 3D Power Spectrum Implementation
# =====================================================================


"""
    power_spectrum_3D(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; k_norm = 1.0)

Computes the 3D power spectrum of the given fields. The fields are expected to be 3D arrays with dimensions (nx, ny, nz).
The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum_3D(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; 
                            cells_per_Mpc::Real = 1.0, verbose::Bool = false)
    
    nx, ny, nz = size(field_x)
    verbose && @info "Calculating 3D power spectrum with dimensions: ($nx, $ny, $nz)"
    
    verbose && @info "Calculating frequencies for 3D grid..."
    # Calculate frequencies for a 3D grid
    kx = rfftfreq(nx, cells_per_Mpc)
    ky = fftfreq(ny, cells_per_Mpc)
    kz = fftfreq(nz, cells_per_Mpc)
    
    # Use the fundamental wavenumber for bin widths
    kmin = kx[2] 
    
    verbose && @info "Performing 3D RFFT on fields..."
    # Perform 3D RFFT (defaults to transforming all 3 dimensions)
    field_x_fft = rfft(field_x)
    field_y_fft = rfft(field_y)
    field_z_fft = rfft(field_z)
    
    nx_fft, ny_fft, nz_fft = size(field_x_fft)
    
    # Determine maximum k (the furthest corner of the 3D k-space box)
    kmax = sqrt(kx[end]^2 + maximum(abs, ky)^2 + maximum(abs, kz)^2)
    k_bins = collect(0:kmin:(kmax - kmin/10))
    nbins = length(k_bins)
    
    amp = zeros(nbins)
    
    # Precompute whether nx is even for the reality condition
    nx_is_even = (nx % 2 == 0)

    P = Progress(nz_fft, desc="Binning Fourier modes...")
    
    # Single-pass loop fusion for 3D spherical shells
    # Loop order corresponds to memory layout: innermost loop processes the 1st dimension.
    @inbounds for k in 1:nz_fft
        kz_val = kz[k]
        
        for j in 1:ny_fft
            ky_val = ky[j]
            
            for i in 1:nx_fft
                kx_val = kx[i]
                
                # 3D Wavenumber magnitude
                k_mag = sqrt(kx_val^2 + ky_val^2 + kz_val^2)
                
                # Compute bin index mathematically
                bin_idx = floor(Int, k_mag / kmin) + 1
                
                if 1 <= bin_idx <= nbins
                    # Reality condition for rfft over the FIRST dimension:
                    # DC (i=1) and Nyquist (if nx is even and i=nx_fft) are unique.
                    # All other modes represent two conjugate frequencies, so we double their energy.
                    nf = (i == 1 || (nx_is_even && i == nx_fft)) ? 1.0 : 2.0
                    
                    # Calculate energy on the fly
                    energy = (abs2(field_x_fft[i, j, k]) + 
                              abs2(field_y_fft[i, j, k]) + 
                              abs2(field_z_fft[i, j, k])) * nf
                    
                    amp[bin_idx] += energy
                end
            end
        end

        verbose && next!(P)
    end
    
    return k_bins, amp
end


"""
    power_spectrum_3D(field::AbstractArray; cells_per_Mpc = 1.0)
Computes the 3D power spectrum of a single scalar field. 
The field is expected to be a 3D array with dimensions (nx, ny, nz). 
The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum_3D(field::AbstractArray; cells_per_Mpc::Real = 1.0, verbose::Bool = false)
    nx, ny, nz = size(field)

    verbose && @info "Calculating 3D power spectrum with dimensions: ($nx, $ny, $nz)"
    
    verbose && @info "Calculating frequencies for 3D grid..."
    # Frequencies: Default 3D rfft() halves the FIRST dimension (nx)
    kx = rfftfreq(nx, cells_per_Mpc)
    ky = fftfreq(ny, cells_per_Mpc)
    kz = fftfreq(nz, cells_per_Mpc)
    kmin = kx[2] 
    
    verbose && @info "Performing 3D RFFT on field..."
    # 3D RFFT
    field_fft = rfft(field)
    nx_fft, ny_fft, nz_fft = size(field_fft)
    
    # Setup bins
    kmax = sqrt(kx[end]^2 + maximum(abs, ky)^2 + maximum(abs, kz)^2)
    k_bins = collect(0:kmin:(kmax - kmin/10))
    nbins = length(k_bins)
    
    amp = zeros(nbins)
    nx_is_even = (nx % 2 == 0)
    
    P = Progress(nz_fft, desc="Binning Fourier modes...")
    # Single-pass loop fusion (Column-major iteration)
    @inbounds for k in 1:nz_fft
        kz_val = kz[k]
        for j in 1:ny_fft
            ky_val = ky[j]
            for i in 1:nx_fft
                kx_val = kx[i]
                
                k_mag = sqrt(kx_val^2 + ky_val^2 + kz_val^2)
                bin_idx = floor(Int, k_mag / kmin) + 1
                
                if 1 <= bin_idx <= nbins
                    # Reality condition for rfft over the FIRST dimension
                    nf = (i == 1 || (nx_is_even && i == nx_fft)) ? 1.0 : 2.0
                    
                    energy = abs2(field_fft[i, j, k]) * nf
                    amp[bin_idx] += energy
                end
            end
        end
        verbose && next!(P)
    end
    
    return k_bins, amp
end


# =====================================================================
# Main functions that dispatch to 2D or 3D implementations based on input dimensions
# =====================================================================


"""
    power_spectrum(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; cells_per_Mpc = 1.0)

Computes the power spectrum of the given fields. 
The fields are expected to be 3D arrays with dimensions (nx, ny, nz) or 2D arrays with dimensions (nx, ny). 
The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum(field_x::AbstractArray, field_y::AbstractArray, field_z::AbstractArray; 
                        cells_per_Mpc = 1.0, verbose::Bool = false)
    if ndims(field_x) == 2
        return power_spectrum_2D(field_x, field_y, field_z; cells_per_Mpc, verbose)
    elseif ndims(field_x) == 3
        return power_spectrum_3D(field_x, field_y, field_z; cells_per_Mpc, verbose)
    else
        error("Unsupported number of dimensions: $(ndims(field_x)). Only 2D and 3D arrays are supported.")
    end 
end

"""
    power_spectrum(field::AbstractArray; cells_per_Mpc = 1.0)

Computes the power spectrum of a single scalar field. 
The field is expected to be a 3D array with dimensions (nx, ny, nz) 
or a 2D array with dimensions (nx, ny). 
The `cells_per_Mpc` parameter is used to scale the frequencies appropriately.
"""
function power_spectrum(field::AbstractArray; cells_per_Mpc = 1.0, verbose::Bool = false)
    if ndims(field) == 2
        return power_spectrum_2D(field; cells_per_Mpc, verbose)
    elseif ndims(field) == 3
        return power_spectrum_3D(field; cells_per_Mpc, verbose)
    else
        error("Unsupported number of dimensions: $(ndims(field)). Only 2D and 3D arrays are supported.")
    end 
end
