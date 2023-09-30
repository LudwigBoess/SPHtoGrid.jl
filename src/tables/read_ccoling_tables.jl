
"""
    Functions to read cooling tables
"""

using DelimitedFiles
using Downloads

function download_cooling_tables(tables_path)
    @info "Downloading band.data"
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/band.data", joinpath(tables_path, "band.data"))
    @info "Downloading dLambda.bin.data"
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/dLambda.bin.data", joinpath(tables_path, "dLambda.bin.data"))
    @info "Downloading T0_400.data"
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/T0_400.data", joinpath(tables_path, "T0_400.data"))
    @info "Downloading Z.data"
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/Z.data", joinpath(tables_path, "Z.data"))
    @info "Downloads complete!"
end

function read_metalicity_table(tables_path)

    # read the data
    z_data = readdlm(joinpath(tables_path, "Z.data"))

    # assign to arrays
    Zmetal = z_data[:, 1]
    ne_nH  = z_data[:, 2]
    mue    = z_data[:, 5]

    return Zmetal, ne_nH, mue
end

function read_temperature_table(tables_path)
    return readdlm(joinpath(tables_path, "T0_400.data"))[:, 1]
end

function read_energy_band_table(tables_path)
    data = readdlm(joinpath(tables_path, "band.data")) .* 1.e-3
    return data[:, 1], data[:, 2]
end

function read_dLcool(tables_path)
    f = open(joinpath(tables_path, "dLambda.bin.data"), "r")
    dummy = read(f, Int32)
    dLcool = read!(f, Array{Float64}(undef, 200, 400, 400))
    close(f)
    return dLcool
end

"""
    read_cooling_tables()

Reads the cooling tables into Arrays.
"""
function read_cooling_tables()

    if !isfile(joinpath(tables_path, "T0_400.data"))
        @info "Downloading cooling tables"
        download_cooling_tables(tables_path)
    end

    # read temperature table for interpolation
    Temp = read_temperature_table(tables_path)

    # Read Lambda.bin.data
    dLcool = read_dLcool(tables_path)

    # read metalicty table 
    Zmetal, ne_nH, mue = read_metalicity_table(tables_path)

    # read energy band table 
    E_low, E_high = read_energy_band_table(tables_path)

    return Temp, Zmetal, ne_nH, mue, dLcool, E_low, E_high
end