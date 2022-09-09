module SPHtoGrid

using Distributed
using SPHKernels
using Unitful
using Printf
using ProgressMeter
using DSAModels
using SpecialFunctions
using LinearAlgebra
using Rotations
using Base.Threads
using Healpix

function output_time(t1, t2)
    return @sprintf("%0.3e", Float64((t2 - t1)) * 1.e-9)
end

# shared functionality
include("shared/parameters.jl")
include("shared/reconstruct_grid.jl")
include("shared/indices.jl")
include("shared/distances.jl")
include("shared/periodic_mapping.jl")
include("shared/rotate_particles.jl")
include("shared/io.jl")

# multi-core functionality 
include("parallel/domain_decomp.jl")

# functions for smac
include("smac/smac1_utility.jl")
include("smac/smac2_utility.jl")

# tsc interpolation
include("tsc_interpolation/tsc_interpolation.jl")

# cic interpolation 
include("cic_interpolation/filter_shift.jl")
include("cic_interpolation/reduce_image.jl")
include("cic_interpolation/cic_shared.jl")
include("cic_interpolation/cic_2D.jl")
include("cic_interpolation/cic_3D.jl")
include("cic_interpolation/weight_functions.jl")
include("cic_interpolation/cic_interpolation.jl")

# splash interpolation
include("splash_interpolation/reduce_image.jl")
include("splash_interpolation/splash_shared.jl")
include("splash_interpolation/splash_2D.jl")
include("splash_interpolation/splash_3D.jl")
include("splash_interpolation/splash_interpolation.jl")

# healpix interpolation
include("healpix_interpolation/shared.jl")
include("healpix_interpolation/main.jl")

# effect functions
include("effects/constants.jl")
include("effects/density.jl")
include("effects/synchrotron_analytic.jl")
include("effects/synchrotron_spectrum.jl")
include("effects/radio_beam.jl")
include("effects/sz_effect.jl")
include("effects/x_ray.jl")

export mappingParameters,                         # parameters for SPH mapping
    sphMapping,                                # main function for mapping 
    map_it,
    allsky_map,
    reduce_image_healpix,
    filter_particles_in_image,                 # helper function to preselect particles
    get_map_grid_2D,
    get_map_grid_3D,
    read_smac1_binary_info,
    read_smac1_binary_image,
    read_smac1_fits_info,
    read_smac1_fits_image,
    read_smac2_image,
    read_smac2_info,
    write_smac1_par,
    write_smac2_par,
    # rotate particles
    rotate_3D, rotate_3D!,
    project_along_axis,
    project_along_axis!,
    # Weight functions
    part_weight_one,
    part_weight_physical,
    part_weight_emission,
    part_weight_XrayBand,
    part_weight_spectroscopic,
    # effect functions
    density_2D,
    comptonY,
    kinetic_SZ,
    thermal_SZ,
    x_ray_emission,
    get_T_keV,
    analytic_synchrotron_emission,
    analytic_synchrotron_GS,
    spectral_synchrotron_emission,
    convert_Pnu_map_to_mJy_beam,
    # IO
    read_fits_image,
    write_fits_image


end # module
