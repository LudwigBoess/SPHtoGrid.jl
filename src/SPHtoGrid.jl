module SPHtoGrid

using Distributed
using SPHKernels
using Unitful
using Printf
using ProgressMeter
using DiffusiveShockAccelerationModels
using SpecialFunctions
using LinearAlgebra
using Rotations
using Base.Threads
using Healpix
using NearestNeighbors
using Distances
using Statistics
using FFTW

function output_time(t1, t2)
    return @sprintf("%0.3e", Float64((t2 - t1)) * 1.e-9)
end

const global tables_path = joinpath(@__DIR__, "tables")

# read cooling tables for Xray emission 
include("tables/read_ccoling_tables.jl")

# shared functionality
include("shared/parameters.jl")
include("shared/reconstruct_grid.jl")
include("shared/indices.jl")
include("shared/distances.jl")
include("shared/periodic_mapping.jl")
include("shared/rotate_particles.jl")
include("shared/rotate_parameters.jl")
include("shared/io.jl")
include("shared/vtk.jl")

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
include("healpix_interpolation/filter_particles.jl")
include("healpix_interpolation/constributing_pixels.jl")
include("healpix_interpolation/pixel_weights.jl")
include("healpix_interpolation/main.jl")

# distributed mapping
include("distributed_mapping/main.jl")
include("distributed_mapping/healpix.jl")
include("distributed_mapping/cic.jl")

# effect functions
include("effects/constants.jl")
include("effects/utils.jl")
include("effects/density.jl")
include("effects/sz_effect.jl")
include("effects/x_ray.jl")
include("effects/gamma.jl")
include("effects/RM.jl")
include("effects/synchrotron_shared.jl")
include("effects/synchrotron_GS.jl")
include("effects/synchrotron_Hoeft.jl")
include("effects/synchrotron_Longair.jl")
include("effects/synchrotron_LMB.jl")
include("effects/mass_density.jl")

# functions for existing images 
include("image_functions/radio_beam.jl")
include("image_functions/stokes_parameters.jl")
include("image_functions/synchrotron_luminosity.jl")
include("image_functions/surface_brightness.jl")

# functions for powerspectrum computation
include("powerspectrum/powerspectrum.jl")

# precompile step 
include("precompile.jl")

export mappingParameters,                         # parameters for SPH mapping
    sphMapping,                                # main function for mapping 
    map_it,
    healpix_map,
    distributed_cic_map,
    distributed_allsky_map,
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
    x_ray_emissivity,
    get_T_keV,
    rotation_measure,
    analytic_synchrotron,
    analytic_synchrotron_HB07,
    analytic_synchrotron_GS,
    analytic_synchrotron_Longair,
    λγ_PE04, jγ_PE04,
    gamma_luminosity_pions_PE04,
    gamma_flux_pions_PE04,
    mass_density,
    get_smac_T_rho_filter,
    apply_smac_T_rho_filter!,
    # image functions
    beam_in_kpc,
    convert_Pnu_map_to_mJy_beam,
    polarisation_fraction,
    polarisation_angle,
    surface_brightness_to_luminosity,
    synchrotron_SB_to_luminosity,
    total_synch_luminosity_from_SB,
    # IO
    read_fits_image,
    read_allsky_fits_image,
    write_fits_image,
    write_vtk_image,
    # power spectrum
    power_spectrum


end # module
