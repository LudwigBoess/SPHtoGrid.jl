module SPHtoGrid

    using Distributed
    using SPHKernels
    using Unitful
    using Printf

    function output_time(t1, t2)
        return @sprintf("%0.3e", Float64((t2-t1))*1.e-9)
    end
    
    include("parameters.jl")
    include("filter_shift.jl")
    include("domain_decomp.jl")
    include("reduce_image.jl")
    include("reconstruct_grid.jl")
    include("mapping_functions.jl")
    include("smac1_utility.jl")
    include("smac2_utility.jl")
    include("rotate_particles.jl")
    include("cic_interpolation.jl")
    include("tsc_interpolation.jl")
    include("weight_functions.jl")
    include("effects.jl")
    include("io.jl")
    


    export mappingParameters,                         # parameters for SPH mapping
           sphMapping,                                # main function for mapping 
           filter_particles_in_image,                 # helper function to preselect particles
           get_map_grid_2D,
           get_map_grid_3D,
           read_smac1_binary_info,
           read_smac1_binary_image,
           write_smac1_par,
           write_smac2_par,
           # rotate particles
           rotate_3D, rotate_3D!,
           project_along_axis,
           # Weight functions
           part_weight_one,
           part_weight_physical, 
           part_weight_emission,
           part_weight_XrayBand,
           # effect functions
           density_2D, 
           kinetic_SZ,
           thermal_SZ,
           # IO
           read_fits_image,
           write_fits_image

           


end # module
