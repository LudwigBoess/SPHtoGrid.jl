var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CurrentModule = SPHtoGrid\nDocTestSetup = quote\n    using SPHtoGrid\nend","category":"page"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Exported-Functions","page":"API reference","title":"Exported Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHtoGrid]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"api/#SPHtoGrid.filter_particles_in_image-Tuple{Array{#s31,N} where N where #s31<:Real,Array{#s22,N} where N where #s22<:Real,mappingParameters}","page":"API reference","title":"SPHtoGrid.filter_particles_in_image","text":"filter_particles_in_image(x, hsml, param::mappingParameters)\n\nChecks if a particle is contained in the image and returns an array of Bool.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_map_grid_2D-Tuple{mappingParameters}","page":"API reference","title":"SPHtoGrid.get_map_grid_2D","text":"get_map_grid_2D(par::mappingParameters)\n\nReconstruct the 2D grid used for mapping.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_map_grid_3D-Tuple{mappingParameters}","page":"API reference","title":"SPHtoGrid.get_map_grid_3D","text":"get_map_grid_3D(par::mappingParameters)\n\nReconstruct the 3D grid used for mapping.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.project_along_axis","page":"API reference","title":"SPHtoGrid.project_along_axis","text":"project_along_axis(x::AbstractArray, projection_axis::Integer=3)\n\nProjects and array of 3D along one of the principle axes. projection_axis ∈ {1, 2, 3} => x, y, z axis.\n\n\n\n\n\n","category":"function"},{"location":"api/#SPHtoGrid.read_smac1_binary_image-Tuple{String}","page":"API reference","title":"SPHtoGrid.read_smac1_binary_image","text":"read_smac1_binary_image(filename::String)\n\nReads a binary image file from Smac and returns a Matrix with the pixel values.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.read_smac1_binary_info-Tuple{String}","page":"API reference","title":"SPHtoGrid.read_smac1_binary_info","text":"read_smac1_binary_info(filename::String)\n\nReturns the image info in a Smac1ImageInfo struct.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.rotate_3D!-Tuple{Array{#s97,N} where N where #s97<:Real,Real,Real,Real}","page":"API reference","title":"SPHtoGrid.rotate_3D!","text":"rotate_3D!(x::AbstractArray, alpha::Real, beta::Real, gamma::Real)\n\nRotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively. α, β and γ need to be given in degrees.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.rotate_3D-Tuple{Array{#s97,N} where N where #s97<:Real,Real,Real,Real}","page":"API reference","title":"SPHtoGrid.rotate_3D","text":"rotate_3D(x::AbstractArray, alpha::Real, beta::Real, gamma::Real)\n\nRotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively. α, β and γ need to be given in degrees.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.sphMapping","page":"API reference","title":"SPHtoGrid.sphMapping","text":"function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, \n                    M::Array{<:Real}, ρ::Array{<:Real}, \n                    Bin_Quant::Array{<:Real},\n                    Weights::Array{<:Real}=ρ;\n                    param::mappingParameters,\n                    kernel::SPHKernel [,\n                    show_progress::Bool=true,\n                    parallel::Bool=false,\n                    filter_particles::Bool=true,\n                    dimensions::Int=2])\n\nMaps the data in Bin_Quant to a grid. Parameters of mapping are supplied in param and the kernel to be used in kernel.\n\nArguments\n\nPos: Array with particle positions.\nHSML: Array with particle hsml.\nM: Array with particle masses.\nρ: Array with particle densities.\nBin_Quant: Array with particle quantity to be mapped.\nWeights: Array with weights. Defaults to density-weighted.\nkernel::SPHKernel: Kernel object to be used.\nshow_progress::Bool=true: Show progress bar.\nparallel::Bool=true: Run on multiple processors.\nfilter_particles::Bool=true: Find the particles that are actually contained in the image.\ndimensions::Int=2: Number of mapping dimensions (2 = to grid, 3 = to cube).\n\n\n\n\n\n","category":"function"},{"location":"api/#SPHtoGrid.write_smac1_par","page":"API reference","title":"SPHtoGrid.write_smac1_par","text":"write_smac1_par([...])\n\nWrites a Smac parameter file. Not all relevant parameters implemented yet!\n\n\n\n\n\n","category":"function"},{"location":"api/#SPHtoGrid.write_smac2_par","page":"API reference","title":"SPHtoGrid.write_smac2_par","text":"write_smac2_par([...])\n\nWrites a P-Smac2 parameter file. Not all relevant parameters implemented yet!\n\n\n\n\n\n","category":"function"},{"location":"api/#Exported-Types","page":"API reference","title":"Exported Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHtoGrid]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"api/#SPHtoGrid.mappingParameters","page":"API reference","title":"SPHtoGrid.mappingParameters","text":"struct mappingParameters\n    x_lim::Vector{Float64}\n    y_lim::Vector{Float64}\n    z_lim::Vector{Float64}\n    center::Vector{Float64}\n    pixelSideLength::Float64\n    pixelArea::Float64\n    Npixels::Vector{Int64}\n    x_size::Float64\n    y_size::Float64\n    z_size::Float64\nend\n\nConstructor:\n\nmappingParameters(;x_lim::Vector{Float64}   = [-1.0, -1.0],\n                   y_lim::Vector{Float64}   = [-1.0, -1.0],\n                   z_lim::Vector{Float64}   = [-1.0, -1.0],\n                   center::Vector{Float64}  = [-1.0, -1.0, -1.0],\n                   x_size::Float64          =  -1.0,\n                   y_size::Float64          =  -1.0,\n                   z_size::Float64          =  -1.0,\n                   pixelSideLength::Float64 =  -1.0,\n                   Npixels::Int64           =   0)\n\nParameter object for sph to grid mapping. Define either *_lim, or center and *_size.  Resolution is defined by pixelSideLength or Npixels.\n\n\n\n\n\n","category":"type"},{"location":"api/#Private-Functions","page":"API reference","title":"Private Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHtoGrid]\nPublic = false\nOrder = [:function]","category":"page"},{"location":"api/#SPHtoGrid.calculate_index_2D-Tuple{Integer,Integer,Integer}","page":"API reference","title":"SPHtoGrid.calculate_index_2D","text":"function calculate_index_2D(i::Integer, j::Integer, x_pixels::Integer)\n\nCalculates the index of a flattened 2D image array.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.calculate_index_3D-NTuple{5,Integer}","page":"API reference","title":"SPHtoGrid.calculate_index_3D","text":"function calculate_index_3D(i::Integer, j::Integer, x_pixels::Integer)\n\nCalculates the index of a flattened 3D image array.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.calculate_weights_2D-Tuple{Array{#s21,1} where #s21<:Real,Integer,Integer,Integer,Integer,Real,Real,Real,Real,SPHKernels.SPHKernel,Integer}","page":"API reference","title":"SPHtoGrid.calculate_weights_2D","text":"function calculate_weights_2D(  wk::Array{<:Real,1}, \n                                iMin::Integer, iMax::Integer, \n                                jMin::Integer, jMax::Integer,\n                                x::Real, y::Real, hsml::Real, hsml_inv::Real,\n                                kernel::SPHKernel,\n                                x_pixels::Integer )\n\nCalculates the kernel- and geometric weights of the pixels a particle contributes to.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.calculate_weights_3D-Tuple{Array{#s31,1} where #s31<:Real,Integer,Integer,Integer,Integer,Integer,Integer,Real,Real,Real,Real,Real,SPHKernels.SPHKernel,Integer,Integer}","page":"API reference","title":"SPHtoGrid.calculate_weights_3D","text":"function calculate_weights_3D(  wk::Array{<:Real,1}, \n                                iMin::Integer, iMax::Integer, \n                                jMin::Integer, jMax::Integer,\n                                kMin::Integer, kMax::Integer,\n                                x::Real, y::Real, z::Real, \n                                hsml::Real, hsml_inv::Real,\n                                kernel::SPHKernel,\n                                x_pixels::Integer, y_pixels::Integer )\n\nCalculates the kernel- and geometric weights of the pixels a particle contributes to.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.check_center_and_move_particles-Tuple{Array{#s32,N} where N where #s32<:Real,mappingParameters}","page":"API reference","title":"SPHtoGrid.check_center_and_move_particles","text":"check_center_and_move_particles(x, par::mappingParameters)\n\nMapping only works if all coordinates are positive. This function shifts the particles accordingly.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.check_in_image-NTuple{7,Real}","page":"API reference","title":"SPHtoGrid.check_in_image","text":"check_in_image(x::Real, y::Real, z::Real, hsml::Real,\n                            halfXsize::Real, halfYsize::Real, halfZsize::Real)\n\nChecks if a particle is in the image frame.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.domain_decomposition-Tuple{Integer,Integer}","page":"API reference","title":"SPHtoGrid.domain_decomposition","text":"domain_decomposition(N::Int64, N_workers::Int64)\n\nCalculate relevant array slices for each worker. Could be done better!\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.euler_matrix-Tuple{Real,Real,Real}","page":"API reference","title":"SPHtoGrid.euler_matrix","text":"euler_matrix(α::Real, β::Real, γ::Real)\n\nReturns the rotation matrix A based on rotation along the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.find_position_periodic-Tuple{Array{#s34,N} where N where #s34<:Real,Integer,Real}","page":"API reference","title":"SPHtoGrid.find_position_periodic","text":"find_position_periodic( pos::Array{<:Real}, k::Integer, boxsize::Real)\n\nPerforms a periodic mapping of the particle position.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_d_hsml_2D-Tuple{Real,Real,Real}","page":"API reference","title":"SPHtoGrid.get_d_hsml_2D","text":"get_d_hsml_2D( dx::Real, dy::Real,\n               hsml_inv::Real )\n\nComputes the distance in 2D to the pixel center in units of the kernel support.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_d_hsml_3D-NTuple{4,Real}","page":"API reference","title":"SPHtoGrid.get_d_hsml_3D","text":"get_d_hsml_3D( dx::Real, dy::Real, dz::Real,\n               hsml_inv::Real )\n\nComputes the distance in 3D to the pixel center in units of the kernel support.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_dxyz-Tuple{Real,Real,Integer}","page":"API reference","title":"SPHtoGrid.get_dxyz","text":"function get_dxyz(x::Real, hsml::Real, i::Integer)\n\nCalculates the extent of the current particle size in units of pixels.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_ijk_min_max-Tuple{Real,Real,Integer}","page":"API reference","title":"SPHtoGrid.get_ijk_min_max","text":"function get_ijk_min_max( x::Real, hsml::Real,\n                          x_pixels::Integer)\n\nCalculates the minimum and maximum pixel to which a particle contributes.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.get_xyz-Tuple{Array{#s33,1} where #s33<:Real,Real,Integer,Real,Integer,Integer,Integer,Real,Bool,Real,Real,Real,Integer}","page":"API reference","title":"SPHtoGrid.get_xyz","text":"function get_xyz( pos::Vector{<:Real}, hsml::Real, k::Integer,\n                  len2pix::Real, x_pixels::Integer, y_pixels::Integer, z_pixels::Integer,\n                  boxsize::Real, periodic::Bool,\n                  halfXsize::Real, halfYsize::Real, halfZsize::Real,\n                  Ndim::Integer)\n\nCalculates x, y, z position in units of pixels and performs periodic mapping, if required.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.reduce_futures-Tuple{Array{#s33,N} where N where #s33<:Tuple}","page":"API reference","title":"SPHtoGrid.reduce_futures","text":"function reduce_futures(fut::Array{<:Tuple})\n\nReduces the touple returned by the Array of Futures to image arrays. \n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.reduce_image_2D-Tuple{Array{#s31,N} where N where #s31<:Real,Array{#s22,N} where N where #s22<:Real,Int64,Int64}","page":"API reference","title":"SPHtoGrid.reduce_image_2D","text":"function reduce_image_2D( image::Array{<:Real}, w_image::Array{<:Real},\n                                        x_pixels::Int64, y_pixels::Int64)\n\nUnflattens an image array to a 2D array.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.reduce_image_3D-Tuple{Array{#s31,N} where N where #s31<:Real,Int64,Int64,Int64}","page":"API reference","title":"SPHtoGrid.reduce_image_3D","text":"function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},\n                                        x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)\n\nUnflattens an image array to a 3D array.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.rotate_3D_quantity-Tuple{Array{#s103,N} where N where #s103<:Real,Real,Real,Real}","page":"API reference","title":"SPHtoGrid.rotate_3D_quantity","text":"rotate_3D_quantity(x, α, β, γ)\n\nRotates a 3D vector along the x-axis with angle α, y-axis with angle β and z-axis with angle γ.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.rotate_to_xz_plane-Tuple{Array{#s103,N} where N where #s103<:Real}","page":"API reference","title":"SPHtoGrid.rotate_to_xz_plane","text":"rotate_to_xz_plane(x::AbstractArray)\n\nRotates an array of 3D positions into the xz-plane.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.rotate_to_yz_plane-Tuple{Array{#s103,N} where N where #s103<:Real}","page":"API reference","title":"SPHtoGrid.rotate_to_yz_plane","text":"rotate_to_yz_plane(x::AbstractArray)\n\nRotates an array of 3D positions into the yz-plane.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.sphMapping_2D-Tuple{Array{#s89,N} where N where #s89<:Real,Array{#s90,N} where N where #s90<:Real,Array{#s91,N} where N where #s91<:Real,Array{#s92,N} where N where #s92<:Real,Array{#s93,N} where N where #s93<:Real,Array{#s94,N} where N where #s94<:Real}","page":"API reference","title":"SPHtoGrid.sphMapping_2D","text":"sphMapping2D( Pos::Array{<:Real}, HSML::Array{<:Real},                    M::Array{<:Real}, Rho::Array{<:Real},                    BinQ::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));                   param::mappingParameters, kernel::SPHKernel,                   show_progress::Bool=false )\n\nUnderlying function to map SPH data to a 2D grid.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.sphMapping_3D-Tuple{Array{#s18,N} where N where #s18<:Real,Array{#s17,N} where N where #s17<:Real,Array{#s16,N} where N where #s16<:Real,Array{#s15,N} where N where #s15<:Real,Array{#s14,N} where N where #s14<:Real,Array{#s13,N} where N where #s13<:Real}","page":"API reference","title":"SPHtoGrid.sphMapping_3D","text":"sphMapping2D( Pos::Array{<:Real}, HSML::Array{<:Real},                    M::Array{<:Real}, Rho::Array{<:Real},                    BinQ::Array{<:Real}, Weights::Array{<:Real}=ones(length(Rho));                   param::mappingParameters, kernel::SPHKernel,                   show_progress::Bool=false )\n\nUnderlying function to map SPH data to a 3D grid.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHtoGrid.update_image-NTuple{5,Real}","page":"API reference","title":"SPHtoGrid.update_image","text":"function update_image( image::Real, w_image::Real, \n                       wk::Real, bin_q::Real, \n                       geometry_norm::Real )\n\nApplies the different contributions to the image and the weight image.\n\n\n\n\n\n","category":"method"},{"location":"api/#Private-Types","page":"API reference","title":"Private Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHtoGrid]\nPublic = false\nOrder = [:type]","category":"page"},{"location":"api/#SPHtoGrid.Smac1ImageInfo","page":"API reference","title":"SPHtoGrid.Smac1ImageInfo","text":"struct Smac1ImageInfo\n\nStores the information in a Smac binary image header.\n\n\n\n\n\n","category":"type"},{"location":"external/#External-Programs","page":"External Programs","title":"External Programs","text":"","category":"section"},{"location":"external/","page":"External Programs","title":"External Programs","text":"SPHtoGrid.jl provides helper function for two external sph mapping Codes: Smac and P-Smac2.","category":"page"},{"location":"external/#P-Smac2","page":"External Programs","title":"P-Smac2","text":"","category":"section"},{"location":"external/","page":"External Programs","title":"External Programs","text":"P-Smac2 by Julius Donnert is an advanced mapping code for a multitude of different quantities. To run a mapping and plotting loop from a Julia script you need to update the parameter files on the fly. The function write_smac2_par provides this functionality.","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"write_smac2_par(x, y, z,\n                euler_angle_0, euler_angle_1, euler_angle_2,\n                xy_size, z_depth, xy_pix::Int64,\n                input_file, output_file, path,\n                effect_module::Int64=0, effect_flag::Int64=0)","category":"page"},{"location":"external/#Smac","page":"External Programs","title":"Smac","text":"","category":"section"},{"location":"external/","page":"External Programs","title":"External Programs","text":"Smac is a SPH mapping Code by Klaus Dolag and others. The implementation is described in Dolag et al. 2005.","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"Smac isn't public unfortunately. So these functions are mainly for my personal use. If you do have access to Smac, here's a reference to what you can do.","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"SPHtoGrid.jl provides some functions to read the binary output of Smac, as I personally prefer that over the FITS output. To get the binary format you need to set FILE_FORMAT = 1 in the parameter file.","category":"page"},{"location":"external/#Reading-image-information","page":"External Programs","title":"Reading image information","text":"","category":"section"},{"location":"external/","page":"External Programs","title":"External Programs","text":"If you set FILE_HEADER = 1 in the Smac parameter file you can read the information of the image header into a Smac1ImageInfo object like so:","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"info = read_smac1_binary_info(filename)","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"The Smac1ImageInfo object contains the following information","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"struct Smac1ImageInfo\n\n    snap::Int32                 # number of input snapshot\n    z::Float32                  # redshift of snapshot\n    m_vir::Float32              # virial mass of halo\n    r_vir::Float32              # virial radius of halo\n    xcm::Float32                # x coordinate of image center\n    ycm::Float32                # y coordinate of image center\n    zcm::Float32                # z coordinate of image center\n    z_slice_kpc::Float32        # depth of the image in kpc\n    boxsize_kpc::Float32        # xy-size of the image in kpc\n    boxsize_pix::Float32        # xy-size of the image in pixels\n    pixsize_kpc::Float32        # size of one pixel in kpc\n    xlim::Array{Float64,1}      # x limits of image\n    ylim::Array{Float64,1}      # y limits of image\n    zlim::Array{Float64,1}      # z limits of image\n    units::String               # unitstring of image\n\nend\n","category":"page"},{"location":"external/#Reading-the-image","page":"External Programs","title":"Reading the image","text":"","category":"section"},{"location":"external/","page":"External Programs","title":"External Programs","text":"The image itself can be read with","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"image = read_smac1_binary_image(filename)","category":"page"},{"location":"external/","page":"External Programs","title":"External Programs","text":"This will return an Array{Float32,2} with the pixel values. You can pass this to any imshow function of your favorite plotting package.","category":"page"},{"location":"mapping/#Mapping-SPH-Data","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"CurrentModule = SPHtoGrid\nDocTestSetup = quote\n    using SPHtoGrid\nend","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"You can map SPH data to a grid using the function sphMapping:","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, \n                    M::Array{<:Real}, ρ::Array{<:Real}, \n                    Bin_Quant::Array{<:Real},\n                    Weights::Array{<:Real}=ρ;\n                    param::mappingParameters,\n                    kernel::SPHKernel [,\n                    show_progress::Bool=true,\n                    parallel::Bool=false,\n                    filter_particles::Bool=true,\n                    dimensions::Int=2])\n\n\n    [...]\n\nend","category":"page"},{"location":"mapping/#Define-parameters-for-mapping","page":"Mapping SPH Data","title":"Define parameters for mapping","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"To map the data you need to define the mapping parameters via the mappingParameters object. One way to set this up is by defining the limits of the map as","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"par = mappingParameters(xlim=[xmin, xmax],\n                        ylim=[ymin, ymax],\n                        zlim=[zmin, zmax],\n                        Npixels=200)","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"or give a center position and the size in each direction","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"par = mappingParameters(center=[x0, y0, z0], \n                        x_size=x_size, \n                        y_size=y_size,\n                        z_size=z_size,\n                        Npixels=200)","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"Instead of Npixels you can also give the keyword argument pixelSideLength if you prefer to define your image that way.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"If you are mapping a periodic box you also can give the keyword boxsize to enable periodic mapping.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"par = mappingParameters(center=[x0, y0, z0], \n                        x_size=x_size, \n                        y_size=y_size,\n                        z_size=z_size,\n                        boxsize=boxsize,\n                        Npixels=200)","category":"page"},{"location":"mapping/#Select-Kernel","page":"Mapping SPH Data","title":"Select Kernel","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"You also need to choose the kernel you used in the simulation. For this you need to install the package SPHKernels.jl. You can currently use these kernels:","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"k = Cubic()\nk = Quintic()\nk = WendlandC4()\nk = WendlandC6()","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"Please see the SPHKernels docs for more details.","category":"page"},{"location":"mapping/#Mapping","page":"Mapping SPH Data","title":"Mapping","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"With the setup done you can now map (e.g.) density of your data using the function above as:","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"image = sphMapping(x, hsml, m, rho, rho, param=par, kernel=k)","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"Replacing the second rho with any other quantity would map that quantity of course. Please note: This function doesn't do any unit conversion for you, so you need to convert to the desired units beforehand. You can do this e.g. with GadgetUnits.jl.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"Image now contains a 2D array with the binned data and can easily be plotted with imshow() from any plotting package of your choosing.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"The keyword parallel = true causes the run to use multiple processors. For this you need to start julia with julia -p <N> where <N> is the number of processors in your machine, or define","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"using Distributed\naddprocs(8)\n\n# now you can load SPHtoGrid\nusing SPHtoGrid","category":"page"},{"location":"mapping/#Conserved-quantities","page":"Mapping SPH Data","title":"Conserved quantities","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"With the latest release you can map the particles to a grid while also conserving the particle volume, following the algorithm described in Dolag et. al. 2006.","category":"page"},{"location":"mapping/#Weight-functions","page":"Mapping SPH Data","title":"Weight functions","text":"","category":"section"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"With the mapping you may decide to use a specivic weighting function. For this you can pass the optional variable Weights in sphMapping.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"You can either use your own weight functions or use one of the built-in ones:","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"part_weight_one just returns an Array of ones.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"part_weight_physical converts from pixel- to physical units.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"part_weight_emission weights the contribution due to density and temperature of the particle.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"part_weight_spectroscopic gives spectroscopic weighting, see Mazotta+ 04.","category":"page"},{"location":"mapping/","page":"Mapping SPH Data","title":"Mapping SPH Data","text":"part_weight_XrayBand weights the particle due to its Xray emission in the defined energy band.","category":"page"},{"location":"install/#Install","page":"Install","title":"Install","text":"","category":"section"},{"location":"install/","page":"Install","title":"Install","text":"As usual with Julia you can install the package via the internal package manager","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"julia> ]\npkg> add https://github.com/LudwigBoess/SPHtoGrid.jl","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"If you want to get the latest version use","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"julia> ]\npkg> add https://github.com/LudwigBoess/SPHtoGrid.jl#development","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"Now you should be good to go!","category":"page"},{"location":"#Table-of-Contents","page":"Table of Contents","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"Pages = [ \"index.md\", \n          \"install.md\", \n          \"mapping.md\",\n          \"weights.md\",\n          \"rotating.md\",\n          \"external.md\", \n          \"api.md\"]\nDepth = 3","category":"page"},{"location":"rotating/#Rotating-images","page":"Rotating Images","title":"Rotating images","text":"","category":"section"},{"location":"rotating/","page":"Rotating Images","title":"Rotating Images","text":"By default sphMapping only maps the xy-plane. To change the mapping you have two options, Project along axis and Define Euler Angle.","category":"page"},{"location":"rotating/#Project-along-axis","page":"Rotating Images","title":"Project along axis","text":"","category":"section"},{"location":"rotating/","page":"Rotating Images","title":"Rotating Images","text":"If you only want to change the axis along which you want to project the data you can use the wrapper function project_along_axis. If you for example want to project along the x-axes, so in the yz-plane use","category":"page"},{"location":"rotating/","page":"Rotating Images","title":"Rotating Images","text":"axis = 1\npos_new = (pos_old, axis)","category":"page"},{"location":"rotating/#Define-Euler-Angle","page":"Rotating Images","title":"Define Euler Angle","text":"","category":"section"},{"location":"rotating/","page":"Rotating Images","title":"Rotating Images","text":"If projection along one of the principle axis is too crude for you, you can define individual angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively and use the function rotate_3D. These angles have to be given in degrees. So to rotate a 3D quantity 45 degrees around the x-axis you can use:","category":"page"},{"location":"rotating/","page":"Rotating Images","title":"Rotating Images","text":"pos_new = rotate_3D(pos_old, 45.0, 0.0, 0.0)","category":"page"}]
}
