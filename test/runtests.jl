using Distributed
addprocs(2)

@everywhere using SPHtoGrid, Test, DelimitedFiles, GadgetIO


@testset "SPH Kernels" begin

    k = Cubic()
    # < 0.5
    d = kernel_value_2D(k, 0.4, 0.5)
    @test d ≈ 0.26992678348385446
    d = kernel_value_3D(k, 0.4, 0.5)
    @test d ≈ 0.13496339174192723
    # < 1.0
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.15915494309189535
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.07957747154594767
    # > 1.0
    d = kernel_value_2D(k, 1.5, 0.5)
    @test d == 0.0
    d = kernel_value_3D(k, 1.5, 0.5)
    @test d == 0.0


    k = Quintic()
    # < 1/3
    d = kernel_value_2D(k, 0.3, 0.5)
    @test d ≈ 0.32700352517410614
    d = kernel_value_3D(k, 0.3, 0.5)
    @test d ≈ 0.2791208661307549
    # 2/3
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.07767855829318415
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.06630419797168217
    # < 1.0
    d = kernel_value_2D(k, 0.8, 0.5)
    @test d ≈ 0.0008155658657050455
    d = kernel_value_3D(k, 0.8, 0.5)
    @test d ≈ 0.0006961437210839494
    # > 1.0
    d = kernel_value_2D(k, 1.5, 0.5)
    @test d == 0.0
    d = kernel_value_3D(k, 1.5, 0.5)
    @test d == 0.0



    k = WendlandC4()
    # < 1.0
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.30960610023345264
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.26606774238812336
    # > 1.0
    d = kernel_value_2D(k, 1.5, 0.5)
    @test d == 0.0
    d = kernel_value_3D(k, 1.5, 0.5)
    @test d == 0.0




    k = WendlandC6()
    # < 1.0
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.052822211162893276
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.05055250677698771
    # > 1.0
    d = kernel_value_2D(k, 1.5, 0.5)
    @test d == 0.0
    d = kernel_value_3D(k, 1.5, 0.5)
    @test d == 0.0

end

@testset "Smac utility" begin

    @test_nowarn write_smac1_par("", 0, "", "", "", "",
                                    0, 0, 0, 4, 3,
                                    20.0, 10.0, 1, 1,
                                    24, 1.0, 1.e6, 10,
                                    1, 0.0, 0.0, 0.0)

    # @test_throws ErrorException("Read error: Incorrect image format!") read_smac1_binary_image(joinpath(dirname(@__FILE__), "snap_050"))

    filename = joinpath(dirname(@__FILE__), "Smac1.pix")
    info = read_smac1_binary_info(filename)

    @test info.snap == 140

    image = read_smac1_binary_image(filename)
    @test length(image[:,1]) == 128
    @test image[1,1] ≈ 0.000120693

    @test_nowarn write_smac1_par("", 0, "", "", "", "",
                                    0, 0, 0, 0, 0, 0.0, 0.0, 
                                    1, 1, 0, 1.0, 1.e6, 10, 0, 0.0, 0.0, 0.0)

    @test_nowarn write_smac2_par(1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1024,
                                    "", "", "")
end

@testset "SPH mappingParameters" begin

    @test_throws ErrorException("Giving a center position requires extent in x, y and z direction.") mappingParameters()
    
    @test_throws ErrorException("Please specify pixelSideLenght or number of pixels!") mappingParameters(center=[0.0, 0.0, 0.0],
                                                                                       x_lim = [-1.0, 1.0],
                                                                                       y_lim = [-1.0, 1.0],
                                                                                       z_lim = [-1.0, 1.0])

    @test_nowarn mappingParameters(center=[0.0, 0.0, 0.0],
                                    x_lim = [-1.0, 1.0],
                                    y_lim = [-1.0, 1.0],
                                    z_lim = [-1.0, 1.0],
                                    Npixels=100)

    @test_nowarn mappingParameters(center=[0.0, 0.0, 0.0],
                                    x_lim = [-1.0, 1.0],
                                    y_lim = [-1.0, 1.0],
                                    z_lim = [-1.0, 1.0],
                                    pixelSideLength=0.2)
end

@testset "Filter particles" begin
    
    par = mappingParameters(center = [3.0, 3.0, 3.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 500)

    x = zeros(2, 3)
    x[1,:] = [  1.0,  1.0, 1.0]
    x[2,:] = [ -3.0, -1.0, 1.0]

    hsml = [0.5, 0.5]

    p_in_image = filter_particles_in_image(x, hsml, par)

    @test p_in_image[1] == true
    @test p_in_image[2] == false

end

@testset "Shift particles" begin
    
    par = mappingParameters(center = [0.0, 0.0, 0.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 500)

    x = zeros(2, 3)
    x[1,:] = [  1.0,  1.0, 1.0]
    x[2,:] = [ -3.0, -1.0, 1.0]

    hsml = [0.5, 0.5]

    x, par = SPHtoGrid.check_center_and_move_particles(x, par)

    @test x[1, :] ≈ [ 4.0, 4.0, 4.0]
    @test x[2, :] ≈ [ 0.0, 2.0, 4.0]

    @test par.center ≈ [ 3.0, 3.0, 3.0 ]

end

@testset "SPH Mapping" begin

    @info "SPH Mapping tests take a while..."

    @info "Data read-in."
    
    snap_file = joinpath(dirname(@__FILE__), "snap_050")

    h = head_to_obj(snap_file)
    info = read_info(snap_file, verbose=false)

    bin_quantity = read_block_by_name(snap_file, "RHO",
                                        info=info[getfield.(info, :block_name) .== "RHO"][1],
                                        parttype=0)

    x = read_block_by_name(snap_file, "POS",
                            info=info[getfield.(info, :block_name) .== "POS"][1],
                            parttype=0)

    rho = read_block_by_name(snap_file, "RHO",
                        info=info[getfield.(info, :block_name) .== "RHO"][1],
                        parttype=0)

    hsml = read_block_by_name(snap_file, "HSML",
                            info=info[getfield.(info, :block_name) .== "HSML"][1],
                        parttype=0)

    m = read_block_by_name(snap_file, "MASS", parttype=0)

    kernel = WendlandC6()

    par = mappingParameters(center = [3.0, 3.0, 3.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 500)

    @info "2D"
    
    @info "Single core, no unit conservation."
    d = sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=false,
						  parallel = false,
						  show_progress=false)


    ideal_file = joinpath(dirname(@__FILE__), "image.dat")
    d_ideal = readdlm(ideal_file)

    @test d[  1,  1] ≈ d_ideal[1, 1]
    @test d[ 30, 32] ≈ d_ideal[30, 32]
    @test d[117, 92] ≈ d_ideal[117, 92]


    # @test_nowarn sphMapping(x, hsml, m, rho, bin_quantity,
	# 					  param=par, kernel=kernel,
	# 					  conserve_quantities=false,
	# 					  parallel = false,
    #                       show_progress=true)

    par = mappingParameters(center = [3.0, 3.0, 3.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 50)

    @info "Single core, unit conservation."
    @test_nowarn sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=true,
						  parallel = false,
                          show_progress=false)

                          
    
    @info "Multicore, unit conservation."
    @test_nowarn sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=false,
						  parallel = true,
                          show_progress=false)
                          
    @info "3D"

    par = mappingParameters(center = [3.0, 3.0, 3.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 20)

    @info "Single core, no unit conservation."
    @test_nowarn  sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=false,
						  parallel = false,
                          show_progress=false,
                          dimensions=3)

    par = mappingParameters(center = [3.0, 3.0, 3.0],
							x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 20)

    @info "Single core, unit conservation."
    @test_nowarn  sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=true,
						  parallel = false,
                          show_progress=false,
                          dimensions=3)

    @info "Multi core, no unit conservation."
    @test_nowarn  sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=false,
						  parallel = true,
                          show_progress=false,
                          dimensions=3)

end


