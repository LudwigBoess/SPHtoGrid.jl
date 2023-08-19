using Downloads

@info "downloading test data..."
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_sedov", "./snap_sedov")
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_cutout_072", "./snap_cutout_072")

@info "done"

using Distributed
addprocs(2)

@everywhere using SPHtoGrid, Test, DelimitedFiles, SPHKernels, GadgetIO, GadgetUnits, DiffusiveShockAccelerationModels

@testset "SPHtoGrid" begin

    @testset "Smac utility" begin

        @test_nowarn write_smac1_par("./")

        filename = joinpath(dirname(@__FILE__), "Smac1.pix")
        info = read_smac1_binary_info(filename)

        @test info.snap == 140

        image = read_smac1_binary_image(filename)
        @test length(image[:, 1]) == 128
        @test image[1, 1] ≈ 0.000120693

        @test_nowarn write_smac1_par("./")


        @test_nowarn write_smac2_par(1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1024,
            "", "", "")
    end

    @testset "SPH mappingParameters" begin

        @test_throws ErrorException("Giving a center position requires extent in x, y and z direction.") mappingParameters()

        @test_throws ErrorException("Please specify pixelSideLength or number of pixels!") mappingParameters(center = [0.0, 0.0, 0.0],
            x_lim = [-1.0, 1.0],
            y_lim = [-1.0, 1.0],
            z_lim = [-1.0, 1.0])

        @test_nowarn mappingParameters(center = [0.0, 0.0, 0.0],
            x_lim = [-1.0, 1.0],
            y_lim = [-1.0, 1.0],
            z_lim = [-1.0, 1.0],
            Npixels = 100)

        @test_nowarn mappingParameters(center = [0.0, 0.0, 0.0],
            x_lim = [-1.0, 1.0],
            y_lim = [-1.0, 1.0],
            z_lim = [-1.0, 1.0],
            pixelSideLength = 0.2)
    end

    @testset "Filter particles" begin

        par = mappingParameters(center = [3.0, 3.0, 3.0],
            x_size = 6.0, y_size = 6.0, z_size = 6.0,
            Npixels = 500)

        x = zeros(3, 2)
        x[:, 1] = [1.0, 1.0, 1.0]
        x[:, 2] = [-5.0, -1.0, 1.0]

        hsml = [0.5, 0.5]

        p_in_image = filter_particles_in_image(x, par)

        @test p_in_image[1] == true
        @test p_in_image[2] == false

    end

    @testset "Shift particles" begin

        par = mappingParameters(center = [1.0, 1.0, 1.0],
            x_size = 6.0, y_size = 6.0, z_size = 6.0,
            Npixels = 500)

        x = zeros(3, 2)
        x[:, 1] = [1.0, 1.0, 1.0]
        x[:, 2] = [-3.0, -1.0, 1.0]

        hsml = [0.5, 0.5]

        x, par2 = SPHtoGrid.center_particles(x, par)

        @test x[:, 1] ≈ [0.0, 0.0, 0.0]
        @test x[:, 2] ≈ [-4.0, -2.0, 0.0]

        @test par.center == [1.0, 1.0, 1.0]
        @test par2.center == [0.0, 0.0, 0.0]

    end

    @testset "Rotate particles" begin

        # matrix no rotation
        x_rand = rand(3, 10)
        x_out = rotate_3D(x_rand, 0.0, 0.0, 0.0)
        @test x_out ≈ x_rand

        # inplace
        rotate_3D!(x_out, 0.0, 0.0, 0.0)
        @test x_out ≈ x_rand

        # project along axis
        x_in = [1.0 1.0
                1.0 1.0
                0.0 0.0]

        # along z-axis should not change anything
        x_out = project_along_axis(x_in, 3)

        @test x_out ≈ x_in

        # along y-axis
        x_out = project_along_axis(x_in, 2)

        # @test x_out ≈ copy(transpose([ 1.0 0.0 1.0
        #                                 1.0 0.0 1.0 ]))

        @test x_out ≈ [1.0 1.0
            0.0 0.0
            1.0 1.0]

        # along x-axis
        x_out = project_along_axis(x_in, 2)

        @test x_out ≈ copy(transpose([1.0 0.0 1.0
            1.0 0.0 1.0]))
    end

    @testset "Indices" begin 

        Npixels = 128

        @testset "2D" begin
            indices = Vector{Int64}(undef, Npixels^2)
            count = 1
            for i = 0:Npixels-1, j = 0:Npixels-1
                idx = SPHtoGrid.calculate_index(i, j, Npixels)
                indices[idx] = count
                count += 1
            end

            @test indices == collect(1:Npixels^2)
        end

        @testset "3D" begin
            indices = Vector{Int64}(undef, Npixels^3)
            count = 1
            for i = 0:Npixels-1, j = 0:Npixels-1, k = 0:Npixels-1
                idx = SPHtoGrid.calculate_index(i, j, k, Npixels, Npixels)
                indices[idx] = count
                count += 1
            end
            
            @test indices == collect(1:Npixels^3)
        end
    end

    @testset "SPH Mapping" begin

        @info "SPH Mapping tests take a while..."

        @info "Data read-in."

        # reference images
        T_ref, par, snap, units   = read_fits_image("sedov_T_reference.fits")
        rho_ref, par, snap, units = read_fits_image("sedov_rho_reference.fits")

        # sedov reference snapshot
        fi = "snap_sedov"

        # kernel definition
        k = WendlandC4(2)

        # read header and find units
        h = read_header(fi)
        GU = GadgetPhysical(xH=0.752)

        blocks = ["POS", "MASS", "HSML", "RHO", "U"]
        data = Dict(block => read_block(fi, block, parttype=0) for block ∈ blocks)

        # convert to physical code units for mapping
        pos  = data["POS"]  .* GU.x_physical
        hsml = data["HSML"] .* GU.x_physical
        rho  = data["RHO"]  .* GU.rho_physical
        mass = data["MASS"] .* GU.m_physical

        rho_gcm3 = data["RHO"]  .* GU.rho_cgs
        m_cgs    = data["MASS"] .* GU.m_cgs
        T = data["U"] .* GU.T_K


        Npart = length(m_cgs)

        center = ones(3) .* 0.5h.boxsize * GU.x_physical
        xy_size = 0.9h.boxsize * GU.x_physical
        z_size = 0.9h.boxsize * GU.x_physical

        # define mapping parameters
        par = mappingParameters(center=center .* GU.x_physical,
            x_size=xy_size,
            y_size=xy_size,
            z_size=z_size,
            Npixels=256,
            boxsize=h.boxsize * GU.x_physical)

        snap = 50

        @testset "2D" begin
            
            kernel = WendlandC6(2)

            @testset "Single core" begin
                @testset "rho" begin 
                    image_prefix = "sedov_rho"
                    weights = part_weight_physical(length(hsml), par, GU.x_cgs)
                    map_it(pos, hsml, mass, rho, rho_gcm3, weights, kernel=k, units="g/cm^2", param=par,
                        reduce_image=false, parallel=false;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_rho.xy.fits")

                    @test image ≈ rho_ref
                end

                @testset "T" begin
                    image_prefix = "sedov_T"
                    map_it(pos, hsml, mass, rho, T, rho, kernel=k, units="K", param=par,
                        reduce_image=true, parallel=false;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_T.xy.fits")

                    @test image ≈ T_ref
                end
            end


            @testset "Multi core" begin
            
                @testset "rho" begin
                    image_prefix = "sedov_rho"
                    weights = part_weight_physical(length(hsml), par, GU.x_cgs)
                    map_it(pos, hsml, mass, rho, rho_gcm3, weights, kernel=k, units="g/cm^2", param=par,
                        reduce_image=false, parallel=true;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_rho.xy.fits")

                    @test image ≈ rho_ref
                end

                @testset "T" begin
                    image_prefix = "sedov_T"
                    map_it(pos, hsml, mass, rho, T, rho, kernel=k, units="K", param=par,
                        reduce_image=true, parallel=true;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_T.xy.fits")

                    @test image ≈ T_ref
                end

            end

        end

        @testset "3D" begin

            # sedov reference snapshot
            fi = "snap_sedov"

            # read header and find units
            h = read_header(fi)
            GU = GadgetPhysical()

            blocks = ["POS", "MASS", "HSML", "RHO", "U"]
            data = Dict(block => read_block(fi, block, parttype=0) for block ∈ blocks)

            pos = data["POS"] .* GU.x_physical
            hsml = data["HSML"] .* GU.x_physical
            rho = data["RHO"] .* GU.rho_physical
            mass = data["MASS"] .* GU.m_physical

            kernel = WendlandC6(3)

            par = mappingParameters(center = [3.0, 3.0, 3.0],
                            x_size = 6.0, y_size = 6.0, z_size = 6.0,
                            Npixels = 10)
            
            @testset "Single core" begin
                d = sphMapping(pos, hsml, mass, rho, rho, ones(length(rho)),
                                    param=par, kernel=kernel,
                                    parallel = false,
                                    show_progress=true,
                                    dimensions=3)

                @test !isnan(d[1,1,1])
                @test d[1, 1, 1] ≈ 0.002470815089339027
            end

            pos = data["POS"] .* GU.x_physical

            @testset "Multi core" begin
                d = sphMapping(pos, hsml, mass, rho, rho, ones(length(rho)),
                                    param=par, kernel=kernel,
                                    parallel = true,
                                    show_progress=false,
                                    dimensions=3)

                @test !isnan(d[1, 1, 1])
                @test d[1, 1, 1] ≈ 0.002470815089339027
            end

        end # 3D

    end

    # @testset "TSC Mapping" begin

    #     fi = joinpath(dirname(@__FILE__), "x.txt")
    #     x = Float32.(copy(transpose(readdlm(fi))))

    #     fi = joinpath(dirname(@__FILE__), "bin_q.txt")
    #     bin_quantity = Float32.(readdlm(fi))


    #     par = mappingParameters(center = [3.0, 3.0, 3.0],
    #         x_size = 6.0, y_size = 6.0, z_size = 6.0,
    #         Npixels = 200,
    #         boxsize = 6.0)


    #     @testset "2D" begin 
    #         # d = sphMapping( x, bin_quantity, 
    #         #             param=par, show_progress=true)

    #         # @test !isnan(d[1,1])

    #         # @test_nowarn sphMapping( x, bin_quantity, 
    #         #                 param=par, show_progress=false)
    #     end

    #     @testset "3D" begin
    #         # @test_nowarn sphMapping( x, bin_quantity, 
    #         #                 param=par, show_progress=false,
    #         #                 dimensions=3)
    #     end

    # end

    @testset "HealPix mapping" begin
        
        snap_base = "snap_cutout_072"

        kernel = WendlandC4(2)

        h  = GadgetIO.read_header(snap_base)
        GU = GadgetPhysical(h, xH=0.752)

        Nside = 128
        center = 0.5h.boxsize .* ones(3) .* GU.x_physical

        hsml = read_block(snap_base, "HSML", parttype=0) .* GU.x_physical
        rho  = read_block(snap_base, "RHO", parttype=0)  .* GU.rho_physical
        mass = read_block(snap_base, "MASS", parttype=0) .* GU.m_physical

        T_K   = read_block(snap_base, "U", parttype=0) .* GU.T_K

        pos  = read_block(snap_base, "POS", parttype=0) .* GU.x_physical

        allsky_map, weight_map = healpix_map(pos, hsml, mass, rho, T_K, rho,
                            output_from_all_workers=true;
                            center, kernel, Nside)

        @inbounds for i ∈ eachindex(allsky_map)
            if !isnan(weight_map[i]) && !iszero(weight_map[i]) && !isinf(weight_map[i])
                allsky_map[i]  /= weight_map[i]
            end
        end

        # check if all pixels are filled
        @test length(findall(iszero.(allsky_map))) == 0

        # check min and max of map 
        @test minimum(allsky_map) ≈ 5312.211023105265
        @test maximum(allsky_map) ≈ 7.478005891665593e7

        # check sum of all pixels 
        @test sum(allsky_map[:]) ≈ 3.877477851438817e11

    end

    @testset "FITS io" begin

        # sedov reference snapshot
        fi = "snap_sedov"

        # kernel definition
        k = WendlandC4(2)

        # read header and find units
        h = read_header(fi)
        GU = GadgetPhysical(xH=0.752)

        blocks = ["POS", "MASS", "HSML", "RHO", "U"]
        data = Dict(block => read_block(fi, block, parttype=0) for block ∈ blocks)

        # convert to physical code units for mapping
        pos = data["POS"] .* GU.x_physical
        hsml = data["HSML"] .* GU.x_physical
        rho = data["RHO"] .* GU.rho_physical
        mass = data["MASS"] .* GU.m_physical

        T = data["U"] .* GU.T_K

        Npart = length(T)

        center = ones(3) .* 0.5h.boxsize * GU.x_physical
        xy_size = 0.9h.boxsize * GU.x_physical
        z_size = 0.9h.boxsize * GU.x_physical

        # define mapping parameters
        par = mappingParameters(center=center .* GU.x_physical,
            x_size=xy_size,
            y_size=xy_size,
            z_size=z_size,
            Npixels=256,
            boxsize=h.boxsize * GU.x_physical)

        d = sphMapping(pos, hsml, mass, rho, T, rho,
            param = par, kernel = k,
            parallel = false,
            show_progress = false)

        # store image in a file
        fits_file = joinpath(dirname(@__FILE__), "image.fits")

        @test_nowarn write_fits_image(fits_file, d, par)

        # read image back into memory and compare
        image, fits_par, snap = read_fits_image(fits_file)

        @test image ≈ d 
        @test par.boxsize == fits_par.boxsize
        @test par.center == fits_par.center
    end

    @testset "Reconstructing Grid" begin
        par = mappingParameters(center = [3.0, 3.0, 3.0],
            x_size = 6.0, y_size = 6.0, z_size = 6.0,
            Npixels = 200,
            boxsize = 6.0)

        @test_nowarn get_map_grid_2D(par)
        @test_nowarn get_map_grid_3D(par)
    end

    @testset "Weight functions" begin

        weight = part_weight_one(1)
        @test weight[1] == 1.0

        par = mappingParameters(center = [3.0, 3.0, 3.0],
                                x_size = 6.0, y_size = 6.0, z_size = 6.0,
                                Npixels = 200,
                                boxsize = 6.0)

        weight = part_weight_physical(1, par)
        @test weight[1] ≈ par.pixelSideLength * 3.085678e21

        weight = part_weight_physical(1)
        @test weight[1] ≈ 3.085678e21

        weight = part_weight_emission([0.5, 0.5], [0.5, 0.5])
        @test weight[1] ≈ 0.1767766952966369

        weight = part_weight_spectroscopic([0.5, 0.5], [0.5, 0.5])
        @test weight[1] ≈ 0.4204482076268573

        weight = part_weight_XrayBand([0.5, 0.5], 0.5, 1.5)
        @test weight[1] ≈ 0.0


    end

    @testset "Effect functions" begin

        @testset "Density" begin
            @test density_2D(1.0, 1.0) ≈ 6.769911178294544e-22
        end

        @testset "SZ-effect" begin
            @test SPHtoGrid.Tcmb(0.0) ≈ 2.728
            @test SPHtoGrid.Tcmb(10.0) ≈ 30.008000000000003

            @test kinetic_SZ([1.0], [1.0])[1] ≈ -2.2190366589946296e-35

            @test thermal_SZ([1.0], [1.0])[1] ≈ 3.876947202614958e-34
        end

        @testset "X-Ray" begin
            T_keV = [10.0]
            rho_cgs = [1.e-28]
            metalicity = [0.3]

            # in spectral range
            @test x_ray_emissivity(T_keV, rho_cgs, metalicity, E0=0.0, E1=Inf)[1] ≈ 3.092778619918078e-31
            # in spectral range
            @test x_ray_emissivity(T_keV, rho_cgs, metalicity)[1] ≈ 6.291391279343498e-32
            # cooling function without metals
            @test x_ray_emissivity(T_keV, rho_cgs, cooling_function=true)[1] ≈ 1.8440514708842878e-32
            # cooling function with metals
            @test x_ray_emissivity(T_keV, rho_cgs, metalicity, cooling_function=true)[1] ≈ 1.877220382159635e-32
            # cooling function with metals bolometric
            @test x_ray_emissivity(T_keV, rho_cgs, metalicity, cooling_function=true, E0=0.0, E1=Inf)[1] ≈ 5.741669221015086e-32
        end

        @testset "gamma" begin
            # set up reference values
            rho = 1.e-28
            T_K = 1.e7
            α_p = 2.235
            m_cgs = 1.989e40
            Mpc = 3.085678e24

            @test λγ_PE04(rho, T_K, α_p) ≈ 2.0962584497769528e-31
            @test jγ_PE04(rho, T_K, α_p, 1.0) ≈ 3.610788768607327e-32
            @test gamma_luminosity_pions_PE04(rho, m_cgs, T_K, α_p) ≈ 3.549481187521353e37
            @test gamma_flux_pions_PE04(rho, m_cgs, T_K, α_p, Mpc) ≈ 3.484725208525935e-13
        end

        @testset "Synchrotron" begin

            @testset "DSA model selection" begin
                # integer based selection
                @test SPHtoGrid.select_dsa_model(0) == Kang07()
                @test SPHtoGrid.select_dsa_model(1) == KR13()
                @test SPHtoGrid.select_dsa_model(2) == Ryu19()
                @test SPHtoGrid.select_dsa_model(3) == CS14()
                @test SPHtoGrid.select_dsa_model(4) == P16()
                # pass-through
                @test SPHtoGrid.select_dsa_model(Kang07()) == Kang07()
                # error handling 
                @test_throws ErrorException("Invalid DSA model selection!") SPHtoGrid.select_dsa_model(10)
            end

            @testset "Default" begin
                @test analytic_synchrotron([2.e-20], [5.0e-6], [3.0], dsa_model=1)[1] ≈ 6.606912895256343e-51
                @test analytic_synchrotron([2.e-20], [5.0e-6], [3.0], [π / 4], dsa_model=1)[1] ≈ 3.3069307682942575e-51
                @test analytic_synchrotron([2.e-20], [5.0e-6], [3.0], dsa_model=1, integrate_pitch_angle=false)[1] ≈ 9.556961692634717e-51
                @test analytic_synchrotron([2.e-20], [5.0e-6], [3.0], dsa_model=1, polarisation=true)[1] ≈ 4.784316234516861e-51
            end

            @testset "Ginzburg & Syrovatskii (1965)" begin
                @test analytic_synchrotron_GS([1.e-28], [5.0e-6], [1.0e8], [3.0], dsa_model=1)[1] ≈ 1.1397808537425788e-54
                @test analytic_synchrotron_GS([1.e-28], [5.0e-6], [1.0e8], [3.0], [π / 4], dsa_model=1)[1] ≈ 5.698904268712894e-55
            end

            @testset "Longair (2011)" begin
                @test analytic_synchrotron_Longair([1.e-28], [5.0e-6], [1.0e8], [3.0], dsa_model=1)[1] ≈ 1.4322908627279946e-53
                @test analytic_synchrotron_Longair([1.e-28], [5.0e-6], [1.0e8], [3.0], [π / 4], dsa_model=1)[1] ≈ 7.161454313639973e-54
            end

            @testset "Hoeft&Brüggen (2007)" begin
                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], dsa_model=0)[1] ≈ 8.710630889115241e-38

                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], dsa_model=1)[1] ≈ 4.6891522065905656e-39

                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], dsa_model=2)[1] ≈ 5.554619187930485e-39

                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], dsa_model=3)[1] ≈ 2.3445761032952828e-39

                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], dsa_model=4)[1] ≈ 3.907626838825512e-37

                @test analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
                    [5.0e-6], [8.618352059925092], [3.0], [π / 4])[1] ≈ 2.3445761032952828e-39
            end

    
        end
    end

    @testset "Image Functions" begin

        @testset "Surface Brightness to Luminosity" begin
            par = mappingParameters(center=[0.0, 0.0, 0.0],
                x_size=1.0,
                y_size=1.0,
                z_size=1.0,
                Npixels=10)

            jnu_image = 1.e-20 .* ones(10, 10)

            @test surface_brightness_to_luminosity(jnu_image, par)[1] ≈ 9.521408719683998e20
            @test surface_brightness_to_luminosity(jnu_image, par.pixelSideLength)[1] ≈ 9.521408719683998e20
        end

        @testset "Synchrotron Surface Brightness to Luminosity" begin
            par = mappingParameters(center=[0.0, 0.0, 0.0],
                x_size=1.0,
                y_size=1.0,
                z_size=1.0,
                Npixels=10)

            jnu_image = 1.e-20 .* ones(10, 10)

            @test synchrotron_SB_to_luminosity(jnu_image, par)[1] ≈ 9.521408719683998e13
            @test synchrotron_SB_to_luminosity(jnu_image, par.pixelSideLength)[1] ≈ 9.521408719683998e13
        end

        @testset "Total Synchrotron Luminosity" begin
            par = mappingParameters(center=[0.0, 0.0, 0.0],
                x_size=1.0,
                y_size=1.0,
                z_size=1.0,
                Npixels=10)

            jnu_image = 1.e-20 .* ones(10, 10)

            @test total_synch_luminosity_from_SB(jnu_image, par) ≈ 9.521408719683998e15

        end

        @testset "Synchrotron Polarisation" begin
            # set up reference images 
            Npixels = 128
            Q_image = zeros(Npixels, Npixels)
            U_image = Matrix{Float64}(undef, Npixels, Npixels)
            U_image .= 3.964929157902007e-28
            Iν_image = Matrix{Float64}(undef, Npixels, Npixels)
            Iν_image .= 5.4753081210232675e-28

            @testset "Fraction" begin
                Π = polarisation_fraction(Q_image, U_image, Iν_image)
                @test Π[1] .≈ 0.7241472206245468
            end

            @testset "Angle" begin
                ψ = polarisation_angle(Q_image, U_image, Iν_image)
                @test ψ[1] .≈ 45.0
            end
        end
    end
end


rm("snap_sedov")
rm("snap_cutout_072")
rm("image.fits")
rm("sedov_rho.xy.fits")
rm("sedov_T.xy.fits")