using Downloads

# @info "downloading test data..."
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_sedov", "./snap_sedov")

# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.0", "./snap_002.0")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.1", "./snap_002.1")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.2", "./snap_002.2")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.3", "./snap_002.3")

# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.0.key", "./snap_002.0.key")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.1.key", "./snap_002.1.key")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.2.key", "./snap_002.2.key")
# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.3.key", "./snap_002.3.key")

# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_002.key.index", "./snap_002.key.index")

# Downloads.download("http://www.usm.uni-muenchen.de/~lboess/SPHtoGrid/snap_cutout_072", "./snap_cutout_072")

# @info "done"

using Distributed
addprocs(2)

@everywhere using SPHtoGrid, Test, DelimitedFiles, SPHKernels, GadgetIO, GadgetUnits

@testset "SPHtoGrid" begin

    @testset "Smac utility" begin

        @test_nowarn write_smac1_par("./")


        # @test_throws ErrorException("Read error: Incorrect image format!") read_smac1_binary_image(joinpath(dirname(@__FILE__), "snap_050"))

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

        @test_throws ErrorException("Please specify pixelSideLenght or number of pixels!") mappingParameters(center = [0.0, 0.0, 0.0],
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

    @testset "SPH Mapping" begin

        @info "SPH Mapping tests take a while..."

        @info "Data read-in."

        # reference images
        T_ref, par, snap, units   = read_fits_image("sedov_T_reference.fits")
        rho_ref, par, snap, units = read_fits_image("sedov_rho_reference.fits")

        # sedov reference snapshot
        fi = "snap_050"

        # kernel definition
        k = WendlandC4(2)

        # read header and find units
        h = read_header(fi)
        GU = GadgetPhysical(xH=0.752)

        blocks = ["POS", "MASS", "HSML", "RHO", "U"]
        data = Dict(block => read_snap(fi, block) for block in blocks)

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

                    image, par, snap, units = read_fits_image("sedov_rho.fits")

                    @test image ≈ rho_ref
                end

                @testset "T" begin
                    image_prefix = "sedov_T"
                    map_it(pos, hsml, mass, rho, T, rho, kernel=k, units="T", param=par,
                        reduce_image=true, parallel=false;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_T.fits")

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

                    image, par, snap, units = read_fits_image("sedov_rho.fits")

                    @test image ≈ rho_ref
                end

                @testset "T" begin
                    image_prefix = "sedov_T"
                    map_it(pos, hsml, mass, rho, T, rho, kernel=k, units="T", param=par,
                        reduce_image=true, parallel=true;
                        snap, image_prefix)

                    image, par, snap, units = read_fits_image("sedov_T.fits")

                    @test image ≈ T_ref
                end

            end

        end

        # @testset "3D" begin

        #     kernel = WendlandC6(3)

        #     par = mappingParameters(center = [3.0, 3.0, 3.0],
        #                     x_size = 6.0, y_size = 6.0, z_size = 6.0,
        #                     Npixels = 10,
        #                     boxsize = 6.0)
            
        #     @testset "Single core" begin
        #         d = sphMapping(x, hsml, m, rho, bin_quantity, ones(Float32, size(rho,1)),
        #                             param=par, kernel=kernel,
        #                             parallel = false,
        #                             show_progress=true,
        #                             dimensions=3)

        #         @test !isnan(d[1,1,1])
        #     end

        #     @testset "Multi core" begin
        #         # @test_nowarn sphMapping(x, hsml, m, rho, bin_quantity, ones(Float32, size(rho,1)),
        #         #                     param=par, kernel=kernel,
        #         #                     parallel = true,
        #         #                     show_progress=false,
        #         #                     dimensions=3)
        #     end

        # end # 3D

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

    # @testset "HealPix mapping" begin
        
    #     snap_base = "snap_cutout_072"

    #     kernel = WendlandC4(2)

    #     h  = GadgetIO.read_header(snap_base)
    #     GU = GadgetPhysical(h, xH=0.752)

    #     Nside = 128
    #     center = 0.5h.boxsize .* ones(3) .* GU.x_physical

    #     hsml = read_block(snap_base, "HSML", parttype=0) .* GU.x_physical
    #     rho  = read_block(snap_base, "RHO", parttype=0)  .* GU.rho_physical
    #     mass = read_block(snap_base, "MASS", parttype=0) .* GU.m_physical

    #     T_K   = read_block(snap_base, "U", parttype=0) .* GU.T_K

    #     pos  = read_block(snap_base, "POS", parttype=0) .* GU.x_physical

    #     allsky_map, weight_map = healpix_map(pos, hsml, mass, rho, T_K, rho;
    #                         center, kernel, Nside)

    #     @inbounds for i ∈ eachindex(allsky_map)
    #         if !isnan(weight_map[i]) && !iszero(weight_map[i]) && !isinf(weight_map[i])
    #             allsky_map[i]  /= weight_map[i]
    #         end
    #     end

    #     # check if all pixels are filled
    #     @test length(findall(iszero.(allsky_map))) == 0

    #     # check min and max of map 
    #     @test minimum(allsky_map) ≈ 4766.125101327221
    #     @test maximum(allsky_map) ≈ 5.937353111693359e7

    #     # check sum of all pixels 
    #     @test sum(allsky_map[:]) ≈ 3.637350200286616e11

    # end

    # @testset "FITS io" begin

    #     # map data
    #     fi = joinpath(dirname(@__FILE__), "bin_q.txt")
    #     bin_quantity = Float32.(readdlm(fi))[:, 1]

    #     fi = joinpath(dirname(@__FILE__), "x.txt")
    #     x = copy(transpose(Float32.(readdlm(fi))))

    #     fi = joinpath(dirname(@__FILE__), "rho.txt")
    #     rho = Float32.(readdlm(fi))[:, 1]

    #     fi = joinpath(dirname(@__FILE__), "hsml.txt")
    #     hsml = Float32.(readdlm(fi))[:, 1]

    #     fi = joinpath(dirname(@__FILE__), "m.txt")
    #     m = Float32.(readdlm(fi))[:, 1]

    #     kernel = WendlandC6(2)

    #     par = mappingParameters(center = [3.0, 3.0, 3.0],
    #         x_size = 6.0, y_size = 6.0, z_size = 6.0,
    #         Npixels = 200,
    #         boxsize = 6.0)

    #     d = sphMapping(x, hsml, m, rho, bin_quantity, ones(Float32, size(rho, 1)),
    #         param = par, kernel = kernel,
    #         parallel = false,
    #         show_progress = false)

    #     # store image in a file
    #     fits_file = joinpath(dirname(@__FILE__), "image.fits")

    #     @test_nowarn write_fits_image(fits_file, d, par)

    #     # read image back into memory and compare
    #     # image, fits_par, snap = read_fits_image(fits_file)

    #     # @test image ≈ d 
    #     # @test par.boxsize == fits_par.boxsize
    #     # @test par.center == fits_par.center


    # end

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

        weight = part_weight_emission([0.5, 0.5], [0.5, 0.5])
        @test weight[1] ≈ 0.1767766952966369

        weight = part_weight_spectroscopic([0.5, 0.5], [0.5, 0.5])
        @test weight[1] ≈ 0.4204482076268573

        weight = part_weight_XrayBand([0.5, 0.5], 0.5, 1.5)
        @test weight[1] ≈ 0.0


    end

    # @testset "Effect functions" begin

    #     @testset "Density" begin
    #         @test density_2D(1.0, 1.0) ≈ 6.769911178294544e-22
    #     end

    #     @testset "SZ-effect" begin
    #         @test SPHtoGrid.Tcmb(0.0) ≈ 2.728
    #         @test SPHtoGrid.Tcmb(10.0) ≈ 30.008000000000003

    #         @test kinetic_SZ([1.0], [1.0])[1] ≈ -2.2190366589946296e-35

    #         @test thermal_SZ([1.0], [1.0])[1] ≈ 3.876935843260665e-34
    #     end

    #     @testset "X-Ray" begin
    #         # in spectral range
    #         @test x_ray_emission([10.0], [1.989e38], [1.e-28])[1] ≈ 1.9479441169462751e34
    #         # bolometric
    #         @test x_ray_emission([10.0], [1.989e38], [1.e-28], E0=0.0, E1=Inf)[1] ≈ 9.575878609659925e34
    #     end

    #     @testset "Synchrotron" begin

    #         @testset "Analytic" begin
    #             @test analytic_synchrotron_emission([1.0], [1.0], [1.0], [10.0])[1] ≈ 6.424386277144697e-25

    #             @test analytic_synchrotron_emission([1.0], [1.0], [1.0], [10.0], convert_to_mJy = true)[1] ≈ 6.424386277144697e1

    #             @test analytic_synchrotron_emission([1.0], [1.0], [1.0], [1.0])[1] == 0.0

    #             @test_throws ErrorException("Invalid DSA model selection!") analytic_synchrotron_emission([1.0], [1.0], [1.0], [1.0], dsa_model = 10)
    #         end

    #         @testset "Spectrum" begin

    #             # # test thermal energy density - redundant!
    #             # # acc_function = SPHtoGrid.KR13_acc
    #             # # ϵ_th = SPHtoGrid.EpsNtherm(1.52606e-30, 4.75088e+08, xH=0.76)
    #             # # ϵ_cr0 = 0.01 * SPHtoGrid.get_rel_energy_density(10.0, acc_function) * ϵ_th
    #             # # @test ϵ_cr0 / ϵ_th ≈ 0.00244270822665505

    #             # @testset "ϵ_th" begin
    #             #     # with pitch angle integration
    #             #     j_ν = spectral_synchrotron_emission(1.52606e-30, 5.0e-6, 4.75088e+08, 5.0, dsa_model = 2,
    #             #         integrate_pitch_angle = true)

    #             #     @test j_ν ≈ 1.8643638246341477e-55

    #             #     # without pitch angle integration
    #             #     j_ν = spectral_synchrotron_emission(1.52606e-30, 5.0e-6, 4.75088e+08, 5.0, dsa_model = 2,
    #             #         integrate_pitch_angle = false)

    #             #     @test j_ν ≈ 4.451252238469209e-55

    #             #     # conversion to mJy/cm
    #             #     j_ν = spectral_synchrotron_emission(1.52606e-30, 5.0e-6, 4.75088e+08, 5.0, dsa_model = 2,
    #             #         integrate_pitch_angle = false,
    #             #         convert_to_mJy = true)

    #             #     @test j_ν ≈ 4.4512522384692095e-29

    #             # end

    #             # @testset "pre-defined" begin
    #             #     # test for pre-defined spectrum
    #             #     Nbins = 128
    #             #     q0 = 4.166666666666667

    #             #     # define spectrum
    #             #     bounds = 10.0 .^ LinRange(-1.0, 6.0, Nbins + 1)
    #             #     norm = Vector{Float64}(undef, Nbins)
    #             #     # reference from previous test
    #             #     norm[1] = 2.3141104241756675e-30
    #             #     for Nbin = 2:Nbins-1
    #             #         norm[Nbin] = norm[Nbin-1] * (bounds[Nbin] / bounds[Nbin-1])^(-q0)
    #             #     end

    #             #     j_ν = spectral_synchrotron_emission(norm, bounds, 5.e-6,
    #             #         integrate_pitch_angle = true)

    #             #     @test j_ν ≈ 1.6530988606617404e-55

    #             #     j_ν = spectral_synchrotron_emission(norm, bounds, 5.e-6,
    #             #         integrate_pitch_angle = false)

    #             #     @test j_ν ≈ 4.2742001929052135e-55

    #             #     j_ν = spectral_synchrotron_emission(norm, bounds, 5.e-6,
    #             #         integrate_pitch_angle = false,
    #             #         convert_to_mJy = true)

    #             #     @test j_ν ≈ 4.2742001929052134e-29
    #             # end # pre-fedined
    #         end # Spectrum
    #     end
    # end
end


# rm("snap_sedov")

# rm("snap_002.0")
# rm("snap_002.1")
# rm("snap_002.2")
# rm("snap_002.3")
# rm("snap_002.0.key")
# rm("snap_002.1.key")
# rm("snap_002.2.key")
# rm("snap_002.3.key")
# rm("snap_002.key.index")

# rm("snap_cutout_072")