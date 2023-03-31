using SnoopPrecompile    # this is a small dependency

@precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.

    kernels = [WendlandC4(2), WendlandC6(2)]

    center = [247.980, 245.480, 255.290] .* 1.e3
    hp_pos = Matrix{Float64}(undef, 3, 7)
    hp_pos[:, 1] = [225160.875, 256677.5625, 242031.765625] # Centaurus
    hp_pos[:, 2] = [245040.453125, 327781.84375, 246168.6875] # Coma
    hp_pos[:, 3] = [252454.0, 238890.296875, 233772.890625] # Fornax
    hp_pos[:, 4] = [202408.359375, 245067.234375, 249614.09375] # Norma
    hp_pos[:, 5] = [176696.22, 266631.8, 315976.38] # Ophiuchus
    hp_pos[:, 6] = [307184.78125, 247627.078125, 230736.484375] # Perseus
    hp_pos[:, 7] = [244450.578125, 255851.78125, 253668.1875] # Virgo

    hp_hsml = [1771.7005615234375, 2177.5625, 714.1318359375, 798.90478515625,
        1067.0657, 1813.0419921875, 1711.311767578125]

    hp_mass = ones(7)
    hp_rho = ones(7)

    Nside = 128

    cic_pos = 15.0 .* (rand(3, 100) .- 0.5)
    cic_hsml = 2.0 .* rand(100)
    cic_mass = rand(100)
    cic_rho = rand(100)
    cic_T = 1.e8 .* rand(100)
    snap = 0

    image_prefix = "dummy"
    # define mapping parameters
    par = mappingParameters(center=[0.0, 0.0, 0.0],
        x_size=10.0,
        y_size=10.0,
        z_size=10.0,
        Npixels=256)
    

    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

        @info "Pre-compiling mapping functions"

        for kernel ∈ kernels 
            # cic 
            map_it(cic_pos, cic_hsml, cic_mass, cic_rho, cic_T, cic_rho, units="T", param=par,
                reduce_image=true, parallel=false, show_progress=true;
                snap, image_prefix, kernel)

            # healpix
            healpix_map(hp_pos, hp_hsml, hp_mass, hp_rho, hp_rho, hp_rho, show_progress=true;
                center, kernel, Nside)

            # delete dummy file 
            rm("dummy.fits")
        end
    end
end