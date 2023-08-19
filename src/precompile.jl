using PrecompileTools    # this is a small dependency

@setup_workload begin
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
    
    Npixels = 128
    Q_image = zeros(Npixels, Npixels)
    U_image = Matrix{Float64}(undef, Npixels, Npixels)
    U_image .= 3.964929157902007e-28
    Iν_image = Matrix{Float64}(undef, Npixels, Npixels)
    Iν_image .= 5.4753081210232675e-28

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

        for kernel ∈ kernels 
            # cic 
            map_it(cic_pos, cic_hsml, cic_mass, cic_rho, cic_T, cic_rho, units="T", param=par,
                reduce_image=true, parallel=false, show_progress=false;
                snap, image_prefix, kernel)

            # healpix
            healpix_map(hp_pos, hp_hsml, hp_mass, hp_rho, hp_rho, hp_rho, show_progress=false;
                center, kernel, Nside)

        end

        # effects 
        kinetic_SZ([1.0], [1.0])
        thermal_SZ([1.0], [1.0])
        x_ray_emissivity([10.0], [1.e-28])
        gamma_luminosity_pions_PE04(1.e-28, 1.989e40, 1.e7, 2.235)
        gamma_flux_pions_PE04(1.e-28, 1.989e40, 1.e7, 2.235, 3.085678e24)
        analytic_synchrotron([2.e-20], [5.0e-6], [3.0], dsa_model=1)
        analytic_synchrotron([2.e-20], [5.0e-6], [3.0], [π / 4], dsa_model=1)
        analytic_synchrotron_HB07([1.e-28], [1.9890000000000002e39], [6.171355999999999e22],
            [5.0e-6], [8.618352059925092], [3.0], dsa_model=1)

        # image functions 
        polarisation_fraction(Q_image, U_image, Iν_image)
        polarisation_angle(Q_image, U_image, Iν_image)
    end

    if isfile("dummy.xy.fits")
        # delete dummy file 
        rm("dummy.xy.fits")
    end
end