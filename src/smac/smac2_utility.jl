
"""
    write_smac2_par(x, y, z,
                    euler_angle_0, euler_angle_1, euler_angle_2,
                    xy_size, z_depth, xy_pix::Integer,
                    input_file, output_file, path,
                    effect_module::Integer=0, effect_flag::Integer=0,
                    ν_obs::Real=1.44e9, cosmology::Integer=0,
                    CR_pmin::Real=10.0, CR_pmax::Real=1.e7)

Writes a P-Smac2 parameter file. Not all relevant parameters implemented yet!
"""
function write_smac2_par(x, y, z,
                         euler_angle_0, euler_angle_1, euler_angle_2,
                         xy_size, z_depth, xy_pix::Integer,
                         input_file, output_file, path,
                         effect_module::Integer=0, effect_flag::Integer=0,
                         ν_obs::Real=1.44e9, cosmology::Integer=0,
                         CR_pmin::Real=10.0, CR_pmax::Real=1.e7)

    open(path * "smac2.par", "w") do f
        write(f,
        """%% PSmac Parameter File %%

        Cosmology $cosmology         % set to 0 for non-cosmological sims.  Cosmo.h = 1

        %% Image Properties, comoving %%

        Center_X $x
        Center_Y $y
        Center_Z $z

        Use_Barycenter 0    % 1 - set to, 2 - set relative to barycenter

        XYSize $xy_size        % set these two ==0 to project whole box
        ZDepth $z_depth
        XYPix $xy_pix

        %% Projection, see description at EOF %%
        Euler_Angle_0 $euler_angle_0    % [deg] along x : (0,90,90)
        Euler_Angle_1 $euler_angle_1    % [deg] along y : (0,90,0)
        Euler_Angle_2 $euler_angle_2     % [deg] along z : (0,0,0)

        %% HEALPIX %%

        NSide 32
        Rmin 1000
        Rmax 1000000

        %% I/O Options %%

        Input_File $input_file
        N_IOTasks 1

        NoClobber 0
        Output_File $output_file


        %% Effect Selection & Options %%

        Effect_Module $effect_module
        Effect_Flag $effect_flag

        E_min 5e4             % Energy Range [eV]: Xray, Gamma, Synchro
        E_max 1e10

        Freq $ν_obs          % Frequency [Hz]: Sz,Synchro

        a_cr 2.375          % CR spectral index
        X_cr 0.01           % CR Normalisation rel. thermal
        IntrRM 0            % Toggle Intrinsic RM
        PitchAngInt 1       % Toggle Pitch angle integration

        t_turb 1e7          % Reacceleration Timescale Cassano 05
        eta_t 0.26          % Fraction turb. Energy in magnetosonic waves
        turb_scale 100      % scale to extrapolate turbulent velocity to

        CR_Emin $CR_pmin       % Minimum energy for BP_REAL_CRs, as in OpenGadget3 parameterfile
        CR_Emax $CR_pmax       % Maximum energy for BP_REAL_CRs, as in OpenGadget3 parameterfile

        %Gadget
        UnitLength_in_cm 			3.085678e21        %  1.0 kpc
        UnitMass_in_g 				1.989e43           %  1.0e10 solar masses
        UnitVelocity_in_cm_per_s 	1e5                %  1 km/sec
        """
        )
    end
end




"""
    read_smac2_info(filename::String)

Returns the image info in a `mappingParameters` struct.
"""
function read_smac2_info(filename::String)

    f = FITS(filename)

    center = [read_key(f[1], 23)[2], 
              read_key(f[1], 24)[2], 
              read_key(f[1], 25)[2] ]

    xy_size = read_key(f[1], 27)[2]
    z_size  = read_key(f[1], 28)[2]
    Npixels = read_key(f[1], 29)[2]

    close(f)

    return mappingParameters(center=center, 
									x_size=xy_size,
									y_size=xy_size,
									z_size=z_size,
									Npixels=Npixels)

end


"""
    read_smac2_image(filename::String)

Returns the image of a Smac2 FITS file.
"""
function read_smac2_image(filename::String, num_image::Integer=1)

    f = FITS(filename)

    image = read(f[1])[:,:,num_image]

    close(f)

    return image

end