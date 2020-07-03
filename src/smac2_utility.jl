
"""
    write_smac2_par([...])

Writes a P-Smac2 parameter file. Not all relevant parameters implemented yet!
"""
function write_smac2_par(x, y, z,
                         euler_angle_0, euler_angle_1, euler_angle_2,
                         xy_size, z_depth, xy_pix::Int64,
                         input_file, output_file, path,
                         effect_module::Int64=0, effect_flag::Int64=0,
                         CR_pmin::Float64=10.0, CR_pmax::Float64=1.e7)

    open(path * "smac2.par", "w") do f
        write(f,
        """%% PSmac Parameter File %%

        Cosmology 0         % set to 0 for non-cosmological sims.  Cosmo.h = 1

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

        Freq 1.4e9          % Frequency [Hz]: Sz,Synchro

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
