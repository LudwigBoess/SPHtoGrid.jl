"""
        Helper function for Smac reading of binary images
"""

"""
    struct Smac1ImageInfo

Stores the information in a Smac binary image header.
"""
struct Smac1ImageInfo

    snap::Int32                 # number of input snapshot
    z::Float32                  # redshift of snapshot
    m_vir::Float32              # virial mass of halo
    r_vir::Float32              # virial radius of halo
    xcm::Float32                # x coordinate of image center
    ycm::Float32                # y coordinate of image center
    zcm::Float32                # z coordinate of image center
    z_slice_kpc::Float32        # depth of the image in kpc
    boxsize_kpc::Float32        # xy-size of the image in kpc
    boxsize_pix::Int32        # xy-size of the image in pixels
    pixsize_kpc::Float32        # size of one pixel in kpc
    xlim::Array{Float64,1}      # x limits of image
    ylim::Array{Float64,1}      # y limits of image
    zlim::Array{Float64,1}      # z limits of image
    units::String               # unitstring of image

    function Smac1ImageInfo(snap::Int32, z::Float32, m_vir::Float32, r_vir::Float32,
                       xcm::Float32, ycm::Float32, zcm::Float32,
                       z_slice_kpc::Float32,
                       boxsize_kpc::Float32, boxsize_pix::Int32, pixsize_kpc::Float32,
                       units::String)

        xlim = [xcm - boxsize_kpc/2.0, xcm + boxsize_kpc/2.0]
        ylim = [ycm - boxsize_kpc/2.0, ycm + boxsize_kpc/2.0]
        zlim = [zcm - z_slice_kpc/2.0, zcm + z_slice_kpc/2.0]

        new(snap, z, m_vir, r_vir,
            xcm, ycm, zcm,
            z_slice_kpc,
            boxsize_kpc, boxsize_pix, pixsize_kpc,
            xlim, ylim, zlim, units)
    end

end

"""
    read_smac1_binary_image(filename::String)

Reads a binary image file from Smac and returns a Matrix with the pixel values.
"""
function read_smac1_binary_image(filename::String)
    f = open(filename)
    first = read(f, Int32)
    n_pixels = Int64(sqrt.(first/4.0))
    a = read!(f, Array{Float32,2}(undef,n_pixels,n_pixels))
    last = read(f, Int32)
    close(f)
    if first == last
        return collect(transpose(a))
    else
        error("Read error: Incorrect image format!")
    end
end

"""
    read_smac1_binary_info(filename::String)

Returns the image info in a `Smac1ImageInfo` struct.
"""
function read_smac1_binary_info(filename::String)

    f = open(filename)

    # skip image
    image = read(f, Int32)
        seek(f, position(f) + image + 8)

    # read info block
    snap = read(f, Int32)
        seek(f, position(f) + 8)
    z = read(f, Float32)
        seek(f, position(f) + 8)
    m_vir = read(f, Float32)
        seek(f, position(f) + 8)
    r_vir = read(f, Float32)
        seek(f, position(f) + 8)
    xcm = read(f, Float32)
        seek(f, position(f) + 8)
    ycm = read(f, Float32)
        seek(f, position(f) + 8)
    zcm = read(f, Float32)
        seek(f, position(f) + 8)
    z_slice_kpc = read(f, Float32)
        seek(f, position(f) + 8)
    boxsize_kpc = read(f, Float32)
        seek(f, position(f) + 8)
    boxsize_pix = read(f, Int32)
        seek(f, position(f) + 8)
    pixsize_kpc = read(f, Float32)
        seek(f, position(f) + 4)
    # read length of unit string
    unitstring_length = read(f, Int32)
    # read unti string
    unit = String(Char.((read!(f, Array{Int8,1}(undef,unitstring_length)))))
    close(f)

    return Smac1ImageInfo(snap, z, m_vir, r_vir,
                    xcm, ycm, zcm,
                    z_slice_kpc,
                    boxsize_kpc, boxsize_pix, pixsize_kpc,
                    unit)
end


"""
    write_smac1_par( path="./"; kwargs...)

Writes a Smac parameter file. Add keyword arguments according to the smac1 parameter file.
"""
function write_smac1_par( path="./"; kwargs...)

    open(path * "smac1.par", "w") do f
        write(f,
        """#=============================================================
        # Input file
        #=============================================================
        
        #**** Flag of the data format (integer - reals):
        #         (1) Navarro,
        #         (2) GADGET
        #             intial redshift, (case (1))
        #             omega matter,    (case (1))
        #             Hubble parameter (case (1))
        #             image redshift (set to 0. if obtained from the input data)
        #        (3) TOY-MODELS
        #        (4) GADGET2
        SIM_FORMAT = $( haskey(kwargs, :SIM_FORMAT) ? kwargs[:SIM_FORMAT] : 4 )
        
        #**** Input file
        USE_KEYS = $( haskey(kwargs, :USE_KEYS) ? kwargs[:USE_KEYS] : 0 )
        OUTPUT_DIR = $( haskey(kwargs, :OUTPUT_DIR) ? kwargs[:OUTPUT_DIR] : "" )
        
        SNAP_BASE = $( haskey(kwargs, :SNAP_BASE) ? kwargs[:SNAP_BASE] : "snap_" )
        SNAP_FILE = $( haskey(kwargs, :SNAP_FILE) ? kwargs[:SNAP_FILE] : "snap_" )
        
        
        #**** Halo ID
        HALO_ID = $( haskey(kwargs, :HALO_ID) ? kwargs[:HALO_ID] : "a" )
        
        #**** Flag of the loop on the input files (integers):
        #        disabled if SNAP_START = -1 else enabled
        #        first input
        #        last input  (=first_input if not given)
        SNAP_START = $( haskey(kwargs, :SNAP_START) ? kwargs[:SNAP_START] : -1 )
        SNAP_END  = $( haskey(kwargs, :SNAP_END) ? kwargs[:SNAP_END] : -1 )
        
        #**** Intial redshift and image redshift (set IMAGE_Z=0. to obtain it from the input data)
        INITIAL_Z = $( haskey(kwargs, :INITIAL_Z) ? kwargs[:INITIAL_Z] : 0. )
        IMAGE_Z   = $( haskey(kwargs, :IMAGE_Z) ? kwargs[:IMAGE_Z] : 0. )
        
        
        #=============================================================
        # KERNEL TO USE
        #=============================================================
        # Flag of kernel used for mapping (integer):
        #         (0) cubic-spline
        #         (1) quintic-spline
        #         (2) WendlandC4
        #         (3) WendlandC6
        KERNEL_TYPE = $( haskey(kwargs, :KERNEL_TYPE) ? kwargs[:KERNEL_TYPE] : 3 )
        
        #=============================================================
        # DISTRIBUTION SCHEME
        #=============================================================
        
        #**** Flag of particle distribution scheme on the image (integer):
        #         (0) NGC       ->  not implemented
        #         (1) CIC       ->  2D-map
        #         (2) HealPix   ->  full-sky
        #         (3) TSC       ->  2D-map
        PART_DISTR = $( haskey(kwargs, :PART_DISTR) ? kwargs[:PART_DISTR] : 1 )
        
        #=============================================================
        # 2D-MAP SETUP
        #=============================================================
        
        #**** Flag of the units of the side of image (integer-real):
        #         (1) arcmin
        #         (2) kpc
        #         side of image in the chosen units
        IMG_XY_UNITS = $( haskey(kwargs, :IMG_XY_UNITS) ? kwargs[:IMG_XY_UNITS] : 2 )
        IMG_XY_SIZE  = $( haskey(kwargs, :IMG_XY_SIZE) ? kwargs[:IMG_XY_SIZE] : 1000.0 )
        
        #**** Flag of the third dimension of the image:
        #        (1) cubic,  same value of the preciding point
        #	 (2) if the value is chosen by the man
        #         side of the third dimension in kpc (case(2))
        #	  minimum filling factor of the carrot [0<->100] (its a percentile)
        IMG_Z_UNITS = $( haskey(kwargs, :IMG_Z_UNITS) ? kwargs[:IMG_Z_UNITS] : 2 )
        IMG_Z_SIZE  = $( haskey(kwargs, :IMG_Z_SIZE) ? kwargs[:IMG_Z_SIZE] : 1000.0 )
        FILL_CARROT = $( haskey(kwargs, :FILL_CARROT) ? kwargs[:FILL_CARROT] : 0.0 )
        
        #**** Flag of the number of image side pixels (integers):
        #         nr of pixel
        #         if (=-1) -> given by smoothing
        IMG_SIZE = $( haskey(kwargs, :IMG_SIZE) ? kwargs[:IMG_SIZE] : 1024 )
        
        #**** Flag of the smoothing (FWHM) angle (integer - real):
        #         (0) no smoothing
        #         (1) in arcmin
        #         (2) in kpc
        #         (3) given by npart
        #       smoothing lenght in the choosen units (cases (1), (2))
        SMOOTH_UNITS = $( haskey(kwargs, :SMOOTH_UNITS) ? kwargs[:SMOOTH_UNITS] : 0 )
        SMOOTH_FWHM  = $( haskey(kwargs, :SMOOTH_FWHM) ? kwargs[:SMOOTH_FWHM] : 5.0 )
        
        #=============================================================
        # FULL-SKY MAP SETUP
        #=============================================================
        
        #**** nside HEALPix variable (Mauro)
        #           This variable must be a power of 2 and <= 8192
        #           nside = ( 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192)
        NSIDE = $( haskey(kwargs, :NSIDE) ? kwargs[:NSIDE] : 8 )
        
        #**** min/max radius of sphere arround x0,y0,z0 to be used [kpc]
        MIN_DIST = $( haskey(kwargs, :MIN_DIST) ? kwargs[:MIN_DIST] : 5000.0 )
        MAX_DIST = $( haskey(kwargs, :MAX_DIST) ? kwargs[:MAX_DIST] : 420000.0 )
        
        #=============================================================
        # OUTPUT EFFECT
        #=============================================================
        
        #**** Flag of the type of output map (integers):
        #         (0)  3D electron density
        #         (1)  2D elctron Density
        #         (2)  Temperature
        #           (0) mass weighted
        #           (1) Emission weighted temperature (e.g. former (3))
        #               (set also Energy cutoff bands below!)
        #           (2) Emission weighted temperature (with tabulated cooling function)
        #           (3) Spectroscopic (e.g. former (13))
        #           (4) Emission temperature (by Elena) (e.g. former (14))
        #         (4)  Volume of particles
        #         (5)  age of particles
        #         (6)  X-ray SB of particles
        #           (0) simple sqrt(T)
        #              or  with cooling function from tables
        #           (1)       - energy band [0.1-10] keV (bolometric)
        #           (2)       - energy band [0.1-2.4] keV
        #           (3)       - energy band [0.3-3.5] keV
        #           (4)       - user defined energy band [ XRAY_E0 - XRAY_E1 ] keV (define below)
        #         (7)  Compton y-parameter: tSZ effect
        #             flag of the subtype of output map:
        #           (0) adimensional [DI/I]
        #           (1) adimensional [DT/T]
        #           (2) adimensinal [DT/T] with relativistic correction
        #           (3) n_e T^2, e.g. (1)*T
        #         (8)  Kinematical Sunjaev Zeldovich Effect of particles: kSZ efect
        #             flag of the subtype of output map:
        #           (0) adimensional [DI/I]
        #           (1) adimensional [DT/T]
        #         (9)  Q & U Stokes parameters for transvers particles motion: kpSZ effect
        #	 	      flag of the subtype of output map:
        #           (0) adimensional [Q_I & U_I] - noramlized on the CMB intensity
        #           (1) adimensional [Q_T & U_T] - noramlized on the CMB temperature
        #        (10) Mass weighted velocity
        #        (11) Mass weighted squared velocity
        #        (12) Magnetic field related maps
        #            flag of the subtype of output map:
        #           (0) RM
        #           (1) <B_x,y,z>
        #           (2) <|B|>
        #           (3) Synchrotron Emission
        #          (31) Synchrotron Polarised Emission & Polarisation Angle
        #          (10) RM with relativistic corrections
        #          (40) Gamma Emission [# of photons / cm^2] @ E_gam
        #          (41) Gamma Emission [erg/cm^2  ]			 @ E_gam
        #          (42) Gamma Emission [# of photons / cm^2] integrated over [E_pi_min,E_pi_max]
        #        (15) Metalicity map (of METAL_SPECIES)
        #            flag of the subtype of output map:
        #           (0) Mass weighted
        #           (1) Emission weighted
        #        (16) Rees-Sciama Effect by mass transverse motion [phi_dot]
        #        (17) Mach-Number
        #            flag of the subtype of output map:
        #           (0) Volume weighted
        #           (1) Energy weighted
        #           (2) Mass weighted
        #        (18) Y_X parameter (mass weighted)
        #        (19) Y_X parameter (pasquale weighted)
        #        (20) Projected potential
        #        (21) BP_CR Protons
        #           (0) Pressure
        #           (1) X_crp ( P_crp/P_th )
        #        (22) BP_CR Electrons
        #           (0) Pressure
        #           (1) Number density
        #           (2) Energy density
        #           (3) Synchrotron Emission
        #          (31) Synchrotron Polarised Emission & Polarisation Angle
        #       (100) 3D DM density
        #       (101) 2D DM density
        #       (102) Rees-Sciama Effect by mass transverse motion [projected momentum]
        #             flag of the subtype of output map:
        #           (0) by integral
        #           (1) by FFT
        #       (103) Gravitational lensing: Deflection angle
        OUTPUT_MAP = $( haskey(kwargs, :OUTPUT_MAP) ? kwargs[:OUTPUT_MAP] : 0 )
        OUTPUT_SUB = $( haskey(kwargs, :OUTPUT_SUB) ? kwargs[:OUTPUT_SUB] : 0 )
        
        # Define metal species
        #	 	flag of the subtype of output map:
        #		(0) total
        #               (1) > Fe
        #               (2) C
        #               (3) N
        #               (4) O
        #               (5) Mg
        #               (6) Si
        #               (7) Fe
        METAL_SPECIES = $( haskey(kwargs, :METAL_SPECIES) ? kwargs[:METAL_SPECIES] : 0 )
        
        #**** Unit of output Temperatures
        #        (0) Kelvin
        #        (1) KeV
        TEMP_UNIT = $( haskey(kwargs, :TEMP_UNIT) ? kwargs[:TEMP_UNIT] : 0 )
        
        #**** Ignore particles below this value [in keV]
        TEMP_CUT = $( haskey(kwargs, :TEMP_CUT) ? kwargs[:TEMP_CUT] : 0.0 )
        
        #**** Energy cut-off for X-ray (Flag of the type = 6):
        #        ene_a, (case(6), case(2)-subtype 1) lower energy bound for X-ray image [keV]
        #        ene_b, (case(6), case(2)-subtype 1) upper energy bound for X-ray image [keV]
        XRAY_E0 = $( haskey(kwargs, :XRAY_E0) ? kwargs[:XRAY_E0] : 0.5 )
        XRAY_E1 = $( haskey(kwargs, :XRAY_E1) ? kwargs[:XRAY_E1] : 2.0 )
        
        # Path to the tables for the cooling function
        X_TABLES = $( haskey(kwargs, :X_TABLES) ? kwargs[:X_TABLES] : "./Tables" )
        
        #**** Input observational frequency for SZ maps [GHz] (real):
        FRQ_ZS = $( haskey(kwargs, :FRQ_ZS) ? kwargs[:FRQ_ZS] : 300.0 )
        
        #**** Input for Radio Maps
        #XCRP_TABLE_PATH gives the change of XCRP over radius.
        ECRP_MIN = $( haskey(kwargs, :ECRP_MIN) ? kwargs[:ECRP_MIN] : 1.0e9 )
        FRQ_Pnu = $( haskey(kwargs, :FRQ_Pnu) ? kwargs[:FRQ_Pnu] : 1.4e9 )
        XCRP = $( haskey(kwargs, :XCRP) ? kwargs[:XCRP] : 0.01 )
        XCRP_TABLE_PATH = $( haskey(kwargs, :XCRP_TABLE_PATH) ? kwargs[:XCRP_TABLE_PATH] : "" )
        B_TABLE_PATH = $( haskey(kwargs, :B_TABLE_PATH) ? kwargs[:B_TABLE_PATH] : "" )
        KERNEL_TABLE_PATH = $( haskey(kwargs, :KERNEL_TABLE_PATH) ? kwargs[:KERNEL_TABLE_PATH] : "" )
        GAM_nu = $( haskey(kwargs, :GAM_nu) ? kwargs[:GAM_nu] : 1.25 )
        #XCRP_TABLE_FILE_FMT = (1F8.6,2X,1F14.6)
        
        
        #**** Input for Gamma Maps
        #energies in GeV
        E_gam = $( haskey(kwargs, :E_gam) ? kwargs[:E_gam] : 1.0e1 )
        E_gam_min = $( haskey(kwargs, :E_gam_min) ? kwargs[:E_gam_min] : 1.e-2 )
        E_gam_max = $( haskey(kwargs, :E_gam_max) ? kwargs[:E_gam_max] : 30.0e3 )
        
        
        #**** Input for LMB_SPECTRAL_CRs
        # Number of momentum bins for CR model
        CR_nbins = $( haskey(kwargs, :CR_nbins) ? kwargs[:CR_nbins] : 48 )
        CR_pmin = $( haskey(kwargs, :CR_pmin) ? kwargs[:CR_pmin] : 1.0 )
        CR_pmax = $( haskey(kwargs, :CR_pmax) ? kwargs[:CR_pmax] : 1.e6 )
        CR_subsamples = $( haskey(kwargs, :CR_subsamples) ? kwargs[:CR_subsamples] : 10 )
        CR_DSlope = $( haskey(kwargs, :CR_DSlope) ? kwargs[:CR_DSlope] : 1.0e-6 )
        
        
        #**** Set to 1 if you want additional statistical informations (L_x,T,...).
        #         Printed on screen (integer):
        #         (0) disabled
        #         (1) enabled
        GIVE_MORE_INFO = $( haskey(kwargs, :GIVE_MORE_INFO) ? kwargs[:GIVE_MORE_INFO] : 0 )
        
        #=============================================================
        # REAST FRAME - PROJECTION
        #=============================================================
        
        #**** Flag of the image projection (integer), from - to:
        #         (1) along z, xy plane
        #         (2) along y, xz plane
        #         (3) along x, yz plane
        #         (4) along all 3 axis
        PROJECT       = $( haskey(kwargs, :PROJECT) ? kwargs[:PROJECT] : 1 )
        
        #**** Flag of the cluster baricenter motion (integer):
        #         (0) disabled
        #         (1) enabled: at the particles velocity is sobstituted the cluster baricenter velocity
        #         (2) enabled: paricles velocity in the cluster baricenter rest frame
        CENTER_MOTION = $( haskey(kwargs, :CENTER_MOTION) ? kwargs[:CENTER_MOTION] : 0 )
        
        #**** Flag to subtract the velocity of the local group (integer):
        #         (0) disabled
        #         (1) enabled
        REMOVE_LOCAL_GROUP_VEL = $( haskey(kwargs, :REMOVE_LOCAL_GROUP_VEL) ? kwargs[:REMOVE_LOCAL_GROUP_VEL] : 0 )
        
        #**** Flag of the definition of the center (integer):
        #         (0) barycenter of gas particles
        #         (1) selected by user
        #         (2) read by a file
        CENTER = $( haskey(kwargs, :CENTER) ? kwargs[:CENTER] : 1 )
        
        #**** Center position/Cluster data for 'select by user' (reals [Mpc]):
        CENTER_X = $( haskey(kwargs, :CENTER_X) ? kwargs[:CENTER_X] : 0.0 )
        CENTER_Y = $( haskey(kwargs, :CENTER_Y) ? kwargs[:CENTER_Y] : 0.0 )
        CENTER_Z = $( haskey(kwargs, :CENTER_Z) ? kwargs[:CENTER_Z] : 0.0 )
        
        #**** Flag for PERIODIC boxes (integer):
        #         (0) not periodic
        #         (1) periodic
        PERIODIC = $( haskey(kwargs, :PERIODIC) ? kwargs[:PERIODIC] : 0 )
        
        #**** Flag for building lightcones (integer):
        LIGHTCONE = $( haskey(kwargs, :LIGHTCONE) ? kwargs[:LIGHTCONE] : 0 )
        #**** Position of box in lightcone (reals [Mpc]):
        X_ORIGIN = $( haskey(kwargs, :X_ORIGIN) ? kwargs[:X_ORIGIN] : 0.0 )
        Y_ORIGIN = $( haskey(kwargs, :Y_ORIGIN) ? kwargs[:Y_ORIGIN] : 0.0 )
        Z_ORIGIN = $( haskey(kwargs, :Z_ORIGIN) ? kwargs[:Z_ORIGIN] : 0.0 )
        #**** Opening angle of the lightcone (degree):
        OPEN_ANGLE = $( haskey(kwargs, :OPEN_ANGLE) ? kwargs[:OPEN_ANGLE] : 1.0 )
        
        
        #**** File of the centers of the images for 'read by a file':
        MAIN_PROG = $( haskey(kwargs, :MAIN_PROG) ? kwargs[:MAIN_PROG] : "" )
        
        
        #=============================================================
        # OUTPUT FILE
        #=============================================================
        
        #**** Flag string of output file selection:
        #         (1) binary image (IDL)
        #             flag of the subtype of output map:'
        #             (1) header in txt format'
        #             (2) header in FITS format'
        #         (2) ASCII image 
        #         (3) FITS image
        FILE_FORMAT = $( haskey(kwargs, :FILE_FORMAT) ? kwargs[:FILE_FORMAT] : 3 )
        
        #**** Flag for producing HEADER
        #         (0) no header 
        #         (1) header 
        #             for idl /asi format header will be created as ascii file
        FILE_HEADER = $( haskey(kwargs, :FILE_HEADER) ? kwargs[:FILE_HEADER] : 1 )
        
        #**** Prefix of the output files:
        PREFIX_OUT = $( haskey(kwargs, :PREFIX_OUT) ? kwargs[:PREFIX_OUT] : "test" )
        
        #=============================================================
        # COSMOLOGY (for distances computations)
        #=============================================================
        
        #**** Cosmology to use
        HUBBLE = $( haskey(kwargs, :HUBBLE) ? kwargs[:HUBBLE] : (haskey(kwargs, :h) ? kwargs[:h].h0 : 0.72) )
        OMEGA  = $( haskey(kwargs, :OMEGA) ? kwargs[:OMEGA] : (haskey(kwargs, :h) ? kwargs[:h].omega_m : 0.24) )
        LAMBDA = $( haskey(kwargs, :LAMBDA) ? kwargs[:LAMBDA] : (haskey(kwargs, :h) ? kwargs[:h].omega_l : 0.76) )

        """
        )
    end
end
