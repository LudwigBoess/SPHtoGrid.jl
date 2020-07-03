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
    boxsize_pix::Float32        # xy-size of the image in pixels
    pixsize_kpc::Float32        # size of one pixel in kpc
    xlim::Array{Float64,1}      # x limits of image
    ylim::Array{Float64,1}      # y limits of image
    zlim::Array{Float64,1}      # z limits of image
    units::String               # unitstring of image

    function Smac1ImageInfo(snap::Int32, z::Float32, m_vir::Float32, r_vir::Float32,
                       xcm::Float32, ycm::Float32, zcm::Float32,
                       z_slice_kpc::Float32,
                       boxsize_kpc::Float32, boxsize_pix::Float32, pixsize_kpc::Float32,
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
    boxsize_pix = read(f, Float32)
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
    write_smac1_par([...])

Writes a Smac parameter file. Not all relevant parameters implemented yet!
"""
function write_smac1_par(path, use_keys, out_dir, snap_base, snap_file, image_prefix,
                         halo_id, snap_start, snap_end, file_format, kernel,
                         xy_size, z_size, out_effect, out_subeffect,
                         cr_bins, cr_pmin, cr_pmax, cr_subsamples,
                         projection, x0, y0, z0, Npixels=1024)

    open(path * "smac1.par", "w") do f
        write(f,
        """#=============================================================
        # Input file
        #=============================================================

        SIM_FORMAT = 4

        #**** Input file
        USE_KEYS = $use_keys
        OUTPUT_DIR = $out_dir
        SNAP_BASE = $snap_base

        SNAP_FILE = $snap_file

        #**** Prefix of the output files:
        PREFIX_OUT = $image_prefix


        #**** Halo ID
        HALO_ID = $halo_id

        SNAP_START = $snap_start
        SNAP_END   = $snap_end

        INITIAL_Z = 0.0
        IMAGE_Z   = 0.0

        FILE_FORMAT = $file_format

        FILE_HEADER = 1

        KERNEL_TYPE = $kernel

        PART_DISTR = 1

        IMG_XY_UNITS = 2
        IMG_XY_SIZE  = $xy_size

        IMG_Z_UNITS = 2
        IMG_Z_SIZE  = $z_size
        FILL_CARROT = 0.

        IMG_SIZE = $Npixels

        SMOOTH_UNITS = 0
        SMOOTH_FWHM  = 5.

        NSIDE = 8

        MIN_DIST =   5000.
        MAX_DIST = 420000.

        OUTPUT_MAP = $out_effect
        OUTPUT_SUB = $out_subeffect

        METAL_SPECIES = 0

        TEMP_UNIT = 1

        TEMP_CUT = 0.5

        XRAY_E0 = 0.5
        XRAY_E1 = 2.

        X_TABLES =

        FRQ_ZS = 300.

        #**** Input for Radio Maps
        #XCRP_TABLE_PATH gives the change of XCRP over radius.
        ECRP_MIN = 1e9
        FRQ_Pnu = 1.4e9
        XCRP = 0.01
        XCRP_TABLE_PATH = ~/
        B_TABLE_PATH = ~/
        KERNEL_TABLE_PATH =
        GAM_nu = 1.25
        #XCRP_TABLE_FILE_FMT = (1F8.6,2X,1F14.6)

        #energies in GeV
        E_gam = 1e1
        E_gam_min = 1e-2
        E_gam_max = 30e3

        # Number of momentum bins for CR model
        CR_nbins = $cr_bins
        CR_pmin = $cr_pmin
        CR_pmax = $cr_pmax
        CR_subsamples = $cr_subsamples
        CR_DSlope = 1.0e-6

        GIVE_MORE_INFO = 0

        PROJECT       = $projection

        CENTER_MOTION = 0

        REMOVE_LOCAL_GROUP_VEL = 0

        CENTER = 1

        CENTER_X = $x0
        CENTER_Y = $y0
        CENTER_Z = $z0

        PERIODIC = 1

        LIGHTCONE = 0
        X_ORIGIN = 0.
        Y_ORIGIN = 0.
        Z_ORIGIN = 0.
        OPEN_ANGLE = 1.0

        MAIN_PROG = /afs/mpa/project/hydrosims/Hutt/g676/csf/Post/main_prog.a_gas.gv

        HUBBLE = 0.7
        OMEGA  = 0.3
        LAMBDA = 0.7

        """
        )
    end
end
