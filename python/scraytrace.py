"""
scraytrace
==========

.. module:: scraytrace
   :synopsis: Python interface for the solar corona ray-tracing tools.

Python interface for the solar corona ray-tracing tools.

.. todo::
   * Dynamically find the compiled libraries based on the architechture
   * Unit tests
   
   
   
"""


import ctypes
import _ctypes
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


class scraytrace:
    """Solar Corona Raytrace Python interface"""
    pass

    def __init__(self):
        pass


    def raytrace(self):
        """Run the raytracing"""
        pass


# ---- Define tests here




if __name__ == "__main__" :

    libName = "/home/thernis/work/cpp/fromgit/SCRaytrace/build/src/libraytrace.so"
    s = ctypes.cdll.LoadLibrary(libName)
    libHandle = s._handle
    print("Handle : ", libHandle)


    s.dumpbuildinfo()

    result = s.testParamPass(ctypes.c_int(2))

    # -- Attempt to unload the library
    _ctypes.dlclose(libHandle)



    st = ctypes.cdll.LoadLibrary("/home/thernis/work/cpp/fromgit/SCRaytrace/build/src/libraytracethread.so")
    #result = st.rtthreadtestExt()

    st.testPassPython.argtypes = [np.ctypeslib.ndpointer(np.float32, 
                                 flags='aligned, c_contiguous'),
                                 ctypes.c_int,
                                 ctypes.c_float]
                                 
    A = np.array([1.4, 2.6, 3.0], dtype=np.float32)
    st.testPassPython(A, 3, 4.5)


    # ---- Set the parameters for the call to rtthread
    # -- rtthread prototype
    #int sx,                         //!> image size, x axis [in]
    #int sy,                         //!> image size, y axis [in]
    #float fovpix,                   //!> fov angle of one pixel in rad [in]
    #float *obspos,                  //!> [x, y, z] position of the observer in the Sun basis [in]
    #float *obsang,                  //!> [ax, ay, az] orientation of the observer. z is the optical axis. Rotation order is z, y, then x. [in]
    #float *nepos,                   //!> [x, y, z] position of the Ne reference in the Sun basis [in]
    #float *neang,                   //!> [ax, ay, az] orientation of the Ne [in]
    #int losnbp,                     //!> number of steps for the integration along the LOS [in]
    #float *losrange,                //!> [lstart, lend] range for the integration along the LOS in Rsun. The origin of the LOS is the orthogonal projection of the Sun cntr on that LOS. [in]
    #int modelid,                    //!> model id number [in]
    #float *btot,                    //!> Total brightness image (for Thomson scattering physiscs) [out]
    #float *bpol,                    //!> Polarized brightness image (for Thomson scattering physiscs [out]
    #float *netot,                   //!> Electron density image (for Thomson scattering physiscs) [out]
    #float *pmodparam,               //!> parameters of the model [in]
    #float *crpix,                   //!> [x, y] pixel boresight center of the image
    #int quiet,                      //!> quiet mode if set to 1 [in]
    #int neonly,                     //!> Only compute the electron density if set to 1 [in]
    #float *hlonlat,                 //!> [Hlon, Hlat, Hrot] heliographic lon and lat of the center of the disk, rotation angle corresponding to the projection of the north pole, counterclockwise [in].
    #float occrad,                   //!> Occulter radius [Rsun]. The integration in not performed within that disk [in]
    #float limbdark,                 //!> limb darkening coeff: default 0.58 [in]
    #float *obslonlat,               //!> [lon, lat, height] position of the observer in Carrington coordinate. If set, then obspos is ignored. The optical axis always points toward the Sun center. Use obsang to change telescope orientation. Note that obslonlat=[0,0,215] correspond to obspos=[0,0,215] and obsang=[!pi,0,0]: this means that the Carrington coordinate origin on the Solar sphere (lon,lat,height)=(0,0,1) is located at (x,y,z)=(0,0,1), with Ox pointing to solar north and Oy pointing to (lon,lat)=(3*!pi/2,0) [in]
    #int obslonlatflag,              //!> Set to 1 to use obslonlat to position the observer, instead of obspos [in]
    #unsigned int projtypecode,      //!> Projection type. (see Calabretta and Greisen,  Representations of celestial coordinates in FITS, A&A 395, 1077-1122(2002)), ARC : Zenithal equidistant (default), TAN : Gnomonic, SIN : Slant orthographic, AZP : Zenithal perspective. [in]
    #float pv2_1,                    //!> mu parameter for the AZP projection [in]
    #float *pc,                      //!> wcs pc[2, 2] matrix: default is unit matrix [in]
    #int frontinteg,                 //!> Set to 1 so that the origin of the LOS is taken at the observer: if set, the losrange parameters must both be positive. [in]
    #unsigned int nbthreads,         //!> Number of threads to run in parallel. [in]
    #unsigned int nbchunk,           //!> Number of chuncks. Use with nbthread. If set to a value less than 2, the threads are  launched by lines of sight. If nbchunks >= 2, the threads are launched by chunk of the image. Ballance nbthreads and nbchunks for optimal performances. [in]
    #float *nerotcntr,               //!> [x, y, z] center of rotation of the Ne model, in the Ne basis [in]
    #float *nerotang,                //!> [ax, ay, az] rotation of the Ne model around the nerotcntr, in the Ne basis [in]
    #float *netranslation,           //!> [tx, ty, tz] translation vector of the Ne model, in the Ne basis [in]
    #int *nerotaxis,                 //!> [axid1, axid2, axid3] axis id corresponding to the nerotang rotation angles. 1: X, 2: Y, 3: Z. Default is [3,2,1]. [in]
    #int physics,                    //!> type of physics to perform the raytracing [in]
    #float *phyparam,                //!> Extra parameters required depending on the chosen physics [in]
    #float fracmax,                  //!> Set it to the fraction of the maximum total B per LOS in order to compute the distance to that fraction of brightness. Disabled if set to 0 [default]. The distance is returned in the bpol image, in Rsun. [in]
    #int runDumpInteg,               //!> Set if you want to save all the integration points in the integrand variable. Note that this feature can require the allocation of a large amount of memory, typically a floating array of imsize[0] x imsize[1] x losnbp. Use this feature only if you have enough free memory. [in]
    #float *pIntegrand);             //!> Contains the all the integration points if rundumpinteg is set [out]


    st.rtthread.argtypes = [ctypes.c_int,          # sx
                            ctypes.c_int,           # sy
                            ctypes.c_float,       # fovpix
                            np.ctypeslib.ndpointer(np.float32,    # obspos 
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # obsang 
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # nepos 
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # neang 
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # losnbp
                            np.ctypeslib.ndpointer(np.float32,    # losrange 
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # modelid
                            np.ctypeslib.ndpointer(np.float32,    # btot
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # bpol
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # netot
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # pmodparam
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # crpix
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # quiet
                            ctypes.c_int,           # neonly
                            np.ctypeslib.ndpointer(np.float32,    # hlonlat
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_float,         # occrad
                            ctypes.c_float,         # limbdark
                            np.ctypeslib.ndpointer(np.float32,    # obslonlat
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # obslonlatflag
                            ctypes.c_uint,          # projtypecode
                            ctypes.c_float,         # pv2_1
                            np.ctypeslib.ndpointer(np.float32,    # pc
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # frontinteg
                            ctypes.c_uint,          # nbthreads
                            ctypes.c_uint,          # nbchunks
                            np.ctypeslib.ndpointer(np.float32,    # nerotcntr
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # nerotang
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.float32,    # netranslation
                                                   flags='aligned, c_contiguous'),
                            np.ctypeslib.ndpointer(np.int32,    # nerotaxis
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_int,           # physics
                            np.ctypeslib.ndpointer(np.float32,    # phyparam
                                                   flags='aligned, c_contiguous'),
                            ctypes.c_float,         # fracmax
                            ctypes.c_int,           #runDumpInteg
                            np.ctypeslib.ndpointer(np.float32,    # pIntegrand
                                                   flags='aligned, c_contiguous')]


    # ---- Set all variables before the call

    profsize = 250
    fovedgesundeg = 0.5
    ffovdeg = 160.
    beta = 0.
    DistToSun_AU = 0.16
    DisttoSun_Rsun = 215. * DistToSun_AU


    #int sx,                         //!> image size, x axis [in]
    sx = profsize

    #int sy,                         //!> image size, y axis [in]
    sy = profsize

    #float fovpix,                   //!> fov angle of one pixel in rad [in]
    fovpix = np.deg2rad(ffovdeg / profsize) 

    #float *obspos,                  //!> [x, y, z] position of the observer in the Sun basis [in]
    obspos = np.array([0., 0., -DisttoSun_Rsun], dtype=np.float32)

    #float *obsang,                  //!> [ax, ay, az] orientation of the observer. z is the optical axis. Rotation order is z, y, then x. [in]
    obsang = np.deg2rad(np.array([-fovedgesundeg - ffovdeg / 2., 0, 0], dtype=np.float32))

    #float *nepos,                   //!> [x, y, z] position of the Ne reference in the Sun basis [in]
    nepos = np.array([0, 0, 0], dtype=np.float32)

    #float *neang,                   //!> [ax, ay, az] orientation of the Ne [in]
    neang = np.deg2rad(np.array([0., 0., beta], dtype=np.float32))

    #int losnbp,                     //!> number of steps for the integration along the LOS [in]
    losnbp = ctypes.c_int(100)

    #float *losrange,                //!> [lstart, lend] range for the integration along the LOS in Rsun. The origin of the LOS is the orthogonal projection of the Sun cntr on that LOS. [in]
    losrange = np.array([0. , 215. * 3], dtype=np.float32)

    #int modelid,                    //!> model id number [in]
    modelid = ctypes.c_int(78)

    #float *btot,                    //!> Total brightness image (for Thomson scattering physiscs) [out]
    btot = np.zeros((sy, sx), dtype=np.float32)

    #float *bpol,                    //!> Polarized brightness image (for Thomson scattering physiscs [out]
    bpol = np.zeros((sy, sx), dtype=np.float32)

    #float *netot,                   //!> Electron density image (for Thomson scattering physiscs) [out]
    netot = np.zeros((sy, sx), dtype=np.float32)

    #float *pmodparam,               //!> parameters of the model [in]
    modparam = np.array([690. / 1361., 1., 1.], dtype=np.float32)

    #float *crpix,                   //!> [x, y] pixel boresight center of the image
    crpix = (np.array([sx, sy], dtype=np.float32) / 2) - 0.5

    #int quiet,                      //!> quiet mode if set to 1 [in]
    quiet = ctypes.c_int(0)

    #int neonly,                     //!> Only compute the electron density if set to 1 [in]
    neonly = ctypes.c_int(0)

    #float *hlonlat,                 //!> [Hlon, Hlat, Hrot] heliographic lon and lat of the center of the disk, rotation angle corresponding to the projection of the north pole, counterclockwise [in].
    hlonlat = np.array([0, 0, 0], dtype=np.float32)

    #float occrad,                   //!> Occulter radius [Rsun]. The integration in not performed within that disk [in]
    occrad = ctypes.c_float(0.)

    #float limbdark,                 //!> limb darkening coeff: default 0.58 [in]
    limbdark = ctypes.c_float(0.58)

    #float *obslonlat,               //!> [lon, lat, height] position of the observer in Carrington coordinate. If set, then obspos is ignored. The optical axis always points toward the Sun center. Use obsang to change telescope orientation. Note that obslonlat=[0,0,215] correspond to obspos=[0,0,215] and obsang=[!pi,0,0]: this means that the Carrington coordinate origin on the Solar sphere (lon,lat,height)=(0,0,1) is located at (x,y,z)=(0,0,1), with Ox pointing to solar north and Oy pointing to (lon,lat)=(3*!pi/2,0) [in]
    obslonlat = np.array([0, 0, 0], dtype=np.float32)

    #int obslonlatflag,              //!> Set to 1 to use obslonlat to position the observer, instead of obspos [in]
    obslonlatflag = ctypes.c_int(0)

    #unsigned int projtypecode,      //!> Projection type. (see Calabretta and Greisen,  Representations of celestial coordinates in FITS, A&A 395, 1077-1122(2002)), ARC : Zenithal equidistant (default), TAN : Gnomonic, SIN : Slant orthographic, AZP : Zenithal perspective. [in]
    projtypecode = ctypes.c_uint(1)

    #float pv2_1,                    //!> mu parameter for the AZP projection [in]
    pv2_1 = ctypes.c_float(0)

    #float *pc,                      //!> wcs pc[2, 2] matrix: default is unit matrix [in]
    pc = np.eye(2, dtype=np.float32)

    #int frontinteg,                 //!> Set to 1 so that the origin of the LOS is taken at the observer: if set, the losrange parameters must both be positive. [in]
    frontinteg = ctypes.c_int(1)

    #unsigned int nbthreads,         //!> Number of threads to run in parallel. [in]
    nbthreads = ctypes.c_uint(16)

    #unsigned int nbchunk,           //!> Number of chuncks. Use with nbthread. If set to a value less than 2, the threads are  launched by lines of sight. If nbchunks >= 2, the threads are launched by chunk of the image. Ballance nbthreads and nbchunks for optimal performances. [in]
    nbchunks = ctypes.c_uint(16)

    #float *nerotcntr,               //!> [x, y, z] center of rotation of the Ne model, in the Ne basis [in]
    nerotcntr = np.array([0, 0, 0], dtype=np.float32)

    #float *nerotang,                //!> [ax, ay, az] rotation of the Ne model around the nerotcntr, in the Ne basis [in]
    nerotang = np.array([0, 0, 0], dtype=np.float32)

    #float *netranslation,           //!> [tx, ty, tz] translation vector of the Ne model, in the Ne basis [in]
    netranslation = np.array([0, 0, 0], dtype=np.float32)

    #int *nerotaxis,                 //!> [axid1, axid2, axid3] axis id corresponding to the nerotang rotation angles. 1: X, 2: Y, 3: Z. Default is [3,2,1]. [in]
    nerotaxis = np.array([3, 2, 1], dtype=np.int32)

    #int physics,                    //!> type of physics to perform the raytracing [in]
    physics = ctypes.c_int(5)

    #float *phyparam,                //!> Extra parameters required depending on the chosen physics [in]
    phyparam = np.array([0], dtype=np.float32)

    #float fracmax,                  //!> Set it to the fraction of the maximum total B per LOS in order to compute the distance to that fraction of brightness. Disabled if set to 0 [default]. The distance is returned in the bpol image, in Rsun. [in]
    fracmax = ctypes.c_float(0.)

    #int runDumpInteg,               //!> Set if you want to save all the integration points in the integrand variable. Note that this feature can require the allocation of a large amount of memory, typically a floating array of imsize[0] x imsize[1] x losnbp. Use this feature only if you have enough free memory. [in]
    runDumpInteg = ctypes.c_int(0)

    #float *pIntegrand);             //!> Contains the all the integration points if rundumpinteg is set [out]
    integrand = np.array([0], dtype=np.float32)


    # ---- Call rtthread
    st.rtthread(ctypes.c_int(sx),
                ctypes.c_int(sy),
                fovpix,
                obspos,
                obsang,
                nepos,
                neang,
                losnbp,
                losrange,
                modelid,
                btot,
                bpol,
                netot,
                modparam,
                crpix,
                quiet,
                neonly,
                hlonlat,
                occrad,
                limbdark,
                obslonlat,
                obslonlatflag,
                projtypecode,
                pv2_1,
                pc,
                frontinteg,
                nbthreads,
                nbchunks,
                nerotcntr,
                nerotang,
                netranslation,
                nerotaxis,
                physics,
                phyparam,
                fracmax,
                runDumpInteg,
                integrand)


    bsun2s10 = 1. / 4.33E-16

    elongSim = np.linspace(fovedgesundeg, fovedgesundeg + ffovdeg, profsize)

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(1, 1, 1)

    #minmax = [1e-12, 1e-7]
    minmax = [1e-15, 1e-10]

    norm = mpl.colors.LogNorm(vmin=minmax[0], vmax=minmax[1])
    ax1.imshow(btot, origin='lower', cmap='gist_ncar', norm=norm)

    #ax1.loglog(elongSim, btot[0,:] * bsun2s10)

    #ax1.set_xlabel('Elongation from Sun Center [deg]')
    #ax1.set_ylabel('Brightness [S10]')
    #ax1.grid()

    plt.show()


