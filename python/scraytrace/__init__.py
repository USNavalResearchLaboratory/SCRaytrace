"""
scraytrace
==========

.. module:: scraytrace
   :synopsis: Python interface for the solar corona ray-tracing tools.

Python interface for the solar corona ray-tracing tools.

.. todo::
   * Dynamically find the compiled libraries based on the architechture
   * output the C++cout to python interpreter
   
"""

import ctypes
import numpy as np
from astropy import units as u
from astropy import wcs
from astropy.constants import iau2015
import yaml
import pathlib

import mathutil as mu


__all__ = ['rtrotmat2rxryrz', 'model54_calcLegHeight', 
           'model54_calcLeadingEdgeHeight', 'scraytrace']

__doc__ = "Python interface for the solar corona ray-tracing tools."


def rtrotmat2rxryrz(m):
    """Get the rx,ry,rz angles from a rotation matrix 

; CATEGORY:
;  raytracing, mathematics, 3d, geometry
;
; INPUTS:
;  m : 3x3 rotation matrix
;
; OUTPUTS:
;  [rx,ry,rz] in rad
;
; DESCRIPTION:
;  - See function angle3123 in dragger.pro
;
;  - The rotation matrix is assumed to be computed in that order:
;    R(-lon,x).R(lat,y).R(-roll,z)
;    with R(angle,axis) the rotation matrix around axis "axis" with an
;    angle "angle".
;  - Coordinate system orientation is assumed to follow 
;    raytracewl software package convention, that, I know, is not
;    usual:
;    X : vertical axis
;    Y : horizontal axis
;    Z : depth axis, perpendicular to the plane of the sky
;    
    
    
    """

    ry = np.arcsin(m[2, 0])
    c1 = np.cos(ry)
    if np.abs(c1) < 1.0e-6:
        rz = -np.arctan2(m[1, 2], m[0, 2])
        rx = 0.0
    else:
        rz = -np.arctan2(m[1, 0], m[0, 0])
        rx = -np.arctan2(m[2, 1], m[2, 2])
    
    return np.array([rx, ry, rz])


def model54_calcLegHeight(leadingEdgeHeight, k, alpha):
    """Compute h, the leg height parameter of GCS model 54
    
    See Thernisien, A., "IMPLEMENTATION OF THE GRADUATED CYLINDRICAL SHELL MODEL FOR THETHREE-DIMENSIONAL RECONSTRUCTION OF CORONAL MASS EJECTIONS", ApJ Supplement Series, 194:33  (6pp), 2011 June
    
    args:
        leadingEdgeHeight: height of the CME leading edge
        k: kappa parameter
        alpha: alpha angle in the publication, half croissant half angle
        
    return:
        h
    
    """    
    legHeight = leadingEdgeHeight * (1 - k) * np.cos(alpha) / (1 + np.sin(alpha))
    return legHeight
    

def model54_calcLeadingEdgeHeight(legHeight, k, alpha):
    """Compute the leading edge height of the GCS model 54
    
    See Thernisien, A., "IMPLEMENTATION OF THE GRADUATED CYLINDRICAL SHELL MODEL FOR THETHREE-DIMENSIONAL RECONSTRUCTION OF CORONAL MASS EJECTIONS", ApJ Supplement Series, 194:33  (6pp), 2011 June
    
    args:
        legHeight: parameter h of the model
        k: kappa parameter
        alpha: alpha angle in the publication, half croissant half angle
        
    return:
        eadingEdgeHeight
    
    """
    eadingEdgeHeight = legHeight / ((1 - k) * np.cos(alpha) / (1 + np.sin(alpha)))
    return eadingEdgeHeight


def obsang2crval(obsang):
    """Compute the CRVAL from the obsang"""
    



class scraytraceLib:
    """Access point to the C++ library"""

    # -- read configuration file
    #    Edit that file to point to the compiled scraytrace library
    thisFilesPath = pathlib.Path(__file__).parents[0]
    fnConfig = thisFilesPath.joinpath('RTLocalEnvPath.yaml')
    
    # -- check if local config file exist
    if not fnConfig.exists():
        # !!! DO NOT EDIT THE LINES BELOW. 
        #     EDIT THE RTLocalEnvPath.yaml file instead
        localConfTemplate = """
rt_libpath: /Users/BrandonBonifacio/build/src/
rt_libfile: libraytracethread.dylib
        
        """
        with open(fnConfig, 'w') as file:
            file.write(localConfTemplate)
    
    with open(fnConfig) as f:    
        localConf = yaml.load(f, Loader=yaml.FullLoader)

    soFullFileName = pathlib.Path(localConf['rt_libpath']).joinpath(localConf['rt_libfile'])

    try:
        st = ctypes.cdll.LoadLibrary(soFullFileName)
    except :
        raise Exception('Cannot find the library. Please edit the RTLocalEnvPath.yaml configuration file.')
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




class scraytrace:
    
    projtypecodehash = {'ARC':1,
                        'TAN':2,
                        'SIN':3,
                        'AZP':4}

        
    projtypenamehash = {1:'ARC',
                        2:'TAN',
                        3:'SIN',
                        4:'AZP'}
    
    
    
    """Solar Corona Raytrace Python interface
    
    
    # ---- Set the parameters for the call to rtthread
    # -- rtthread prototype
    #int sx,                         //!> image size, x axis [in]
    #int sy,                         //!> image size, y axis [in]
    #float fovpix,                   //!> fov angle of one pixel in rad [in]
    #float *obspos,                  //!> [x, y, z] position of the observer in the Sun basis [in] (Units: Rsun)
    #float *obsang,                  //!> [ax, ay, az] orientation of the observer. z is the optical axis. Rotation order is z, y, then x. [in] (Units: Radians)
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
    
    """
    
    lib = scraytraceLib()

    def __init__(self,
                 imsize=(250, 250),
                 fovpix=np.deg2rad(10.)/250,
                 obspos=[0., 0., -215.],
                 obsang=[0., 0., 0.],
                 nepos=[0., 0., 0.],
                 neang=[0., 0., 0.],
                 losnbp=100,
                 losrange=[0, 215. * 2],
                 modelid=78,
                 modparam = [690. / 1361., 1.],
                 crpix=None,
                 quiet=False,
                 neonly=False,
                 hlonlat=[0., 0., 0.],
                 occrad=0.,
                 obslonlat=None,
                 obslonlatflag=False,
                 projtypecode=1,
                 pv2_1=0.,
                 pc=np.eye(2, dtype=np.float32),
                 frontinteg=True,
                 nbthreads=16,
                 nbchunks=16,
                 nerotcntr=[0., 0., 0.],
                 nerotang=[0., 0., 0.],
                 netranslation=[0., 0., 0.],
                 nerotaxis=[3, 2, 1],
                 physics=5,
                 phyparam=[0],
                 fracmax=0.,
                 runDumpInteg=False):
        

        if obslonlat is not None:
            obspos = obslonlat[2] * np.array([np.sin(obslonlat[1]),
                                     np.sin(obslonlat[0]) * np.cos(obslonlat[1]),
                                     -np.cos(obslonlat[0]) * np.cos(obslonlat[1])])
            obslonlatflag = 1
        else:
            obslonlat = self._computeObslonlat(obspos)
            obslonlatflag = 0



        self.imsize = np.array(imsize)
        
        #int sx,                         //!> image size, x axis [in]
        self.sx = imsize[0]

        #int sy,                         //!> image size, y axis [in]
        self.sy = imsize[1]

        #float fovpix,                   //!> fov angle of one pixel in rad [in]
        self.fovpix = fovpix

        #float *obspos,                  //!> [x, y, z] position of the observer in the Sun basis [in]
        self.obspos = np.array(obspos, dtype=np.float32)

        #float *obsang,                  //!> [ax, ay, az] orientation of the observer. z is the optical axis. Rotation order is z, y, then x. [in]
        self.obsang = np.array(obsang, dtype=np.float32)

        #float *nepos,                   //!> [x, y, z] position of the Ne reference in the Sun basis [in]
        self.nepos = np.array(nepos, dtype=np.float32)

        #float *neang,                   //!> [ax, ay, az] orientation of the Ne [in]
        self.neang = np.array(neang, dtype=np.float32)

        #int losnbp,                     //!> number of steps for the integration along the LOS [in]
        self.losnbp = losnbp

        #float *losrange,                //!> [lstart, lend] range for the integration along the LOS in Rsun. The origin of the LOS is the orthogonal projection of the Sun cntr on that LOS. [in]
        self.losrange = np.array(losrange, dtype=np.float32)

        #int modelid,                    //!> model id number [in]
        self.modelid = modelid

        #float *btot,                    //!> Total brightness image (for Thomson scattering physiscs) [out]
        self.btot = np.zeros((self.sy, self.sx), dtype=np.float32)

        #float *bpol,                    //!> Polarized brightness image (for Thomson scattering physiscs [out]
        self.bpol = np.zeros((self.sy, self.sx), dtype=np.float32)

        #float *netot,                   //!> Electron density image (for Thomson scattering physiscs) [out]
        self.netot = np.zeros((self.sy, self.sx), dtype=np.float32)

        #float *pmodparam,               //!> parameters of the model [in]
        self.modparam = np.array(modparam, dtype=np.float32)

        #float *crpix,                   //!> [x, y] pixel boresight center of the image
        if crpix is None:
            crpix = self.imsize / 2 - 0.5
        self.crpix = np.array(crpix, dtype=np.float32)

        #int quiet,                      //!> quiet mode if set to 1 [in]
        self.quiet = quiet

        #int neonly,                     //!> Only compute the electron density if set to 1 [in]
        self.neonly = neonly

        #float *hlonlat,                 //!> [Hlon, Hlat, Hrot] heliographic lon and lat of the center of the disk, rotation angle corresponding to the projection of the north pole, counterclockwise [in].
        self.hlonlat = np.array(hlonlat, dtype=np.float32)

        #float occrad,                   //!> Occulter radius [Rsun]. The integration in not performed within that disk [in]
        self.occrad = occrad

        #float limbdark,                 //!> limb darkening coeff: default 0.58 [in]
        # limbdark should now be passed in phyparam[0]
        self.limbdark = 0.

        #float *obslonlat,               //!> [lon, lat, height] position of the observer in Carrington coordinate. If set, then obspos is ignored. The optical axis always points toward the Sun center. Use obsang to change telescope orientation. Note that obslonlat=[0,0,215] correspond to obspos=[0,0,215] and obsang=[!pi,0,0]: this means that the Carrington coordinate origin on the Solar sphere (lon,lat,height)=(0,0,1) is located at (x,y,z)=(0,0,1), with Ox pointing to solar north and Oy pointing to (lon,lat)=(3*!pi/2,0) [in]
        self.obslonlat = np.array(obslonlat, dtype=np.float32)

        #int obslonlatflag,              //!> Set to 1 to use obslonlat to position the observer, instead of obspos [in]
        self.obslonlatflag = obslonlatflag

        #unsigned int projtypecode,      //!> Projection type. (see Calabretta and Greisen,  Representations of celestial coordinates in FITS, A&A 395, 1077-1122(2002)), ARC : Zenithal equidistant (default), TAN : Gnomonic, SIN : Slant orthographic, AZP : Zenithal perspective. [in]
        self.projtypecode = projtypecode

        #float pv2_1,                    //!> mu parameter for the AZP projection [in]
        self.pv2_1 = pv2_1

        #float *pc,                      //!> wcs pc[2, 2] matrix: default is unit matrix [in]
        self.pc = np.array(pc, dtype=np.float32)

        #int frontinteg,                 //!> Set to 1 so that the origin of the LOS is taken at the observer: if set, the losrange parameters must both be positive. [in]
        self.frontinteg = frontinteg

        #unsigned int nbthreads,         //!> Number of threads to run in parallel. [in]
        self.nbthreads = nbthreads

        #unsigned int nbchunk,           //!> Number of chuncks. Use with nbthread. If set to a value less than 2, the threads are  launched by lines of sight. If nbchunks >= 2, the threads are launched by chunk of the image. Ballance nbthreads and nbchunks for optimal performances. [in]
        self.nbchunks = nbchunks

        #float *nerotcntr,               //!> [x, y, z] center of rotation of the Ne model, in the Ne basis [in]
        self.nerotcntr = np.array(nerotcntr, dtype=np.float32)

        #float *nerotang,                //!> [ax, ay, az] rotation of the Ne model around the nerotcntr, in the Ne basis [in]
        self.nerotang = np.array(nerotang, dtype=np.float32)

        #float *netranslation,           //!> [tx, ty, tz] translation vector of the Ne model, in the Ne basis [in]
        self.netranslation = np.array(netranslation, dtype=np.float32)

        #int *nerotaxis,                 //!> [axid1, axid2, axid3] axis id corresponding to the nerotang rotation angles. 1: X, 2: Y, 3: Z. Default is [3,2,1]. [in]
        self.nerotaxis = np.array(nerotaxis, dtype=np.int32)

        #int physics,                    //!> type of physics to perform the raytracing [in]
        self.physics = physics        # 5: vsf vary dist
        #physics = ctypes.c_int(0)       # 0: Thomson scattering

        #float *phyparam,                //!> Extra parameters required depending on the chosen physics [in]
        self.phyparam = np.array(phyparam, dtype=np.float32)

        #float fracmax,                  //!> Set it to the fraction of the maximum total B per LOS in order to compute the distance to that fraction of brightness. Disabled if set to 0 [default]. The distance is returned in the bpol image, in Rsun. [in]
        self.fracmax = fracmax

        #int runDumpInteg,               //!> Set if you want to save all the integration points in the integrand variable. Note that this feature can require the allocation of a large amount of memory, typically a floating array of imsize[0] x imsize[1] x losnbp. Use this feature only if you have enough free memory. [in]
        self.runDumpInteg = runDumpInteg

        #float *pIntegrand);             //!> Contains the all the integration points if rundumpinteg is set [out]
        self.integrand = np.array([0], dtype=np.float32)

        self.crval = np.array([-self.obsang[0], self.obsang[1]])

        self.ctype = ['HPLN-' + self.projtypenamehash[self.projtypecode],
                               'HPLT-' + self.projtypenamehash[self.projtypecode]]



    def _computeObslonlat(self, obspos):
        """Compute the obslonlat parameter from obspos"""
        obslonlat = [np.arctan2(obspos[1], -obspos[2]),
                     np.arcsin(obspos[0] / np.linalg.norm(obspos)),
                     np.linalg.norm(obspos)]
        return obslonlat



    def setParamFromWCS(self, fitsHeader, imsize, rollang=0.,
                        frame='CAR'):
        """Set the parameters using the WCS header  
            Args:
                fitsHeader: Fits header with WCS compliant keywords
                imsize (int, int): Size of the image you want to raytrace, in pixels
                frame (str): 'CAR' for carrington
                             'HAE' for Heliocentric Ares Ecliptic
                
        """
        
        self.sx = imsize[0]
        self.sy = imsize[1]        
        
        w = wcs.WCS(fitsHeader)        

        oneRsun_m = iau2015.R_sun.to_value(u.meter)
        
        if frame.upper() == 'CAR':
            self.obslonlat = np.array([np.deg2rad(fitsHeader['CRLN_OBS']),
                                       np.deg2rad(fitsHeader['CRLT_OBS']),
                                       fitsHeader['DSUN_OBS'] / oneRsun_m], dtype=np.float32)
                              
            self.obslonlatflag = True
        elif frame.upper() == 'HAE':
            # x_RT =  Z_HAE
            # y_RT = -Y_HAE
            # z_RT =  X_HAE
            # From Franz & Harper 2002:
            # XY plane: Earth mean ecliptic at J2000.0
            # +X: First point of Aries: Earth-Sun vector of vernal equinox at J2000.0
            self.obspos = np.array([ fitsHeader['HAEZ_OBS'],
                                    -fitsHeader['HAEY_OBS'],
                                     fitsHeader['HAEX_OBS']]) / oneRsun_m
                                    
            self.obslonlat = self._computeObslonlat(self.obspos)
            self.obslonlatflag = True

        else:
            raise Exception("Frame code {0} not supported".format(frame))
        
        
        imszratio = fitsHeader['NAXIS1'] / self.sx
        
        cdelt1 =  u.Quantity(fitsHeader['CDELT1'], fitsHeader['CUNIT1'])
        self.fovpix = cdelt1.to_value(u.rad) * imszratio
              
        # ---- get the angular coordinates of the center of the image: optical axis
        crval1 = u.Quantity(fitsHeader['CRVAL1'], fitsHeader['CUNIT1']).to_value(u.rad)
        crval2 = u.Quantity(fitsHeader['CRVAL2'], fitsHeader['CUNIT2']).to_value(u.rad)
        self.crval = np.array([crval1, crval2])
        
        self.crpix = mu.piximchangereso(w.wcs.crpix - np.array([1, 1]), 
                                        -np.log10(imszratio) / np.log10(2))
    
    
        # ---- compute the spacecraft attitude
        rmat = np.dot(mu.rotmat(-self.crval[1], 2, retArray=True), mu.rotmat(self.crval[0], 1, retArray=True))  
        rmat = np.dot(rmat, mu.rotmat(rollang, 3, retArray=True))
        self.obsang = rtrotmat2rxryrz(rmat)

        self.pc = w.wcs.pc
        
        # ---- get projection type
        projtypename = w.wcs.ctype[0][-3:]
        
        try:
            self.projtypecode = self.projtypecodehash[projtypename]
        except:
            print("{0} projection not supported yet. Will use the ARC instead.".format(projtypename))
            self.projtypecode = 1
            projtypename = 'ARC'
            
        
        if projtypename == 'AZP':
            self.pv2_1 = fitsHeader['PV2_1']
        else:
            self.pv2_1 = 0
        
        self.ctype = [w.wcs.ctype[0][:5] + projtypename, w.wcs.ctype[1][:5] + projtypename]        
        


    def printRTParam(self):
        """Print the parameters of the simulation"""
        print()
        print('Raytracing parameter summary')
        print('sx : ', self.sx)
        print('sy : ', self.sy)
        print('fovpix : ', self.fovpix)
        print('obspos : ', self.obspos)
        print('obsang : ', self.obsang)
        print('nepos : ', self.nepos)
        print('neang : ', self.neang)
        print('losnbp : ', self.losnbp)
        print('losrange : ', self.losrange)
        print('modelid : ', self.modelid)
        print('modparam : ', self.modparam)
        print('crpix : ', self.crpix)
        print('quiet : ', self.quiet)
        print('neonly : ', self.neonly)
        print('hlonlat : ', self.hlonlat)
        print('occrad : ', self.occrad)
        print('obslonlat : ', self.obslonlat)
        print('obslonlatflag')
        print('projtypecode : ', self.projtypecode)        
        print('pv2_1 : ', self.pv2_1)        
        print('pc : ', self.pc)        
        print('frontinteg : ', self.frontinteg)        
        print('nbthreads : ', self.nbthreads)        
        print('nbchunks : ', self.nbchunks)        
        print('nerotcntr : ', self.nerotcntr)        
        print('nerotang : ', self.nerotang)        
        print('netranslation : ', self.netranslation)        
        print('nerotaxis : ', self.nerotaxis)        
        print('physics : ', self.physics)        
        print('phyparam : ', self.phyparam)        
        print('fracmax : ', self.fracmax)        
        print('runDumpInteg : ', self.runDumpInteg)        

        
        
        
        
    def createWCSHead(self):
        """Create the WCS compliant header based on raytracing parameters"""        
        
        self.wcsrt = wcs.WCS(naxis=2)
        self.wcsrt.wcs.crpix = self.crpix + np.array([1, 1])
        self.wcsrt.wcs.crval = np.rad2deg(self.crval)
        self.wcsrt.wcs.ctype = self.ctype
        self.wcsrt.wcs.pc = self.pc
        self.wcsrt.wcs.cdelt = np.rad2deg([self.fovpix, self.fovpix])
        self.wcsrt.wcs.cunit = ['deg', 'deg']

        self.fitsrt = self.wcsrt.to_header()
        self.fitsrt['NAXIS1'] = self.sx
        self.fitsrt['NAXIS2'] = self.sy
        self.fitsrt['PV2_1'] = self.pv2_1
        


    def raytrace(self):
        """Run the raytracing"""
        
        self.btot = np.zeros((self.sy, self.sx), dtype=np.float32)
        self.bpol = np.zeros((self.sy, self.sx), dtype=np.float32)
        self.netot = np.zeros((self.sy, self.sx), dtype=np.float32)
        self.integrand = np.array([0], dtype=np.float32)


        # -- Trick to make things work:
        #    Transpose the pc matrix before C++ call. Using the built-in numpy transpose makes it crash. Don't know why yet...
        #    The pc matrix is transposed back after the C++ call. This means that there might be a problem in the C code...
#        tmp = self.pc[0, 1]
#        self.pc[0, 1] = self.pc[1, 0]
#        self.pc[1, 0] = tmp
        


        # ---- Call rtthread
        self.lib.st.rtthread(ctypes.c_int(self.sx),
                    ctypes.c_int(self.sy),
                    self.fovpix,
                    np.array(self.obspos, dtype=np.float32),
                    np.array(self.obsang, dtype=np.float32),
                    np.array(self.nepos, dtype=np.float32),
                    np.array(self.neang, dtype=np.float32),
                    ctypes.c_int(self.losnbp),
                    np.array(self.losrange, dtype=np.float32),
                    ctypes.c_int(self.modelid),
                    self.btot,
                    self.bpol,
                    self.netot,
                    np.array(self.modparam, dtype=np.float32),
                    np.array(self.crpix, dtype=np.float32),
                    ctypes.c_int(self.quiet),
                    ctypes.c_int(self.neonly),
                    np.array(self.hlonlat, dtype=np.float32),
                    ctypes.c_float(self.occrad),
                    ctypes.c_float(self.limbdark),
                    np.array(self.obslonlat, dtype=np.float32),
                    ctypes.c_int(self.obslonlatflag),
                    ctypes.c_uint(self.projtypecode),
                    ctypes.c_float(self.pv2_1),
                    np.array(self.pc, dtype=np.float32),
                    ctypes.c_int(self.frontinteg),
                    ctypes.c_uint(self.nbthreads),
                    ctypes.c_uint(self.nbchunks),
                    np.array(self.nerotcntr, dtype=np.float32),
                    np.array(self.nerotang, dtype=np.float32),
                    np.array(self.netranslation, dtype=np.float32),
                    np.array(self.nerotaxis, dtype=np.int32),
                    ctypes.c_int(self.physics),
                    np.array(self.phyparam, dtype=np.float32),
                    ctypes.c_float(self.fracmax),
                    ctypes.c_int(self.runDumpInteg),
                    self.integrand)

        # -- Trick to make things work:
        #    See details before the C++ call
#        tmp = self.pc[0, 1]
#        self.pc[0, 1] = self.pc[1, 0]
#        self.pc[1, 0] = tmp
        




    def getWCS(self):
        """Returns a WCS header for the simulated image."""
        pass


    def dispim(self, minmax=[1e-15, 1e-10], cmap='gist_ncar'):
        """Display the image"""
        
        fig, ax = mu.dispim(self.btot, log=True,
                            minmax=minmax,
                            cmap=cmap)

        return fig, ax
        



