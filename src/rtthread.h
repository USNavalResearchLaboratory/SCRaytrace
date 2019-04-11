
#ifndef RTTHREAD_H
#define RTTHREAD_H

/** \file rtthread.h
 * \brief Raytrace with threads to speed up on multi-core.
 */
 
// All angles are in radian
// All distances are in Rsun
extern "C" int rtthread(int sx,                         //!> image size, x axis [in]
                        int sy,                         //!> image size, y axis [in]
                        float fovpix,                   //!> fov angle of one pixel in rad [in]
                        float *obspos,                  //!> [x, y, z] position of the observer in the Sun basis [in]
                        float *obsang,                  //!> [ax, ay, az] orientation of the observer. z is the optical axis. Rotation order is z, y, then x. [in]
                        float *nepos,                   //!> [x, y, z] position of the Ne reference in the Sun basis [in]
                        float *neang,                   //!> [ax, ay, az] orientation of the Ne [in]
                        int losnbp,                     //!> number of steps for the integration along the LOS [in]
                        float *losrange,                //!> [lstart, lend] range for the integration along the LOS in Rsun. The origin of the LOS is the orthogonal projection of the Sun cntr on that LOS. [in]
                        int modelid,                    //!> model id number [in]
                        float *btot,                    //!> Total brightness image (for Thomson scattering physiscs) [out]
                        float *bpol,                    //!> Polarized brightness image (for Thomson scattering physiscs [out]
                        float *netot,                   //!> Electron density image (for Thomson scattering physiscs) [out]
                        float *pmodparam,               //!> parameters of the model [in]
                        float *crpix,                   //!> [x, y] pixel boresight center of the image
                        int quiet,                      //!> quiet mode if set to 1 [in]
                        int neonly,                     //!> Only compute the electron density if set to 1 [in]
                        float *hlonlat,                 //!> [Hlon, Hlat, Hrot] heliographic lon and lat of the center of the disk, rotation angle corresponding to the projection of the north pole, counterclockwise [in].
                        float occrad,                   //!> Occulter radius [Rsun]. The integration in not performed within that disk [in]
                        float limbdark,                 //!> limb darkening coeff: default 0.58 [in]
                        float *obslonlat,               //!> [lon, lat, height] position of the observer in Carrington coordinate. If set, then obspos is ignored. The optical axis always points toward the Sun center. Use obsang to change telescope orientation. Note that obslonlat=[0,0,215] correspond to obspos=[0,0,215] and obsang=[!pi,0,0]: this means that the Carrington coordinate origin on the Solar sphere (lon,lat,height)=(0,0,1) is located at (x,y,z)=(0,0,1), with Ox pointing to solar north and Oy pointing to (lon,lat)=(3*!pi/2,0) [in]
                        int obslonlatflag,              //!> Set to 1 to use obslonlat to position the observer, instead of obspos [in]
                        unsigned int projtypecode,      //!> Projection type. (see Calabretta and Greisen,  Representations of celestial coordinates in FITS, A&A 395, 1077-1122(2002)), ARC : Zenithal equidistant (default), TAN : Gnomonic, SIN : Slant orthographic, AZP : Zenithal perspective. [in]
                        float pv2_1,                    //!> mu parameter for the AZP projection [in]
                        float *pc,                      //!> wcs pc[2, 2] matrix: default is unit matrix [in]
                        int frontinteg,                 //!> Set to 1 so that the origin of the LOS is taken at the observer: if set, the losrange parameters must both be positive. [in]
                        unsigned int nbthreads,         //!> Number of threads to run in parallel. [in]
                        unsigned int nbchunk,           //!> Number of chuncks. Use with nbthread. If set to a value less than 2, the threads are  launched by lines of sight. If nbchunks >= 2, the threads are launched by chunk of the image. Ballance nbthreads and nbchunks for optimal performances. [in]
                        float *nerotcntr,               //!> [x, y, z] center of rotation of the Ne model, in the Ne basis [in]
                        float *nerotang,                //!> [ax, ay, az] rotation of the Ne model around the nerotcntr, in the Ne basis [in]
                        float *netranslation,           //!> [tx, ty, tz] translation vector of the Ne model, in the Ne basis [in]
                        int *nerotaxis,                 //!> [axid1, axid2, axid3] axis id corresponding to the nerotang rotation angles. 1: X, 2: Y, 3: Z. Default is [3,2,1]. [in]
                        int physics,                    //!> type of physics to perform the raytracing [in]
                        float *phyparam,                //!> Extra parameters required depending on the chosen physics [in]
                        float fracmax,                  //!> Set it to the fraction of the maximum total B per LOS in order to compute the distance to that fraction of brightness. Disabled if set to 0 [default]. The distance is returned in the bpol image, in Rsun. [in]
                        int runDumpInteg,               //!> Set if you want to save all the integration points in the integrand variable. Note that this feature can require the allocation of a large amount of memory, typically a floating array of imsize[0] x imsize[1] x losnbp. Use this feature only if you have enough free memory. [in]
                        float *pIntegrand);             //!> Contains the all the integration points if rundumpinteg is set [out]


int rtthreadtest();

extern "C" void footest();
extern "C" void testPassPython(float* array, int N, float f);
extern "C" int rtthreadtestExt();

//! Raytracing with threads to speed up on multi-core computers
extern "C" int rtthreadidl(int argc, void **argv);



#endif
