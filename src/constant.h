/*! \file constant.h 
 * \brief Defines the constants used in the software.
 *
 *  
 */


#ifndef _CONSTANT_H
#define _CONSTANT_H

#define PI 3.1415926535897932384626433832795028841971693993751058209749
#define TWOPI (2*PI)
#define DTOR (PI/180.)
#define RADEG (1./DTOR)
//! Rsun, in cm
#define RSUN_CM 695500e5
//! Rsun, in m
#define RSUN_M 695500000
// I previously used 696000000 m
//! one astronomical unit, in m
#define ONEAU_M 149597870691
//! one astronomical unit, in Rsun
#define ONEAU_RSUN (ONEAU_M / RSUN_M)

//! Useful to print variables
#define printvar(X) std::cout << #X << " : " << X << std::endl;
#define printvarnoend(X) std::cout << #X << " : " << X << ", " ;



//! Minimum limit for the radius of integration around the sun
#define LIMBLIMIT 1.001
//! Maximum recursion
#define MAXRECURSION 10
//! Rsun, in Rsun (very useful !)
#define RSUN 1.
//! limb darkening
#define U 0.58
//! Thomson scattering integral constant factor 
#define constfactor(u) (1.24878E-25/(1-u/3)*RSUN_CM)
//! Minimum electron density
#define MIN_ELECTRON_DENSITY 1E-1





#endif /* _CONSTANT_H */

