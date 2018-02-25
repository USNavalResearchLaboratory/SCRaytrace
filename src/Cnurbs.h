//
// File: Cnurbs.h
// $Id: Cnurbs.h,v 1.3 2009/02/09 20:51:07 thernis Exp $

#ifndef _CNURBS_H_
#define _CNURBS_H_

#include <vector>
#include "CControlPoint.h"

#define VMIN 0.00001
#define VMAX (1.-VMIN)
#define VDIFF (VMAX-VMIN)

//! Nurbs curves implementation
class Cnurbs
{
	public:
		 ~Cnurbs();
	
		// Cnurbs interface
		//! define all the parameters of the nurbs
		Cnurbs(unsigned int Ncontrol0,std::vector<float> KnotVec0,unsigned int Nbvertices0,std::vector<CControlPoint> cp0);
		// member function declarations...
		//! simple print of the basis functions
		void printBasis(unsigned int Ko);
		//! simple print of the points of the curve
		void printCurve();
		//! get the curve 
		std::vector<Cvec> getCurve();
	
	protected:
		//! run computation of the bspline basis
		void calcBasis();
		//! allocate memory for the bspline basis array
		int allocate3Darray();
		//! free the bspline basis array
		int delete3Darray();
		//! Calculate the nurbs curve points
		void calcCurve();
	
	protected:
		// Cnurbs variables
		unsigned int Ncontrol; //!< Number of control points
		unsigned int Korder; //!< Order of the NURBS

		unsigned int KnotSize; //!< size of the knot vector	
		std::vector<float> KnotVec; //!< knot vector

		unsigned int Nbvertices; //!< Number of vertices for the curve
		std::vector<float> t; //!< parametric coordinate
	
		std::vector<CControlPoint> cp; //!< The vector of control points

		std::vector<Cvec> curve;  //!< the points of the nurbs curve 
	
		float*** Nik; 	//!< array containing the bspline basis
						//!< Nik[Nbvertices][KnotSize][Korder]

};


#endif	//_CNURBS_H_



/*
* $Log: Cnurbs.h,v $
* Revision 1.3  2009/02/09 20:51:07  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
