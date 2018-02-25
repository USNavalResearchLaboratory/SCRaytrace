//
// File: Clos.h
// $Id: Clos.h,v 1.3 2009/02/09 20:51:04 thernis Exp $

#ifndef _CLOS_H_
#define _CLOS_H_

//! Line Of Sight definition
class Clos
{
	public:
		Clos();
		 ~Clos();
	
		// Clos interface
		Clos(int n,float st,float en);
    	Clos(const Clos &a);

	    void setLOS(int n,float st,float en);
	
	public:
		// Clos variables
	  	int nbp;
  		float sstart;
  		float send;
  		float ds;
	
};


#endif	//_CLOS_H_


/*
* $Log: Clos.h,v $
* Revision 1.3  2009/02/09 20:51:04  thernis
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
