/*! \file Clos.h 
 * \brief Line of Sight definition
 *
 *  
 */


#ifndef _CLOS_H_
#define _CLOS_H_

//! Line Of Sight integration definition
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
	  	int nbp;        //!> Number of steps along the LOS
  		float sstart;   //!> Start position distance from reference point (Obs or Plane of Sky)
  		float send;     //!> End position distance from reference point (Obs or Plane of Sky)
  		float ds;       //!> Integration step along the LOS
	
};


#endif	//_CLOS_H_

