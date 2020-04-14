/*! \file CModelBase.h 
 * \brief Base class for all models
 *
 *  
 */


#ifndef _CMODELBASE_H_
#define _CMODELBASE_H_

#include <vector>
#include <string>
#include "Cvec.h"



//! Forward declaration of class Scene
class Scene;



//! Model default parameter description structure
struct moddefparam
{
    moddefparam() {};
    moddefparam(std::string,std::string,std::string,std::string);

    std::string keyword; //!< Name of the parameter
    std::string def; //!< Parameter default value
    std::string desc; //!< Description of the parameter
    std::string units; //!< Units of the parameter
};

//! Base model class, purely virtual
class CModelBase
{
public:
    CModelBase(){};
    virtual ~CModelBase(){};

    //! Initialization of the parameters
    //! \param pparam Points to the vector of parameters of the model.
    virtual void initParam(float* pparam)
    {
        this->pparam=pparam;
    };

    
    /** Compute the density at position v
     \param v Position where to calculate the electron density in the density model coordinate system.
     \param temperature Returns the temperature at position v, when calculated in the model.
     \return The electron density at position v.
    */
    virtual float Density(const Cvec &v,float &temperature)
    {
        temperature=0.;
        return 0.;
    };
    
    
    /** Compute the density at position v
     \param v Position where to calculate the electron density in the density model coordinate system.
     \return The electron density at position v.
    */
    virtual float Density(const Cvec &v)
    {
        return 0.;
    };

    
    /** Save the density parameters in a file for IDL
      \param vp Vector of moddefparam containing the default value and a description of each parameter of the model.
      \param flagcase Flag relative to the following information:
    		\li 0x1 : default parameters undefined.
    		\li 0x2 : model for testing purpose.
    		\li 0x4 : no parameters needed.
    */
    virtual void dumpDefaultParamForIDL(std::vector<moddefparam> & vp,int & flagcase);

    
    /** Link the density to the Scene that's rendered.
     * \param pparentscene Pointer to a Scene object.
     */
    void setParentScene(Scene *pparentscene) {this->pparentscene=pparentscene;}
    
protected:
    float *pparam;
    Scene *pparentscene;
};

//! model selection function
CModelBase* modelselect(int modelid);



#endif	//_CMODELBASE_H_



