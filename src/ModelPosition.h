
#ifndef MODELPOSITION_H
#define MODELPOSITION_H


/** \file ModelPosition.h
 * \brief Container containing the model positioning information
 */

//! Contains all the model positioning information: basis, rotation, translation
class ModelPosition {
public:
  //! Default constructor
  ModelPosition() {modelbasis=Cbasis();rotation=Cbasis();translation=Cvec();};
  //! Copy
  ModelPosition(const ModelPosition &a) {this->modelbasis=a.modelbasis;this->rotation=a.rotation;this->translation=a.translation;};

  //! Accessors 
  void setBasis(const Cbasis &b) {this->modelbasis=b;};
  void setRotation(const Cbasis &r) {this->rotation=r;};
  void setTranslation(const Cvec &t) {this->translation=t;};

public:
  Cbasis modelbasis; //!< Basis of the model, defined in the absolute coordinate system
  Cbasis rotation; //!< Rotation of the model, defined in the modelbasis coordinate system
  Cvec translation; //!< Translation of the model, defined in the modelbasis coordinate system

};


#endif
