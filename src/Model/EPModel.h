//-*-C++-*-
/*!
  \file EPModel.h

  \brief Implementation of Cardiac Electrophysiology. This is the parent
  class for cardiac EP. Monodomain and Bidomain classes will derive from 
  this class
*/

#ifndef __EPModel_h__
#define __EPModel_h__
#include "Model.h"
#include "EPResult.h"

namespace voom{
  class EPModel : public Model {
  protected:
    //! Diffusion tensor data for each element
    vector<Matrix3d>   _diffusion;

  public:
    //! Constructor
    EPModel(Mesh* myMesh, const string inputFile, const uint NodeDof);

    //! Destructor
    virtual ~EPModel() {;}

    //! Pure virtual compute function call
    void compute(EPResult& R) {;}
  };
}

#endif
