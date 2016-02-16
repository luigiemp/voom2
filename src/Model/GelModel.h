//-*-C++-*-
/*!
  \file EPModel.h

  \brief Implementation of Cardiac Electrophysiology. This is the parent
  class for cardiac EP. Monodomain and Bidomain classes will derive from 
  this class
*/

#ifndef __GelModel_h__
#define __GelModel_h__
#include "Model.h"
#include "GelResult.h"

#include "Spring.h"
#include "EntropicSpring.h"
#include "AngleSpring.h"
#include "BrownianRod.h"

#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Crosslink.h"
#include "Motor.h"
#include "Grid.h"
#include "TwoBodyPotential.h"
#include "IntersectionFinder.h"
#include "PeriodicTie.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "PinchForce.h"
#include "NematicProbTable.h"
#include "GelInput.h"


namespace voom{
  class EPModel : public Model {
  protected:
    //! Diffusion tensor data for each element
    vector<Matrix3d>   _diffusion;

  public:
    //! Constructor
    GelModel(Mesh* Filaments, const GelInput inputFile, const uint NodeDof);

    //! Destructor
    virtual ~GelModel() {;}

    //! Pure virtual compute function call
    void compute(GelResult& R) {;}
  };
}

#endif
