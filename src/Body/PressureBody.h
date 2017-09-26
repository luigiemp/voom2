//-*-C++-*-
#ifndef __PressureBody_h__
#define __PressureBody_h__

#include "Body.h"

namespace voom{

  // PressureBody
  class PressureBody: public Body {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    PressureBody(Mesh* myMesh, State *myState, Real Pressure);

    //! Destructor
    ~PressureBody() {};

    //! Initialize field - Only of the nodes related to this body
    // From constant value
    void initializeField(Real fact = 1.0);
    void setPrevField() { _myState->getPhi(_prevPhi); };

    //! Compute function
    void compute(Result* R);

    

  protected:
    //! Field at the beginning of the computation step
    vector<Real > _prevPhi;

    Real _pressure;
  };

} // namespace voom

#endif
