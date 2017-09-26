//-*-C++-*-
#ifndef __PreferredNormalBody_h__
#define __PreferredNormalBody_h__

#include "Body.h"

namespace voom{

  // PreferredNormalBody
  class PreferredNormalBody: public Body {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    PreferredNormalBody(Mesh* myMesh, State *myState, Real Krot);

    //! Destructor
    ~PreferredNormalBody() {};

    //! Initialize field
    void initializeField(Real fact = 1.0);

    //! Initialize normals to initial value
    void initializeNormals();
   
    //! Compute function
    void compute(Result* R);

    

  protected:
    //! Field at the beginning of the computation step
    vector<Vector3d > _initNormals;

    Real _krot;
  };

} // namespace voom

#endif
