//-*-C++-*-
#ifndef __TorsionalBody_h__
#define __TorsionalBody_h__

#include "Body.h"

namespace voom{

  // TorsionalBody
  class TorsionalBody: public Body {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    TorsionalBody(Mesh* myMesh, State *myState, const vector<int > & BodyNodes, Real TorK, Vector3d Center);

    //! Destructor
    ~TorsionalBody() {};

    //! Accessors and Mutators
    void setTorsionalStiffness(Real TorK) { _torK = TorK; };
    Real getTorsionalStiffness() { return _torK; };
    void setCenter(Vector3d Center) { _center = Center; };
    Vector3d getCenter() { return _center; };
    void updateXprev();
 
    //! Compute function
    void compute(Result* R);

    //! Plot function as 2D elements from Center to BodyNodes
    void writeOutputVTK(const string OutputFile, int step);



  protected:
    //! Field at the beginning of the computation step
    const vector<int > & _bodyNodes;
    vector<Vector3d >    _Xprev;
    Real                 _torK;
    Vector3d             _center;

  };

} // namespace voom

#endif
