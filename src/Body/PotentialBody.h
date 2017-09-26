//-*-C++-*-
#ifndef __PotentialBody_h__
#define __PotentialBody_h__

#include "Body.h"
#include "PairMaterial.h"

namespace voom{

  // PotentialBody
  class PotentialBody: public Body {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    PotentialBody(Mesh* myMesh, State *myState, PairMaterial* Potential, const vector<int > & BodyNodes, Real SearchR);

    //! Destructor
    ~PotentialBody() {};

    //! Set particle-particle interaction
    void getUniqueMeshNodes();
    void setInteractions();
    void setNodesToEle();
    void computeNormals();
    void setSearchR(Real SearchR) { _searchR = SearchR; };
    Real getSearchR() { return _searchR; };

    //! Compute function
    void compute(Result* R);

    //! Plot function using getX from mesh
    void writeOutputVTK(const string OutputFile, int step);



  protected:
    //! Field at the beginning of the computation step
    PairMaterial* _potential;
    const vector<int > &    _bodyNodes;
    Real _searchR;
    set<int >               _uniqueNodes;
    map<int, vector<int > > _nodesToEle;
    vector<vector<int > >   _interactions;
    map<int, Vector3d     > _normals;

  };

} // namespace voom

#endif
