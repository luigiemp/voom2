#if !defined(__SCElastic_h__)
#define __SCElastic_h__

#include "ShellMaterial.h"
#include "VoomMath.h"


namespace voom
{

  /*!  Class implementing the simple Helfrich spontaneous curvature
    model for an elastic lipid bilayer.
  */
  
  class SCElastic : public ShellMaterial
  {
  protected:

    double _kC;
    double _kG;
		
    double _C0;
    double _H; // mean curvature
    double _K; // gaussian curvature

  public:

    SCElastic(const double kC, const double kG, const double C0=0.0)
      : _kC(kC), _kG(kG), _C0(C0) {};
    
    void compute(Shellresults & R, const ShellGeometry & defGeom);

    double meanCurvature() const {return _H;}
    double gaussianCurvature() const {return _K;}

    double spontaneousCurvature() const { return _C0; }
    double bendingModulus() const { return _kC; }
    double gaussianModulus() const { return _kG; }

    void setMaterialParameters(const vector<Real > &) {;}
    void setInternalParameters(const vector<Real > &) {;}
    
    vector<Real > getMaterialParameters() {vector<Real > empty; return empty;}
    vector<Real > getInternalParameters() {vector<Real > empty; return empty;}

  };
  
} //namespace voom

#endif //  !defined(__SCElastic_h__)
