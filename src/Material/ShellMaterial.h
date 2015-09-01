#if !defined(__ShellMaterial_h__)
#define __ShellMaterial_h__

#include <vector>
#include "voom.h"
#include "ShellGeometry.h"

namespace voom
{

  /*!  Base Class template for a thin-shell material object.  This
    template can be combined with generic 3-D materials for use with
    thin shells.
  */

  class ShellMaterial
  {
  public:
    struct Shellresults
    { 
      Shellresults() 
      {
	W = 0.0;
	n.resize(2,Vector3d::Zero());
	m.resize(2,Vector3d::Zero());
	request = 0;
      };
    
      double W; // Energy Density 
      vector< Vector3d > n;   // stress resultants
      vector< Vector3d > m;   // moment resultants
      double stretch;
      Matrix3d cauchy; // Cauchy Stress Tensor
      int request;
    
    }; // struct FKresults
    
    //! Constructor
    // No material can be initialized without ID ****** TODO : Must be made private *****
    ShellMaterial(){;};

    ShellMaterial(int MatID): _matID(MatID) {;}

    //! virtual destructor
    virtual ~ShellMaterial() { ; }
    
    //! Clone
    //virtual ShellMaterial* clone() const = 0;
    //! Compute function
    virtual void compute(Shellresults & R, const ShellGeometry & defGeom) = 0;

    //! Consistency Check for all Mecahnics Material Classes
    void checkConsistency(Shellresults & R, const ShellGeometry & defGeom,
			  const Real h = 1.0e-7, const Real tol = 1.0e-6);

    const ShellGeometry& refShellGeometry() const { return _referenceGeometry; }
    void setRefGeometry( const ShellGeometry& g ) {_referenceGeometry = g;}
    
    //! SetMaterialParameters function
    virtual void setMaterialParameters(const vector<Real > &) = 0;
    virtual void setInternalParameters(const vector<Real > &) = 0;
    //virtual void setRegularizationParameters(const vector<Real > &) = 0;

     //! GetMaterialParameters function
    virtual vector<Real > getMaterialParameters() = 0;
    virtual vector<Real > getInternalParameters() = 0;
    //virtual vector<Real > getRegularizationParameters() = 0;

    // Get Material ID
    int getMatID() {
      return _matID;
    }

  protected:
    int _matID;
    ShellGeometry _referenceGeometry;
  };

} //namespace voom

#endif //  !defined(__ShellMaterial_h__)
