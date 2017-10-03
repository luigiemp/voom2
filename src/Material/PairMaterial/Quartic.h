/*! 
  \file Quartic.h
  \brief Interface for a quartic spring potential
*/

#ifndef _QUARTIC_H_
#define _QUARTIC_H_

#include "PairMaterial.h"

namespace voom {
  
  class Quartic : public PairMaterial
  {
  public: 
    // Constructors/destructors:
  Quartic(): _k(1.0), _r0(0.0) {};
  Quartic(Real k, Real R0): _k(k), _r0(R0) {};
    Quartic(Quartic* BaseMaterial): _k(BaseMaterial->_k), _r0(BaseMaterial->_r0) {};

    // Clone
    virtual Quartic* clone() const {
      return new Quartic(*this);
    }

    void setMaterialParameters(const vector<Real > & parameters) {
      _k  = parameters[0];
      _r0 = parameters[1];
    }

    vector<Real > getMaterialParameters() {
      vector<Real > MatProp(2, 0.0);
      MatProp[0] = _k;
      MatProp[1] = _r0;
      return MatProp;
    }
    
    // Operators
    //! Based on particle distance r, calculates state of material
    void compute(PairMresults & R, Real r);

    
    
  private:
    //! Data
    Real _k;
    Real _r0;
	
  }; // class Quartic
  
}
#endif // _QUARTIC_H_
