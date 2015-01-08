//-*-C++-*-
/*!\brief
  A base Ionic material class. Various cell models will derive from this class.
*/

#ifndef __IonicMaterial_h__
#define __IonicMaterial_h__
#include "voom.h"
#include <vector>

// Calling LAPACK Ax=b Solver
extern "C" void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV,
		       double *B, int *LDB, int *INFO);
namespace voom {
  class IonicMaterial{
  protected:
    //! Surface Area to Volume Ratio Xi
    Real        _xi, _capacitance, _gamma;
    vector<Real> _state, _rate;

  public:
    //! Constructor
    IonicMaterial(){;}

    //! Destructor
    ~IonicMaterial(){;}

    //! set Xi - Surface Area to Volume Ratio
    void setXi(const Real Xi){_xi = Xi;}

    //! get State
    vector<Real> getState() const{return _state;}

    //! set State
    void setState(const vector<Real> state){_state = state;}

    //! get Rate
    vector<Real> getRate() const{return _rate;}

    //! set Rate
    Real setRate(const vector<Real> rate){_rate = rate;}

    //! get Capacitance
    Real getCapacitance() const{return _capacitance;}

    //! set Capacitance
    void setCapacitance(const Real capacitance){_capacitance = capacitance;}

    //! Compute Ion virtual function. Return dV/dt to body class
    virtual Real compute(Real Xi, Real C_m, Real dt, 
			     Real Volt, Real istim) = 0;

    /*! 
      Get gamma. Magnitude of active deformation. This depends on the 
      [Ca] concentration.
    */
    virtual Real getGamma() = 0;

    //! Get all internal parameters
    virtual const std::vector<Real>& getInternalParameters(int& nData) const 
    = 0;

    //! Set all internal parameters
    virtual void setInternalParameters(const std::vector<Real>& data) = 0;
  };
}
#endif
