//-*-C++-*_
/*!
  \file LMEShape.h
  
  \brief Mesh free Linear Maximum Entropy shape function
*/

#if !defined(__LMEShape_h__)
#define __LMEShape_h__

#include "MFShape.h"

namespace voom {

  class LMEShape: public MFShape {

  public:
    //! Constructor compute N and DN
    LMEShape(const vector<VectorXd> & Nodes, const VectorXd & Point, 
	     const Real beta, const Real tol, const uint MaxIter);

    //! Update recomputes N and DN at new Point
    void update(const VectorXd & Point);

    // Accessors/Mutators to specific LME class members
    Real getBeta() {return _beta; };
    Real getTol() {return _tol; };
    uint getMaxIter() {return _maxIter; };
    
    void setBeta(Real beta) {_beta = beta; };
    void setTol(Real tol) {_tol = tol; };
    void setMaxIter(uint maxIter) {_maxIter = maxIter; };
      
  private:
    // private functions for LME calculation

    //! Calculate shape function j value
    Real LME_pa(const VectorXd & Point,
		const VectorXd & lambda,
		uint j);
  
    //! Calculate vector r
    VectorXd LME_r(const VectorXd & Point, 
		   const VectorXd & lambda);

    //! Calculate Matrix J
    MatrixXd LME_J(const VectorXd & Point, 
		   const VectorXd & lambda);

  private:
    // extra members wrt classic FE shape functions
    Real _beta;
    Real _tol;
    uint _maxIter;
  };

} // namespace voom

#endif
