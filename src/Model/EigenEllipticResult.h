//-*-C++-*-
#ifndef __EigenEllipticResult_h__
#define __EigenEllipticResult_h__
#include "EllipticResult.h"

namespace voom {
  //! Eigen Elliptic Result is derived from the general EllipticResult class and contains Eigen storage objects. Will be used by Eigen solver class

  struct EigenEllipticResult: public EllipticResult
  {

    //! Constructor and destructor
    EigenEllipticResult(int PbDoF, int NumMatProp): _pbDoF(PbDoF), _numMatProp(NumMatProp)
    {
      // Create results structures
      _stiffness = new SparseMatrix<Real >(_pbDoF, _pbDoF);
      _residual = new VectorXd(_pbDoF);
      
      _Hg = new SparseMatrix<Real >(NumMatProp, NumMatProp);
      _Gradg = new VectorXd(NumMatProp);
    };

    ~EigenEllipticResult() {
      delete _stiffness;
      delete _residual;
      delete _Hg;
      delete _Gradg;
    };

    // Return NumMatProp and PbDoF
    int getNumMatProp() {
      return _numMatProp;
    }

    int getPbDoF() {
      return _pbDoF;
    }
    
    
    //! Interface (Mutators)
    // Set functions
    void setResidual(int ind, Real value) {
      (*_residual)(ind) = value;
    }
    
    // Reset function
    void resetResidualToZero() {
      *_residual = VectorXd::Zero(_pbDoF);
    }

    void resetStiffnessToZero() {
      _stiffness->setZero();
    }

    void resetGradgToZero() {
      *_Gradg = VectorXd::Zero(_numMatProp);
    }

    void resetHgToZero() {
      _Hg->setZero();
    }



    // Add functions
    void addResidual(int ind, Real value) {
      (*_residual)(ind) += value;
    }
    
    void addStiffness(int indRow, int indCol, Real value) {
      _stiffness->coeffRef(indRow, indCol) += value;
    };

    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->makeCompressed();
    };

    void setStiffnessFromTriplets(vector<Triplet<Real > > & B) {
      _stiffness->setFromTriplets(B.begin(), B.end());
    };

    void addGradg(int ind, Real value) {
      (*_Gradg)(ind) += value;
    }

    void addHg(int indRow, int indCol, Real value) {
      _stiffness->coeffRef(indRow, indCol) += value;
    }

    void setHgFromTriplets(vector<Triplet<Real > > & B) {
      _Hg->setFromTriplets(B.begin(), B.end());
    };
  


    //! InterfacesetZero (Accessors)
    Real getResidual(int ind) {
      return (*_residual)(ind);
    }

    Real getStiffness(int indRow, int indCol) {
      return _stiffness->coeff(indRow, indCol);
    };

    Real getGradg(int ind) {
      return (*_Gradg)(ind);
    }

    Real getHg(int indRow, int indCol) {
      return _Hg->coeff(indRow, indCol);
    };


  public:
    
    SparseMatrix<Real > *_stiffness;
    VectorXd *_residual;
    
    SparseMatrix<Real > *_Hg;
    VectorXd *_Gradg;
    
    int _pbDoF;
    int _numMatProp;
  }; // EigenEllipticResult 
  
}; // namespace voom

#endif
