//-*-C++-*-
#ifndef __EigenEllipticResult_h__
#define __EigenEllipticResult_h__
#include "EllipticResult.h"

namespace voom {
  //! Eigen Elliptic Result is derived from the general EllipticResult class and contains Eigen storage objects. Will be used by Eigen solver class

  struct EigenEllipticResult: public EllipticResult
  {

    //! Constructor and destructor
    EigenEllipticResult(int PbDoF, int NumDiffMat): _pbDoF(PbDoF)
    {
      // Create results structures
      _stiffness = new SparseMatrix<Real >(_pbDoF, _pbDoF);
      _residual = new VectorXd(_pbDoF);
      
      _Hg = new SparseMatrix<Real >(NumDiffMat, NumDiffMat);
      _Gradg = new VectorXd(NumDiffMat);
    };

    ~EigenEllipticResult() {
      delete _stiffness;
      delete _residual;
      delete _Hg;
      delete _Gradg;
    };

    // Again residual and stiffness are needed by Solver (not good as above)
    VectorXd * getResidual() {return _residual; };
    SparseMatrix<Real > * getStiffness() {return _stiffness; };
    
    
    //! Interface (Mutators)
    // Set functions
    void setResidual(int ind, Real value) {
      (*_residual)[ind] = value;
      }
    
    // Reset function
      void resetResidualToZero() {
	for (int i = 0; i < _pbDoF; i++)
	  (*_residual)[i] = 0.0;
      }

      void resetStiffnessToZero() {
	_stiffness->setZero();
	// for (int k = 0; k < _stiffness->outerSize(); ++k)
	//   for (SparseMatrix<Real>::InnerIterator it(*_stiffness,k); it; ++it) {
	//       coeffRef(it.row(), it.col()) = 0.0;
	//   }
      }



    // Add functions
      void addResidual(int ind, Real value) {
	(*_residual)[ind] += value;
      }

      void addStiffness(int indRow, int indCol, Real value) {
	_stiffness->coeffRef(indRow, indCol) += value;
      };

    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->makeCompressed();
    };
  


    //! InterfacesetZero (Accessors)
      Real getResidual(int ind) {
	return (*_residual)[ind];
      }

      Real getStiffness(int indRow, int indCol) {
	return _stiffness->coeff(indRow, indCol);
      };


  public:
    
    SparseMatrix<Real > *_stiffness;
    VectorXd *_residual;
      
    SparseMatrix<Real > *_Hg;
    VectorXd *_Gradg;
    
    int _pbDoF;
  }; // EigenEllipticResult 

}; // namespace voom

#endif
