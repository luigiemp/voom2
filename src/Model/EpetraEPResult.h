//-*-C++-*-
#ifndef __EpetraEPResult_h__
#define __EpetraEPResult_h__
#include "EPResult.h"

// Epetra Headers
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

namespace voom {
  struct EpetraEPResult: public EPResult{
    //! Constructor
    EpetraEPResult(Epetra_MpiComm* comm,
		   const vector<int> & localDoF,
		   const vector<int> & ghostDoF): 
      _localDoFNum(localDoF.size()), _ghostDoFNum(ghostDoF.size()),
      _localDoF(localDoF), _ghostDoF(ghostDoF)
    {
      // Create target and sourcemaps
      _targetMap = new Epetra_Map(-1, localDoF.size(), &localDoF[0], 0, *comm);
      _sourceMap = new Epetra_Map(-1, ghostDoF.size(), &ghostDoF[0], 0, *comm);
      
      // Create Epetra vector and Crs matrices
      _diffusionRV       = new Epetra_Vector(*_targetMap);
      _ghostDiffusionRV  = new Epetra_Vector(*_sourceMap);
      _ionicCurrent      = new Epetra_Vector(*_targetMap);
      _ghostIonicCurrent = new Epetra_Vector(*_sourceMap);
      _diffusionLHS      = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      _diffusionRHS      = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      _ionicLHS          = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      _ionicRHS          = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
      
      _ionicCurrent->ExtractView(&_extractedIonicCurrent);
      _ghostIonicCurrent->ExtractView(&_extractedGhostIonicCurrent);
    };
    
    //! Destructor
    ~EpetraEPResult() {
      delete _targetMap;
      delete _sourceMap;
      delete _diffusionRV;       delete _ghostDiffusionRV;
      delete _ionicCurrent;      delete _ghostIonicCurrent;
      delete _diffusionLHS;      delete _diffusionRHS;
      delete _ionicLHS;          delete _ionicRHS;
    };

    // The same Epetra maps are used from Solver to create the solution vector
    Epetra_Map * getTargetMap() {return _targetMap; };
    Epetra_Map * getSourceMap() {return _sourceMap; };

    //! Interface functions. Virtual functions implemented.
    //! Add to Ionic Current Epetra Vector
    void addIonicCurrent(const int localIndex, const Real value) {
      if ( localIndex < _localDoFNum ) 
	_extractedIonicCurrent[localIndex] += value;
      else
	_extractGhostIonicCurrent[localIndex - localDoFNum] += value;
    }
    
    //! Add to diffusionLHS matrix
    void addDiffusionLHSMatrix(const int localRow, const int localCol,
			       const Real value) {
      _addToMatrix(_diffusionLHS, localRow, localCol, value);
    }

    //! Add to diffusionRHS matrix
    void addDiffusionRHSMatrix(const int localRow, const int localCol,
			       const Real value) {
      _addToMatrix(_diffusionRHS, localRow, localCol, value);
    }

    //! Add to ionicLHS matrix
    void addIonicLHSMatrix(const int localRow, const int localCol,
			   const Real value) {
      _addToMatrix(_ionicLHS, localRow, localCol, value);
    }

    //! Add to ionicRHS matrix
    void addIonicRHSMatrix(const int localRow, const int localCol,
			   const Real value) {
      _addToMatrix(_ionicRHS, localRow, localCol, value);
    }

    //! Assemble all encessary matrices
    void finalizeMatrixAssembly() {
      _diffusionLHS->GlobalAssemble();
      _diffusionRHS->GlobalAssemble();
      _ionicLHS->GlobalAssemble();
      _ionicRHS->GlobalAssemble();
    }

  private:
    
    uint _localDoFNum;
    uint _ghostDoFNum;
    const vector<int> & _localDoF;
    const vector<int> & _ghostDoF;

    Epetra_Map * _targetMap;
    Epetra_Map * _sourceMap;
    
    // Intermediate RHS(b in Ax=b) for Diffusion equation solve
    Epetra_Vector * _diffusionRV;
    Epetra_Vector * _ghostDiffusionRV;

    // Ionic currents
    Epetra_Vector* _ionicCurrent;
    Real* _extractedIonicCurrent;
    Epetra_Vector* _ghostIonicCurrent;
    Real* _extractedGhostIonicCurrent;

    Epetra_FECrsMatrix* _diffusionLHS;
    Epetra_FECrsMatrix* _diffusionRHS;
    Epetra_FECrsMatrix* _ionicLHS;
    Epetra_FECrsMatrix* _ionicRHS;

    //! Private function to add entires in Epetra FECrs matrix. 
    void _addToMatrix(Epetra_FECrsMatrix* const A, const int localRow,
		      const int localCol, const Real value) {
      int GlobalRow = 0;
      int GlobalCol = 0;
      
      if ( localRow < _localDoFNum ) GlobalRow = _localDoF[localRow]; 
      else GlobalRow = _ghostDoF[localRow]; 
	
      if ( localCol < _localDoFNum ) GlobalCol = _localDoF[localCol]; 
      else GlobalCol = _ghostDoF[localCol]; 
      
      A->InsertGlobalValues(1, &GlobalRow, 1, &GlobalCol, &value);
    }
  };
}
#endif
