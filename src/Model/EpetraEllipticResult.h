//-*-C++-*-
#ifndef __EpetraEllipticResult_h__
#define __EpetraEllipticResult_h__
#include "EllipticResult.h"

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
  //! Epetra Elliptic Result is derived from the general EllipticResult class and contains Epetra storage objects. Will be used by Epetra solver class

  struct EpetraEllipticResult: public EllipticResult
  {
    //! Constructor and destructor
    EpetraEllipticResult(Epetra_MpiComm* comm,
			 const vector<int> & localDoF,
			 const vector<int> & ghostDoF): 
      _localDoFNum(localDoF.size()), _ghostDoFNum(ghostDoF.size()),
      _localDoF(localDoF), _ghostDoF(ghostDoF)
    {
      // Create target and sourcemaps
      _targetMap = new Epetra_Map(-1, localDoF.size(), &localDoF[0], 0, *comm);
      _sourceMap = new Epetra_Map(-1, ghostDoF.size(), &ghostDoF[0], 0, *comm);

      // Create residual and stiffness matrix from source and target maps
      _residual       = new Epetra_Vector(*_targetMap);
      _ghostResidual  = new Epetra_Vector(*_sourceMap);
      _stiffness      = new Epetra_FECrsMatrix(Copy, *_targetMap, 10);
   
      _residual->ExtractView(&_extractedResidual);
      _ghostResidual->ExtractView(&_extractedGhostResidual);
    };

    ~EpetraEllipticResult() {
      delete _targetMap;
      delete _sourceMap;
      delete _residual;
      delete _ghostResidual;
      delete _stiffness;
    };

    // The same Epetra maps are used from Solver to create the solution vector
    // LEP: this is not good because we give full access to private members outside the class
    Epetra_Map * getTargetMap() {return _targetMap; };
    Epetra_Map * getSourceMap() {return _sourceMap; };

    // Again residual and stiffness are needed by Solver (not good as above)
    Epetra_Vector * getResidual() {return _residual; };
    Epetra_Vector * getGhostResidual() {return _ghostResidual; };
    Epetra_FECrsMatrix* getStiffness() {return _stiffness; };
    
    
    //! Interface (Mutators)
    // Set functions
      void setResidual(int LocalIndex, Real value) {
	if ( LocalIndex < _localDoFNum ) {
	  _extractedResidual[LocalIndex] = value;
	} 
	else {
	  _extractedGhostResidual[LocalIndex - _localDoFNum] = value;
	}
      }
    
    // Reset function
      void resetResidualToZero() {
	for (uint i = 0; i < _localDoFNum; i++)
	  _extractedResidual[i] = 0.0;
	for (uint i = 0; i < _ghostDoFNum; i++)
	  _extractedGhostResidual[i] = 0.0;
      }

      void resetStiffnessToZero() {
	_stiffness->PutScalar(0.0);
      }



    // Add functions
      void addResidual(int LocalIndex, Real value) {
	if ( LocalIndex < _localDoFNum ) {
	  _extractedResidual[LocalIndex] += value;
	} 
	else {
	  _extractedGhostResidual[LocalIndex - _localDoFNum] += value;
	}
      }

      void addStiffness(int LocalRow, int LocalCol, Real value) {
	int GlobalRow = 0;
	int GlobalCol = 0;

	if ( LocalRow < _localDoFNum ) { GlobalRow = _localDoF[LocalRow]; }
	else { GlobalRow = _ghostDoF[LocalRow]; }
	
	if ( LocalCol < _localDoFNum ) { GlobalCol = _localDoF[LocalCol]; }
	else { GlobalCol = _ghostDoF[LocalCol]; }

	_stiffness->InsertGlobalValues(1, &GlobalRow, 1, &GlobalCol, &value);
      };

    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->GlobalAssemble();
    };
  


    //! Interface (Accessors)
      Real getResidual(int LocalIndex) {
	if ( LocalIndex < _localDoFNum ) { return _extractedResidual[LocalIndex]; } 
	else { return _extractedGhostResidual[LocalIndex - _localDoFNum]; }
      }

      Real getStiffness(int LocalRow, int LocalCol) {
	Real *Values = 0;
	int *Indices = 0;
	int numEntries = 0;  
	_stiffness->ExtractMyRowView(LocalRow, numEntries, Values, Indices); 
	
	map<int, int > FilledColumns;
	for(int m = 0; m < numEntries; m++) 
	  FilledColumns.insert(make_pair(Indices[m], m)); // Only some columns are full
	
	if ( FilledColumns.find( LocalCol ) != FilledColumns.end() ) {
	  return Values[FilledColumns[LocalCol] ]; }
	else {
	  return 0.0; }
      };
 


  private:
    
    uint _localDoFNum;
    uint _ghostDoFNum;
    const vector<int> & _localDoF;
    const vector<int> & _ghostDoF;

    Epetra_Map * _targetMap;
    Epetra_Map * _sourceMap;
    
    Epetra_Vector * _residual;
    Real * _extractedResidual;
    Epetra_Vector * _ghostResidual;
    Real * _extractedGhostResidual;

    Epetra_FECrsMatrix* _stiffness;
    
  }; // EpetraEllipticResult 

}; // namespace voom

#endif
