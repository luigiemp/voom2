//-*-C++-*-
#ifndef _EPSolver_h_
#define _EPSolver_h_
#include "voom.h"
#include "Model.h"
#include "EpetraEPResult.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_CrsMatrix.h"

// Aztec gives us Iterative Solvers
#include "AztecOO_config.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

namespace voom{
  //! Specify the Solution technique to be used
  enum Method {EXPLICIT, BACKWARDEULER, CRANKNICHOLSON};

  //! Specify the Lumping Scheme to be used (C = Consistent, L = Lumped)
  enum LumpingScheme {CCC, CCL, CLC, CLL, LCC, LCL, LLC, LLL};
 
  class EPSolver {
  protected:
    //! Pointer to Model being used
    Model* _myModel;

    //! Model Results Struct
    EpetraEPResult* _result;

    // Tolerances and Number of Iterations:
    //! Linear System Solver tolerance
    Real _linearSolverTolerance;

    //! Linear Solver max Iterations
    int _linearSolverMaxIterations;

    //! Time Steps
    Real _ionicTimeStep;
    Real _diffusionTimeStep;

    //! Max NR iteration
    int _maxIter;

    //! Lumping Scheme
    LumpingScheme _lumpingScheme;

    //! Epetra Objects
    Epetra_MpiComm* _comm;
    Epetra_Export*  _exporter;

    //! Diffusive Current
    Epetra_Vector* _diffusiveCurrent;

    //! Ghost Diffusive Current  
    Epetra_Vector* _ghostDiffusiveCurrent;

  public:

    EPSolver(Epetra_MpiComm *comm, const Model& model, EpetraEPResult& Result,
	             const Mesh* mesh);
    
    //! Destructor
    ~EPSolver() {
      delete _exporter;
      delete _diffusiveCurrent;
      delete _ghostDiffusiveCurrent;
    }

    /*! 
      Pure virtual function. Will be implemented in derived classes
      @param dtDiff Diffusion time step. Default is 0.1 ms
     */
    virtual void Solve(const Real dtDiff = 0.1) = 0;

    //! Set Ionic Time Step
    void setIonicTimeStep(const Real dt) { _ionicTimeStep = dt; }

    //! Set Diffusion Time Step
    void setDiffusionTimeStep(const Real dt) {_diffusionTimeStep = dt;}

    /*! 
      Write Output (Binary File)
      @param time Time at which output is written. 
      @param path Directory where output files are written
    */
    void WriteOutput(const Real time, std::string path = "");

    //! Initialize the Voltage Epetra vector to some value
    void InitializeVoltage(const Real value){
      // _Voltage->PutScalar( value );
      // _ghostVoltage->PutScalar( value );
    }

    //! Consistency Check
    void checkConsistency();

    //! Set max iteration value
    void setMaxIteration(const int maxIter) {
      _maxIter = maxIter;
    }

    //! Gateway to EPSolver Class. Use this function to create EPSolver objects.
    static EPSolver* New(Epetra_MpiComm *comm, Model& model, 
			 EpetraEPResult& EPResult,
			 const LumpingScheme LScheme, const Method SolverType);
  };
}
#endif
