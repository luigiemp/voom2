#include "MechanicsSolver.h"

namespace voom {
  template<class nodeType, unsigned int dim>
  MechanicsSolver<nodeType, dim>::MechanicsSolver(Epetra_Mpicomm* comm,
						  Model<nodeType, dim>* mod):
    Model<nodeType,dim>(comm, mod) {
    // Create target and sourcemaps
    _sourceMap = new Epetra_Map(-1, localID.size(), &localID[0], 0, *_comm);
    _targetMap = new Epetra_Map(-1, ghostID.size(), &ghostID[0], 0, *_comm);
    _exporter  = new Epetra_Export(*sourceMap, *targetMap);
    
    //! Create Epetra vectors now
    _U              = rcp(new Epetra_Vector(*_targetMap));
    _ghostU         = rcp(new Epetra_Vector(*_sourceMap));
    _FInternal      = rcp(new Epetra_Vector(*_targetMap));
    _ghostFInternal = rcp(new Epetra_Vector(*_sourceMap));
    _FExternal      = rcp(new Epetra_Vector(*_targetMap));
    _ghostFExternal = rcp(new Epetra_Vector(*_sourceMap));

    // Stiffness matrix
    _K              = rcp(new Epetra_FECrsmatrix(Copy, *_targetMap, 10));
  }

  // Solver the model
  template<class nodeType, unsigned int dim>
  MechanicsSolver<nodeType, dim>::Solver() {
    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& nlParams = *nlParamsPtr.get();
    // Set the nonlinear solver method
    nlParams.set("Nonlinear Solver", "Line Search Based");
    
    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
		    NOX::Utils::OuterIteration + 
		    NOX::Utils::OuterIterationStatusTest + 
		    NOX::Utils::InnerIteration +
		    NOX::Utils::Parameters + 
		    NOX::Utils::Details + 
		    NOX::Utils::Warning);
    
    // Create printing utilities
    NOX::Utils utils(printParams);
    
    // Sublist for line search 
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", "Full Step");
    
    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");
    
    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 800);  
    lsParams.set("Tolerance", 1e-4); 
    lsParams.set("Preconditioner", "None");
    //lsParams.set("Preconditioner", "Ifpack");
    lsParams.set("Max Age Of Prec", 5); 
    
    // Create the interface between the test problem and the nonlinear solver
    
    Teuchos::RCP<Nox_Interface> interface = Teuchos::rcp(noxface);
    
    // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
    // 1. User supplied (Epetra_RowMatrix)
    //Teuchos::RCP<Epetra_RowMatrix> Analytic = Problem.getJacobian();
    // 2. Matrix-Free (Epetra_Operator)
    Teuchos::RCP<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
    // 3. Finite Difference (Epetra_RowMatrix)
    Teuchos::RCP<NOX::Epetra::FiniteDifference> FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln));					
    
    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = MF;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
													      printParams, lsParams,iReq, iJac, MF, noxSoln));
    
    
    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln, linSys)); 
    
    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::NormF> relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-2));
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RCP<NOX::StatusTest::Combo> converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(relresid);
    converged->addStatusTest(wrms);
    converged->addStatusTest(update);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);				
    
    
    Teuchos::RCP<Teuchos::ParameterList> finalParamsPtr = nlParamsPtr;
    
    // Create the method
    solver = NOX::Solver::buildSolver(grp, combo, finalParamsPtr);
    
    NOX::StatusTest::StatusType status = solver->solve();
    
    if (status == NOX::StatusTest::Converged)
      utils.out() << "Test Passed!" << std::endl;
    else {
      if (MyPID==0) 
	utils.out() << "Nonlinear solver failed to converge!" << std::endl;
    }
    
    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();		
    
    //then this should be send to print method
    //or pushed in the nodes before printing 
  }
}
