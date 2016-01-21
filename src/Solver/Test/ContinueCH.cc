#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "CahnHilliard.h"
#include "PhaseModel.h"
#include "continuation/continuation.h"
#include <gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>

using namespace voom;
class Problem{
public:
  Problem(){;}
  Problem(const string Nodes, const string connTable){
    _myMesh = new LoopShellMesh(Nodes, connTable);
    _CH = new CahnHilliard(1.0);
    uint NumMat = _myMesh->getNumberOfElements();
    uint NodeDoF = 1;
    _PbDoF = _myMesh->getNumberOfNodes() + 1; // +1 for lag. multiplier
    vector<CahnHilliard *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(_CH);
    int TotNumMatProp = NumMat*1;
    _myModel = new PhaseModel( _myMesh, materials, NodeDoF);
    _myResults = new EigenEllipticResult(_PbDoF, TotNumMatProp);
  }
  int getProbDoF() { return _PbDoF;}

  ~Problem(){
    delete _myModel;
    delete _myMesh;
    delete _CH;
    delete _myResults;
  }
public:
  uint _PbDoF; 
  LoopShellMesh * _myMesh;
  PhaseModel * _myModel;
  EigenEllipticResult * _myResults;
  CahnHilliard* _CH;
};

Problem* prob;
int FV_wrapper(double* indvar, double* args, gsl_vector* out){
  /*
    args = [ --------- PbDof -----------|epsilon ]
  */
  int PbDoF = prob->getProbDoF();
  for (int j=0; j< PbDoF-1; j++)
    (prob->_myModel)->setField(j, args[j] ); //sets the dof
  (prob->_myModel)->setLagMult(args[PbDoF-1]); //set Lag multiplier
  (prob->_CH)->setEpsilon(args[PbDoF] );
  ComputeRequest myRequest = FORCE;
  (prob->_myResults)->setRequest(myRequest);
  (prob->_myModel)->compute( *(prob->_myResults) );
  for(int i=0; i< PbDoF-1; i++){
    out->data[i]=(*( (prob->_myResults)->_residual))(i);    
  }
  out->data[PbDoF-1] = ((prob->_myResults)->_constr);
  return 0;
}

int main(int argc, char** argv) {

  prob = new Problem("T1sphere_nodes.dat","T1sphere_conn.dat");
  int Tot_len = prob->getProbDoF();
  Continuer c;
  c.ptr_F = FV_wrapper;
  c.setNoOfUnknownVar(Tot_len, Tot_len);
  c.setNoOfInArgs(Tot_len + 1); //1 for material parameter
  c.performQuad = false; 
  c.noOfPara = 1; //epsilon
  c.posOfPara=Tot_len;
  c.correctorType = "quasi-newton";
  c.solver = "LU";
  double* trivial_sol=new double[Tot_len + 1];
  for(int i=0;i<Tot_len;i++) trivial_sol[i]= 0.0;
    //for (int i=0; i<2*(NOE-1); i+=2) trivial_sol[i]=0.0;
  trivial_sol[Tot_len]=8.0;

  /*
   gsl_vector* vec =gsl_vector_calloc(Tot_len+1);
   c.F(trivial_sol,vec);
   printf("b=[");
   gsl_vector_fprintf(stdout,vec,"%le");
   printf("];");
  */

  int choice = 1;
  if (choice==1){
    c.setInitialSolution(trivial_sol);
    c.Continue();
  }
  
  if(choice==2){
    c.BPFileName="secondary_BP_eps";
    c.DataFileName="branch3.br";
    c.stepSize= 0.1;
    c.load_BP("BP_Filename",3);
    c.Continue();
  }

  
      c.setInitialSolution(trivial_sol);

      /*
      //Matrix and Vector definitions needed for SVD
      gsl_matrix* jac_b=gsl_matrix_calloc(c.jacobian->size2,c.jacobian->size2);
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,c.jacobian,c.jacobian,0.0,jac_b);
   
      gsl_matrix* V=gsl_matrix_alloc(jac_b->size1,jac_b->size2);
      gsl_vector* S=gsl_vector_alloc(jac_b->size2);
      gsl_vector* work=gsl_vector_alloc(jac_b->size2);
      gsl_vector* tgt = gsl_vector_alloc(V->size1);

    
      gsl_linalg_SV_decomp(jac_b,V,S,work);	
      gsl_matrix_transpose(V);
      gsl_matrix_get_row(tgt,V,V->size1-1);
      gsl_vector_fprintf(stdout,S,"%le");	
      cout << endl;
      gsl_vector_fprintf(stdout,tgt, "%le");
      //pprintMat(c.jacobian);
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);
      */

      /*
      gsl_matrix* jac_b=gsl_matrix_calloc(c.jacobian->size2,c.jacobian->size2);
      //gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,c.jacobian,c.jacobian,0.0,jac_b);
      gsl_matrix_memcpy(jac_b,c.jacobian);
      //gsl_matrix_set_identity(jac_b);
      //gsl_matrix_set(jac_b,0,0,-2);
      gsl_eigen_symmv_workspace* work =gsl_eigen_symmv_alloc(jac_b->size2);
      gsl_matrix* evec = gsl_matrix_alloc(jac_b->size1,jac_b->size2);
      gsl_vector* eval=gsl_vector_alloc(jac_b->size2);
      gsl_eigen_symmv(jac_b, eval, evec, work);
      gsl_eigen_symmv_sort(eval,evec, GSL_EIGEN_SORT_VAL_DESC);
      gsl_vector_fprintf(stdout,eval,"%le");	
      cout << endl;
      //gsl_vector_fprintf(stdout,tgt, "%le");
      //pprintMat(c.jacobian);
      gsl_eigen_symmv_free(work);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      */
  }

