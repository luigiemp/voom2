#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "LandauBrazovskii.h"
#include "LBModel.h"
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
    Real l = 24.0;
    _LB = new LandauBrazovskii(.0001,sqrt(l*(l+1)), -.1,.05,.1 ); 
    uint NumMat = _myMesh->getNumberOfElements();
    uint NodeDoF = 1;
    _PbDoF = _myMesh->getNumberOfNodes();
    vector<LandauBrazovskii *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(_LB);
    int TotNumMatProp = NumMat*1;
    _myModel = new LBModel( _myMesh, materials, NodeDoF);
    _myResults = new EigenEllipticResult(_PbDoF, TotNumMatProp);
  }
  int getProbDoF() { return _PbDoF;}

  ~Problem(){
    delete _myModel;
    delete _myMesh;
    delete _LB;
    delete _myResults;
  }
public:
  uint _PbDoF; 
  LoopShellMesh * _myMesh;
  LBModel * _myModel;
  EigenEllipticResult * _myResults;
  LandauBrazovskii* _LB;
};

Problem* prob;
int FV_wrapper(double* indvar, double* args, gsl_vector* out){
  /*
    args = [ --------- PbDof -----------|epsilon ]
  */
  (prob->_myModel)->setField( args ); //sets the dof
  (prob->_LB)->setR(args[prob->_PbDoF] );
  ComputeRequest myRequest = FORCE;
  (prob->_myResults)->setRequest(myRequest);
  (prob->_myModel)->compute( *(prob->_myResults) );
  for(int i=0; i< prob->_PbDoF; i++){
    out->data[i]=(*( (prob->_myResults)->_residual))(i);    
  }
  return 0;
}

int main(int argc, char** argv) {

  prob = new Problem("T5sphere_nodes.dat","T5sphere_conn.dat");
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
  int choice = 1;
  
  FILE* pfile;
  pfile=fopen("T5_sph_harm_4.txt","r");
  double* trivial_sol=new double[Tot_len + 1];
  for(int i=0;i<Tot_len+1;i++) {
    double temp;
    fscanf(pfile,"%lf\n", &temp);
    trivial_sol[i]= temp;
  }
  
  /*
  double* trivial_sol=new double[Tot_len + 1];
  for (int i=0; i<Tot_len; i++) trivial_sol[i]=0.0;
  trivial_sol[Tot_len]=0.1;
  */
  
  gsl_vector* vec =gsl_vector_calloc(Tot_len+1);
  c.F(trivial_sol,vec);
  printf("b=[");
   gsl_vector_fprintf(stdout,vec,"%le");
   printf("];\n");
  
  
  if (choice==1){
    c.setInitialSolution(trivial_sol);
    c.Continue();
  }
  
  if (choice == 2){
    c.posOfPara=Tot_len;  
    
    c.BPFileName="secondary_BP_eps";
    c.DataFileName="outfile.txt";
    c.stepSize= .02;
    c.load_BP("BP_Filename",1);
    c.Continue();
  }
  if (choice == 3){
     c.setInitialSolution(trivial_sol);
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
      int counter = 0;
      for (int i=0; i<eval->size; i++){
	if (eval->data[i] <0) counter++;
      }
      cout << "Number of negative eig vals = " << counter << endl;
      
      //gsl_vector_fprintf(stdout,tgt, "%le");
      //pprintMat(c.jacobian);
      gsl_eigen_symmv_free(work);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
  }
    
}

