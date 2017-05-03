#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "GaugeLipid.h"
#include "GaugeModel.h"
#include "continuation/continuation.h"
#include <gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#define M 18 //238
float Z;
using namespace voom;
class Problem{
public:
  Problem(){;}
  Problem(const string Nodes, const string connTable){
    _myMesh = new LoopShellMesh(Nodes, connTable);
    _Mat = new GaugeLipid(2.0,0.0,0.0);
    _Mat->setPressure(-6.0);
    _Mat->setGamma(0.5);
    uint NumMat = _myMesh->getNumberOfElements();
    uint NodeDoF = 3;
    _PbDoF = 3*_myMesh->getNumberOfNodes(); //+1(A)+3(CM)+3(Gauge),for lag.multiplier -1 for rotation constraint
    vector<GaugeLipid *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(_Mat);
    int TotNumMatProp = NumMat*1;
    _myModel = new GaugeModel( _myMesh, materials, NodeDoF);
    _myResults = new EigenEllipticResult(_PbDoF, TotNumMatProp);
  }
  int getProbDoF() { return _PbDoF;}

  ~Problem(){
    delete _myModel;
    delete _myMesh;
    delete _Mat;
    delete _myResults;
  }
public:
  uint _PbDoF; 
  LoopShellMesh * _myMesh;
  GaugeModel * _myModel;
  EigenEllipticResult * _myResults;
  GaugeLipid* _Mat;
};

Problem* prob;
int FV_wrapper(double* indvar, double* args, gsl_vector* out){
  /*
    args = [ --------- PbDof -----------|epsilon ]
  */
  int PbDoF = prob->getProbDoF();
  int j =0;
  for (int cnt=0; cnt < PbDoF; cnt++){
    (prob->_myModel)->setField(cnt, args[j] ); //sets the dof
    j ++;
  }
  (prob->_Mat)->setPressure(args[PbDoF] );
  ComputeRequest myRequest = FORCE;
  (prob->_myResults)->setRequest(myRequest);
  (prob->_myModel)->compute( *(prob->_myResults) );
  int idx = 0;
  for(int cnt=0; cnt< PbDoF; cnt++){
      out->data[idx]= (*( (prob->_myResults)->_residual))(cnt);
      idx ++;
  }
  return 0;
}

int main(int argc, char** argv) {

  const char* file_name = "trivial_nodes.dat" ; 
  prob = new Problem( file_name ,"sphere_212_conn.dat");
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
  

  FILE* pfile;
  pfile=fopen(file_name,"r");
  double x,y,z;
  double temp;
  // two ints at the beginning of the file are "header" info
  fscanf(pfile,"%d", &temp);
  fscanf(pfile,"%d", &temp);
  cout << "Tot_len "<< Tot_len << endl;
  double* trivial_sol=new double[Tot_len + 1];
  int i=0;
  for(int nd=0;nd<(Tot_len)/3;nd++) {
    fscanf(pfile,"%lf ", &x);
    fscanf(pfile,"%lf ", &y);
    fscanf(pfile,"%lf ", &z);
    
      trivial_sol[i]= x;
      trivial_sol[i+1]= y;
      trivial_sol[i+2]= z;
      i=i+3;
  }
  double pressure = -6.99;
  trivial_sol[Tot_len] = pressure;
  
  // for(int i=0;i<Tot_len;i++) {
  //   fscanf(pfile,"%lf ", &x);
  //   trivial_sol[i]= x;///sqrt(x*x+y*y+z*z);
  // }
  
   gsl_vector* vec =gsl_vector_calloc(Tot_len+1);
   c.F(trivial_sol,vec);
   printf("b=[");
   gsl_vector_fprintf(stdout,vec,"%le");
   printf("];");
  
  int choice = 1;
  if (choice==1){
    c.setInitialSolution(trivial_sol);
    c.Continue();
  }
  
  else if(choice==2){
    c.BPFileName="secondary_BP_eps";
    c.DataFileName="branch3.br";
    c.stepSize= 0.1;
    c.load_BP("BP_Filename",3);
    c.Continue();
  }

  
  c.setInitialSolution(trivial_sol);


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
      cout << "size = " << V->size1-1  << endl;
      gsl_vector_fprintf(stdout,tgt, "%le");
      pprintMat(c.jacobian);
      gsl_vector_free(work);
      gsl_vector_free(S);
      gsl_matrix_free(V);

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
      gsl_eigen_symmv_sort(eval,evec, GSL_EIGEN_SORT_ABS_DESC);
      gsl_vector_fprintf(stdout,eval,"%le");	
      cout << endl;
      //gsl_vector_fprintf(stdout,tgt, "%le");
      //pprintMat(c.jacobian);
      gsl_eigen_symmv_free(work);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      */
  }

