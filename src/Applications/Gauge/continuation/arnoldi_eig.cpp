#include "ardsmat.h"
#include "ardgsym.h"
#include "lsymsol.h"
#include"arnoldi_eig.h"

int arnoldi_eig(gsl_matrix* a, gsl_matrix* mass, int no_of_eigs, gsl_vector* eigs, gsl_matrix* evecs, double sigma)
{

  // Defining variables;

  int     n=a->size1;       // Dimension of the problem.
  double* valM;
  double* valA;

  // Creating matrices A and B.
  arpack_matrix_wrapper( valA, valM, a, mass);
  ARdsSymMatrix<double> A(n, valA);
  ARdsSymMatrix<double> M(n, valM);

 // Defining what we need: the four eigenvectors nearest to 0.0.

 ARluSymGenEig<double> dprob('S', no_of_eigs, A,M, sigma);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();
  for(int e=0; e<no_of_eigs; e++){
    // exporting eigenvalues
    eigs->data[e]=dprob.Eigenvalue(e);
    //exporting eigenvectors
    for(int v=0;v<n; v++)
   gsl_matrix_set(evecs,v,e,dprob.RawEigenvector(e)[v]);
  }
  
}
template<class ARFLOAT>
void arpack_matrix_wrapper( ARFLOAT* &A, ARFLOAT* &B, gsl_matrix* a, gsl_matrix* b)
{
  //check for size compatibility of a and b

  int n=a->size1;
  
  int index=0;
  A  = new ARFLOAT[(n*n+n)/2];
  B  = new ARFLOAT[(n*n+n)/2];
  // column wise lower triangular is the default for the package
  for(int col=0; col < n;col++){
    for(int row=col;row<n;row++){
      
      A[index] = gsl_matrix_get(a,row,col);
      B[index] =  gsl_matrix_get(b,row,col);
      //printf("%lf\n",A[index]);
      index++;
    }
  }
}
