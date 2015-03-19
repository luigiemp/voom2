#include "voom.h"
#include <Eigen/Sparse>
// #include <QImage>


typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);

// void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename);

int main(int argc, char** argv)
{
  {
    int n = 300;  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)
    // Assembly:
    std::vector<T> coefficients;            // list of non-zeros coefficients
    Eigen::VectorXd b(m);                   // the right hand side-vector resulting from the constraints
    buildProblem(coefficients, b, n);
    SpMat A(m,m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    // Solving:
    Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
    Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

    // for (uint i = 0; i < x.size(); i++)
    //   cout << x[i] << endl;
    // Export the result to a file:
    // saveAsBitmap(x, n, argv[1]);
  }
  
  {
    int m = 5;
    int n = 5;
    SpMat A(m,m);
    VectorXd x(5), b(5);
    for (int i = 0; i < 5; i++) {
      A.coeffRef(i, i) += 2.0;
      A.coeffRef(i, 1) += 3.0;
      x(i) = 1.0;
    }
    x(0) = -1.0;
    b = A*x;
    cout << "b is = " << endl;
    for (int i = 0; i < 5; i++) {
      cout << b(i) << endl;
    }
  
    // VectorXd r(5); r = A.row(1);
    cout << (A.middleRows(0,1))*x << endl;
  }
  return 0;
}



// Auxiliary functions
void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs,
                       Eigen::VectorXd& b, const Eigen::VectorXd& boundary)
{
  int n = boundary.size();
  int id1 = i+j*n;
        if(i==-1 || i==n) b(id) -= w * boundary(j); // constrained coefficient
  else  if(j==-1 || j==n) b(id) -= w * boundary(i); // constrained coefficient
  else  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n)
{
  b.setZero();
  Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);
  for(int j=0; j<n; ++j)
  {
    for(int i=0; i<n; ++i)
    {
      int id = i+j*n;
      insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j,    4, coefficients, b, boundary);
    }
  }
}

// void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename)
// {
//   Eigen::Array<unsigned char,Eigen::Dynamic,Eigen::Dynamic> bits = (x*255).cast<unsigned char>();
//   QImage img(bits.data(), n,n,QImage::Format_Indexed8);
//   img.setColorCount(256);
//   for(int i=0;i<256;i++) img.setColor(i,qRgb(i,i,i));
//   img.save(filename);
// }

