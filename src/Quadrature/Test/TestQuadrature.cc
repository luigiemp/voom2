
#include "HexQuadrature.h"
#include "LineQuadrature.h"
#include "QuadQuadrature.h"

using namespace voom;

int main()
{
  cout << endl << "Testing voom quadrature classes ... " << endl;
  cout << "...................................." << endl << endl;
  
  // -- Hex quad rule --
  HexQuadrature QuadRule_A(1), QuadRule_B(3);

  cout << "QuadRule_A weight[0] = " << QuadRule_A.getQuadWeights()[0] << endl;
  cout << "QuadRule_B weight[6] = " << QuadRule_B.getQuadWeights()[6] << endl;

  const vector<VectorXd > & QuandPoints_A = QuadRule_A.getQuadPoints();
  cout << "QuadRule_A point[0] = " << QuandPoints_A[0](0) << " " << QuandPoints_A[0](1) << " " << QuandPoints_A[0](2) << endl;
  const vector<VectorXd > & QuandPoints_B = QuadRule_B.getQuadPoints();
  cout << "QuadRule_B point[6] = " << QuandPoints_B[6](0) << " " << QuandPoints_B[6](1) << " " << QuandPoints_B[6](2) << endl;
  cout << endl;

  bool PassA = QuadRule_A.check(1);
  bool PassB = QuadRule_B.check(3);

  // -- Lin quad rule --
  LineQuadrature LinQuadRule(3);
  bool PassLineQuad = LinQuadRule.check(5);

  // --Quad quad rule
  QuadQuadrature QuadQuadRule(10);
  bool PassQuadQuad = QuadQuadRule.check(19);
  
  cout << endl << "................................ " << endl;
  cout << "Test of voom quadrature classes completed" << endl;
  return 0;
}
