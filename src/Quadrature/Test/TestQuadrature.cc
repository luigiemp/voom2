
#include "HexQuadrature.h"
#include "TetQuadrature.h"
#include "TriQuadrature.h"
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

  // -- Tet rule --
  TetQuadrature TetRule_A(1), TetRule_B(2);

  cout << "TetRule_A weight[0] = " << TetRule_A.getQuadWeights()[0] << endl;
  cout << "TetRule_B weight[0] = " << TetRule_B.getQuadWeights()[0] << endl;

  const vector<VectorXd > & TetQuadPoints_A = TetRule_A.getQuadPoints();
  cout << "TetRule_A point[0] = " << TetQuadPoints_A[0](0) << " " << TetQuadPoints_A[0](1) << " " << TetQuadPoints_A[0](2) << endl;
  const vector<VectorXd > & TetQuadPoints_B = TetRule_B.getQuadPoints();
  cout << "TetRule_B point[0] = " << TetQuadPoints_B[0](0) << " " << TetQuadPoints_B[0](1) << " " << TetQuadPoints_B[0](2) << endl;
  cout << endl;

  bool PassTetRuleA = TetRule_A.check(1);
  bool PassTetRuleB = TetRule_B.check(2);

  // -- 2D tri rule --
  TriQuadrature TriRule_A(1), TriRule_B(2);
  bool PassTriRuleA = TriRule_A.check(1);
  bool PassTriRuleB = TriRule_B.check(2);

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
