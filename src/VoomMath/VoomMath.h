/*!
  \file VoomMath.h
  \brief Some basic math functions that are handy to have around
 */
#ifndef _VOOM_MATH_
#define _VOOM_MATH_

#include "voom.h"
#include <math.h>
// WARNING
#include <unsupported/Eigen/MatrixFunctions>
// END WARNING
#include <Eigen/Eigenvalues>
#include "ThirdOrderTensor.h"
#include "FourthOrderTensor.h"

namespace voom
{
  inline Real square(Real a) {return (a*a); };
  uint factorial(uint n);

  Real det(const MatrixXd & A);

  void inv(const MatrixXd & A, MatrixXd & B);

  Matrix3d VoomExpSymmMatrix(const Matrix3d &A);
};

#endif // _VOOM_MATH_
