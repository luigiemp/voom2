/*!
  \file VoomMath.h
  \brief Some basic math functions that are handy to have around
 */
#ifndef _VOOM_MATH_
#define _VOOM_MATH_

#include "voom.h"
#include "ThirdOrderTensor.h"
#include "FourthOrderTensor.h"

namespace voom
{
  inline Real square(Real a) {return (a*a);};
  unsigned int factorial(unsigned int n);

  Real det(const MatrixXd & A);

  void inv(const MatrixXd & A, MatrixXd & B);
};

#endif // _VOOM_MATH_
