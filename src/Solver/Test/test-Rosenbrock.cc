#include <iostream>
#include "ConjugateGradientWSK.h"
#include "Model.h"
#include "Rosenbrock.h"

//using namespace voom;

int main()
{
  voom::Rosenbrock rosie( 4, true );

  std::vector<voom::Body*> bodies;
  bodies.push_back(&rosie);

  voom::Model mo( bodies );

  voom::ConjugateGradientWSK cg(false);

  cg.setModel( &mo );

  cg.solve();

  rosie.printState();

  return 0;
}
