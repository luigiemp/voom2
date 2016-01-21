#include <vector>
#include <iostream>

#include "CompNeoHookean.h"
#include "PassMyoA.h"
#include "Jacobian.h"
#include "IsotropicDiffusion.h"
#include "SCElastic.h"
using namespace voom;

int main()
{
  cout << endl << "Testing mechanics material class ... " << endl;
  cout << ".................................... " << endl << endl;
  
  {
    cout << endl << "Testing CompNeoHookean material. " << endl;
    CompNeoHookean MatMech(0, 1.0, 3.0);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = 15;
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += 0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
    

    MatMech.compute(Rm, F);
    
    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J ) << endl;
    //   }
    // }
    
    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;
 
    MatMech.checkConsistency(Rm,F);
  }



  {
    cout << endl << "Testing PassMyoA material. " << endl;
    // PassMyoA MatMech(1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX);
    PassMyoA MatMech(0, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = 15;
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += pow(-1.0, double(i+J) )*0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();

    Real I4 = N.dot((F.transpose()*F)*N);
    cout << "I4 = " << I4 << endl;

    MatMech.compute(Rm, F, &N);

    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int J = 0; J<3; J++) {
  	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J )  << endl;
      }
    }

    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;

    MatMech.checkConsistency(Rm, F, &N);
  }



  {
    cout << endl << "Testing Jacobian material. " << endl;
    Jacobian MatMech(0);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = 7;
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += 0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
    

    MatMech.compute(Rm, F);
    
    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J ) << endl;
    //   }
    // }
    
    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;
 
    MatMech.checkConsistency(Rm,F);
  }




  cout << endl << "Testing diffusion material class ... " << endl;
  cout << ".................................... " << endl << endl;

  {
    Real k = double(rand())/RAND_MAX;
    IsotropicDiffusion MatDiff(k);
    
    DiffusionMaterial::DiffusionResults Rd;
    MatDiff.compute(Rd);
    cout << "Conductivity = " << k << endl;
    cout << "A = " << Rd.A << endl;
    
  }

  {
    cout << endl << "Testing curvature elasticity material class ... " << endl;
    cout << ".................................... " << endl << endl;
    Real t = M_PI/2; Real p = M_PI/4;
    vector<Vector3d> a(2,Vector3d::Zero());
    a[0] << cos(p)*cos(t), sin(p)*cos(t), -sin(t);  //a_theta
    a[1] << -sin(p)*sin(t), cos(p)*sin(t), 0;       //a_phi
    
    vector<Vector3d> aPartials(3,Vector3d::Zero()); //<a_tt, a_pp, a_tp
    aPartials[0] << -cos(p)*sin(t), -sin(p)*sin(t), -cos(t); // a_tt
    aPartials[1] << -cos(p)*sin(t), -sin(p)*sin(t), 0;      // a_pp
    aPartials[2] << -sin(p)*cos(t), cos(p)*cos(t), 0;        // a_tp
    
    ShellGeometry sphere(a, aPartials);

    SCElastic membrane(1,0,0);
    ShellMaterial::Shellresults res; res.request = ENERGY;
    membrane.compute(res, sphere);
    cout << res.W <<endl;
  }
    cout << endl << "....................................... " << endl;
    cout << "Test of voom material classes completed " << endl;
  
  
  return 0;
}
