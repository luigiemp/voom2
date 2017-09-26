#include <vector>
#include <iostream>

#include "CompNeoHookean.h"
#include "PassMyoA.h"
#include "Jacobian.h"
#include "Holzapfel.h"
#include "Guccione.h"
#include "IsotropicDiffusion.h"
#include "SCElastic.h"
#include "Humphrey.h"
#include "Humphrey_Compressible.h"
#include "LinYinActive.h"
#include "LinYinActive_Compressible.h"
#include "ImposedKinematics.h"
using namespace voom;

int main()
{
  cout << endl << "Testing mechanics material class ... " << endl;
  
  {
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing CompNeoHookean material. " << endl;
    CompNeoHookean MatMech(0, 1.0, 3.0);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
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
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing PassMyoA material. " << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();
    vector<Vector3d > Fibers;
    Fibers.push_back(N);
    PassMyoA MatMech(0, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, Fibers);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += pow(-1.0, double(i+J) )*0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
   
    Real I4 = N.dot((F.transpose()*F)*N);
    cout << "I4 = " << I4 << endl;

    MatMech.compute(Rm, F);

    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J )  << endl;
    //   }
    // }

    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;

    MatMech.checkConsistency(Rm, F);
  }



  {
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing Jacobian material. " << endl;
    Jacobian MatMech(0);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS);
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
    cout << ".................................... " << endl << endl;
    cout << endl << "Testing Holzapfel material. " << endl;
    Vector3d f, s;
    f << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    s << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    f /= f.norm();
    s /= s.norm();
    vector<Vector3d > Fibers;
    Fibers.push_back(f);
    Fibers.push_back(s);
    Holzapfel MatMech(0, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 
		      1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, Fibers);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | DMATPROP);
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += pow(-1.0, double(i+J) )*0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
   

    MatMech.compute(Rm, F);

    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J )  << endl;
    //   }
    // }

    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;

    MatMech.checkConsistency(Rm, F);
  }



 {
   cout << ".................................... " << endl << endl;
    cout << endl << "Testing Guccione material. " << endl;
    Vector3d f, c, r;
    f << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    c << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    r << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    f /= f.norm();
    c /= c.norm();
    r /= r.norm();
    vector<Vector3d > Fibers;
    Fibers.push_back(f);
    Fibers.push_back(c);
    Fibers.push_back(r);
    Guccione MatMech(0, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, Fibers);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += pow(-1.0, double(i+J) )*0.1*(double(rand())/RAND_MAX);
    cout << "determinant(F) = " << F.determinant() << endl;
   

    MatMech.compute(Rm, F);

    cout << "Energy     = " << Rm.W << endl;
    cout << "P(2,2)     = " << Rm.P(2,2) << endl;
    cout << "K[0,0,0,0] = " << Rm.K.get(0,0,0,0) << endl;
    // for (unsigned int i = 0; i<3; i++) {
    //   for (unsigned int J = 0; J<3; J++) {
    // 	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << " " << (Rm.Dmat).get( 1, i, J )  << endl;
    //   }
    // }

    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;

    MatMech.checkConsistency(Rm, F);
  }

 {
   cout << ".................................... " << endl << endl;
   cout << endl << "Testing diffusion material class ... " << endl;
 
    Real k = double(rand())/RAND_MAX;
    IsotropicDiffusion MatDiff(k);
    
    DiffusionMaterial::DiffusionResults Rd;
    MatDiff.compute(Rd);
    cout << "Conductivity = " << k << endl;
    cout << "A = " << Rd.A << endl;
  }

  {
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing curvature elasticity material class ... " << endl;
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

  {
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing Humphrey material. " << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();

    vector<Vector3d> DirVecs;
    DirVecs.push_back(N);
    DirVecs.push_back(N);
    DirVecs.push_back(N);

    Humphrey MatMech(0, 1.0, 3.0, 5.0, 7.0, 9.0, DirVecs);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
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

    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int J = 0; J<3; J++) {
	cout << i << " " << J << " " << (Rm.Dmat).get( 0, i, J ) << endl;
      }
    }
    
    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;
 
    MatMech.checkConsistency(Rm,F);
  }

  {
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing Lin Yin Active material. " << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();

    vector<Vector3d> DirVecs;
    DirVecs.push_back(N);
    DirVecs.push_back(N);
    DirVecs.push_back(N);

    LinYinActive MatMech(0, 0.0, -13.03, 36.65, 35.42, 15.52, 1.62, DirVecs);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS);
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += 0.1*(double(rand())/RAND_MAX);
    // F << 1.06331, 0.0511631, .00747381, 0.081262, 1.0603, 0.0943724, 0.0544078, 0.0547286, 1.05075;
    cout << "determinant(F) = " << F.determinant() << endl;
    cout << "F = \n " << F << endl;
    
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
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing Lin Yin Incompressible Active material. " << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();

    vector<Vector3d> DirVecs;
    DirVecs.push_back(N);
    DirVecs.push_back(N);
    DirVecs.push_back(N);

    LinYinActive_Compressible MatMech(0, -13.03, 36.65, 35.42, 0.1, 4.*0.45*0.1/(1. - 2. * 0.45), DirVecs);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS);
    Matrix3d F;
    F << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    srand(time(NULL));
    for (unsigned int i = 0; i<3; i++) 
      for (unsigned int J = 0; J<3; J++) 
	F(i,J) += 0.1*(double(rand())/RAND_MAX);
    // F << 1.06331, 0.0511631, .00747381, 0.081262, 1.0603, 0.0943724, 0.0544078, 0.0547286, 1.05075;
    cout << "determinant(F) = " << F.determinant() << endl;
    cout << "F = \n " << F << endl;
    
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
 
    MatMech.checkConsistency(Rm,F, 1.0e-6);
  }

  {
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing Humphrey Compressible material. " << endl;
    Vector3d N;
    N << (double(rand())/RAND_MAX), double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    N /= N.norm();

    vector<Vector3d> DirVecs;
    DirVecs.push_back(N);
    DirVecs.push_back(N);
    DirVecs.push_back(N);

    Humphrey_Compressible MatMech(0, 1.0, 3.0, 5.0, 7.0, 9.0, 0.1, 4.*0.45*0.1/(1. - 2. * 0.45), DirVecs );
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
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
    cout << endl << ".................................... " << endl << endl;
    cout << endl << "Testing ImposedKinematics material. " << endl;
    
    int NumAlphas = 5;
    int MatID = 0;
    vector<Real > Alphas;
    vector<Real > Stretches;
    vector<Vector3d > Directions;
    Real Beta = double(rand())/RAND_MAX;
    
    for (int a = 0; a < NumAlphas; a++) {
      Alphas.push_back(double(rand())/RAND_MAX);
      Stretches.push_back(double(rand())/RAND_MAX);
      Vector3d N; N << double(rand())/RAND_MAX, double(rand())/RAND_MAX, double(rand())/RAND_MAX;
      N /= N.norm();
      Directions.push_back(N);
    }

    ImposedKinematics MatMech(MatID, Alphas, Stretches, Directions, Beta);
    
    MechanicsMaterial::FKresults Rm;
    Rm.request = (ENERGY | FORCE | STIFFNESS | DMATPROP);
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
    cout << endl << "Material ID = " << MatMech.getMatID() << endl << endl;
 
    MatMech.checkConsistency(Rm,F);
 }




    cout << endl << "....................................... " << endl;
    cout << "Test of voom material classes completed " << endl;
  
  

  return 0;
}
