#include <vector>
#include <iostream>

#include "Potential.h"
#include "APForceVelPotential.h"
#include "BlankPotential.h"
#include "HillForceVelPotential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"

using namespace voom;

int main()
{
	// BlankViscousPotential ViscPotential;
	NewtonianViscousPotential ViscPotential(100, 0);

	double deltaT = 0.01;
	

	Matrix3d F;
	F << 1., 0.0, 0.0, 0.0, 1., 0.0, 0.0, 0.0, 1.;
	// F << 1.0008, .0015, .0100, .0023, 1.0083, .0008, .0091, .0054, 1.0044;
		
	srand(time(NULL));
	for (unsigned int i = 0; i<3; i++) 
		for (unsigned int J = 0; J<3; J++) 
			F(i,J) += 0.1*(double(rand())/RAND_MAX); 
	// cout << "determinant(F) = " << F.determinant() << endl;
	
	cout << "F = " << endl << F << endl;

	double phi = ViscPotential.phi(Matrix3d::Identity(), F, deltaT);	
	cout << "phi = " << phi << endl;

	Matrix3d dphidF = ViscPotential.dphidF(Matrix3d::Identity(), F, deltaT);

	double h = 1E-7;
	
	double first_derivative_error = 0.0;

	Matrix3d dphidFint;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			F(i,j) = F(i,j) + h;
			double phiplus = ViscPotential.phi(Matrix3d::Identity(), F, deltaT);
			
			F(i,j) = F(i,j) - 2 * h;
			double phiminus = ViscPotential.phi(Matrix3d::Identity(), F, deltaT);
			
			F(i,j) = F(i,j) + h;
			dphidFint(i,j) = (phiplus - phiminus)/(2*h);

			first_derivative_error += square(dphidFint(i,j) - dphidF(i,j));
		}
	}

	cout << "dphidF = " << endl << dphidF << endl;
	cout << "dphidFint = " << endl << dphidFint << endl;

	cout << "First derivative error: " << first_derivative_error << endl;


	// Compute Numerical Derivative
	FourthOrderTensor d2phidF2 = ViscPotential.d2phidF2(Matrix3d::Identity(), F, deltaT);

	FourthOrderTensor d2phidF2int(3,3,3,3);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					F(k,l) = F(k,l) + h;
					Matrix3d dphidFplus = ViscPotential.dphidF(Matrix3d::Identity(), F, deltaT);
					F(k,l) = F(k,l) - 2 * h;
					Matrix3d dphidFminus = ViscPotential.dphidF(Matrix3d::Identity(), F, deltaT);
					F(k,l) = F(k,l) + h;

					d2phidF2int.set(i,j,k,l, (dphidFplus(i,j) - dphidFminus(i,j))/(2*h));
				}
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					cout << d2phidF2.get(i,j,k,l) - d2phidF2int.get(i,j,k,l) << endl;
				}
			}
		}
	}

  return 0;
}
