#include "NewtonianViscousPotential.h"

namespace voom
{
	double NewtonianViscousPotential::phi(Matrix3d Fn, Matrix3d Fnp1, double deltaT)
	{
		Matrix3d Fnpalpha = (1. - _alpha) * Fn + _alpha * Fnp1;
		Matrix3d invFnpalpha = Fnpalpha.inverse();
		double Jnpalpha = Fnpalpha.determinant();

		Matrix3d dsymm = 1./deltaT * ((Fnp1 - Fn) * invFnpalpha);
		Matrix3d d = 0.5 * (dsymm + dsymm.transpose());

		Matrix3d ddev = d - 1./3. * d.trace() * Matrix3d::Identity();

		double phi = _eta * Jnpalpha * (ddev.transpose() * ddev).trace();
		return phi;
	}

	Matrix3d NewtonianViscousPotential::dphidF(Matrix3d Fn, Matrix3d Fnp1, double deltaT)
	{
		Matrix3d IdentityMatrix = Matrix3d::Identity();
		
		Matrix3d Fnpalpha = (1. - _alpha) * Fn + _alpha * Fnp1;
		double detFnpalpha = Fnpalpha.determinant();
		Matrix3d invFnpalpha = Fnpalpha.inverse();
		Matrix3d invTFnpalpha = invFnpalpha.transpose();
		double Jnpalpha = Fnpalpha.determinant();

		Matrix3d dsymm = 1./deltaT * ((Fnp1 - Fn) * invFnpalpha);
		Matrix3d d = 0.5 * (dsymm + dsymm.transpose());

		Matrix3d ddev = d - 1./3. * d.trace() * Matrix3d::Identity();

		Matrix3d dJdF = Matrix3d::Zero();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			    dJdF(i,j) = _alpha * detFnpalpha * invFnpalpha(j,i);

		Matrix3d Term1 = dJdF * (ddev.transpose() * ddev).trace();


		// Term 2:
		FourthOrderTensor ddnaFnp1(3,3,3,3);
		FourthOrderTensor Temp1(3,3,3,3);
		FourthOrderTensor Temp2(3,3,3,3);
		FourthOrderTensor Temp3(3,3,3,3);
		FourthOrderTensor Temp4(3,3,3,3);
		FourthOrderTensor Temp5(3,3,3,3);
		FourthOrderTensor Temp6(3,3,3,3);

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 3; k++) {
							Temp2.set(m, n, i, j, Temp2.get(m, n, i, j) - _alpha * Fnp1(m,k) * invFnpalpha(k,i) * invFnpalpha(j,n));
							Temp3.set(m, n, i, j, Temp3.get(m, n, i, j) + _alpha * Fn(m,k) * invFnpalpha(k,i) * invFnpalpha(j,n));
							Temp4.set(m, n, i, j, Temp4.get(m, n, i, j) - _alpha * invFnpalpha(k,i) * invFnpalpha(j,m) * Fnp1(n,k));
							Temp6.set(m, n, i, j, Temp6.get(m, n, i, j) + _alpha * Fn(n,k) * invFnpalpha(k,i) * invFnpalpha(j,m));
						}
						Temp5.set(m, n, i, j, invFnpalpha(j,m) * IdentityMatrix(n,i));
						Temp1.set(m, n, i, j, IdentityMatrix(m,i) * invFnpalpha(j,n));
						ddnaFnp1.set(m, n, i, j, 1./(2.*deltaT) * (Temp1.get(m, n, i, j) + Temp2.get(m, n, i, j) + Temp3.get(m, n, i, j) + Temp4.get(m, n, i, j) + Temp5.get(m, n, i, j) + Temp6.get(m, n, i, j)));
					}
				}
			}
		}

		// Term 2.2:
		Matrix3d temp22_1 = Matrix3d::Zero();
		Matrix3d temp22_2 = Matrix3d::Zero();
		Matrix3d temp22_3 = Matrix3d::Zero();
		Matrix3d tempterm22 = Matrix3d::Zero();
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				temp22_1(i,j) = invFnpalpha(j,i);
				for (int g = 0; g < 3; g++) {
					for (int k = 0; k < 3; k++) {
						temp22_2(i,j) = temp22_2(i,j) - _alpha * Fnp1(g,k) * invFnpalpha(k,i) * invFnpalpha(j,g);
						temp22_3(i,j) = temp22_3(i,j) + _alpha * Fn(g,k) * invFnpalpha(k,i) * invFnpalpha(j,g);
					}
				}
			}
		}
		tempterm22 = temp22_1 + temp22_2 + temp22_3;

		FourthOrderTensor term22(3,3,3,3);
		FourthOrderTensor ddevFnp1(3,3,3,3);
		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						term22.set(m, n, i, j, -1./(3. * deltaT) * tempterm22(i,j) * IdentityMatrix(m,n));
						ddevFnp1.set(m,n,i,j, ddnaFnp1.get(m,n,i,j) + term22.get(m,n,i,j));
					}
				}
			}
		}

		Matrix3d Term2 = Matrix3d::Zero();
		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						Term2(m,n) = Term2(m,n) + ddev(i,j) * ddevFnp1.get(i,j,m,n);
					}
				}
			}
		}
		Term2 = Term2 * 2. * Jnpalpha;

		Matrix3d dphidF = _eta * (Term1 + Term2);

		return dphidF;
	}

	FourthOrderTensor NewtonianViscousPotential::d2phidF2(Matrix3d Fn, Matrix3d Fnp1, double deltaT)
	{
		FourthOrderTensor d2phidF2(3,3,3,3);

		// Compute Numerical Derivative
		/*
		double h = 1.0E-7;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 3; l++) {
						Fnp1(k,l) = Fnp1(k,l) + h;
						Matrix3d dphidFplus = dphidF(Fn, Fnp1, deltaT);
						Fnp1(k,l) = Fnp1(k,l) - 2. * h;
						Matrix3d dphidFminus = dphidF(Fn, Fnp1, deltaT);
						Fnp1(k,l) = Fnp1(k,l) + h;

						d2phidF2.set(i,j,k,l, (dphidFplus(i,j) - dphidFminus(i,j))/(2.*h));
					}
				}
			}
		}
		*/

		// Compute Analytical Second Derivative
		Matrix3d IdentityMatrix = Matrix3d::Identity();
		
		Matrix3d Fnpalpha = (1. - _alpha) * Fn + _alpha * Fnp1;
		double detFnpalpha = Fnpalpha.determinant();
		Matrix3d invFnpalpha = Fnpalpha.inverse();
		Matrix3d invTFnpalpha = invFnpalpha.transpose();
		double Jnpalpha = Fnpalpha.determinant();

		Matrix3d dsymm = 1./deltaT * ((Fnp1 - Fn) * invFnpalpha);
		Matrix3d d = 0.5 * (dsymm + dsymm.transpose());

		Matrix3d ddev = d - 1./3. * d.trace() * Matrix3d::Identity();

		Matrix3d dJdF = Matrix3d::Zero();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			    dJdF(i,j) = _alpha * detFnpalpha * invFnpalpha(j,i);


		FourthOrderTensor dFinvnpalphadFnpalpha(3,3,3,3);
		for (int J = 0; J < 3; J++)
			for (int i = 0; i < 3; i++)
				for (int o = 0; o < 3; o++)
					for (int P = 0; P < 3; P++)
						dFinvnpalphadFnpalpha.set(J,i,o,P, -1 * _alpha * invFnpalpha(J,o) * invFnpalpha(P,i));

		FourthOrderTensor d2JdF2(3,3,3,3);
		for (int i = 0; i < 3; i++)
			for (int J = 0; J < 3; J++)
				for (int m = 0; m < 3; m++)
					for (int N = 0; N < 3; N++)
						d2JdF2.set(i,J,m,N, _alpha*_alpha * (Jnpalpha * invFnpalpha(N,m) * invFnpalpha(J,i) - Jnpalpha * invFnpalpha(J,m) * invFnpalpha(N,i)));


		// Begin: This section of the code has been copied from the stress calculation
		Matrix3d Term1 = dJdF * (ddev.transpose() * ddev).trace();


		// Term 2:
		FourthOrderTensor ddnaFnp1(3,3,3,3);
		FourthOrderTensor Temp1(3,3,3,3);
		FourthOrderTensor Temp2(3,3,3,3);
		FourthOrderTensor Temp3(3,3,3,3);
		FourthOrderTensor Temp4(3,3,3,3);
		FourthOrderTensor Temp5(3,3,3,3);
		FourthOrderTensor Temp6(3,3,3,3);

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 3; k++) {
							Temp2.set(m, n, i, j, Temp2.get(m, n, i, j) - _alpha * Fnp1(m,k) * invFnpalpha(k,i) * invFnpalpha(j,n));
							Temp3.set(m, n, i, j, Temp3.get(m, n, i, j) + _alpha * Fn(m,k) * invFnpalpha(k,i) * invFnpalpha(j,n));
							Temp4.set(m, n, i, j, Temp4.get(m, n, i, j) - _alpha * invFnpalpha(k,i) * invFnpalpha(j,m) * Fnp1(n,k));
							Temp6.set(m, n, i, j, Temp6.get(m, n, i, j) + _alpha * Fn(n,k) * invFnpalpha(k,i) * invFnpalpha(j,m));
						}
						Temp5.set(m, n, i, j, invFnpalpha(j,m) * IdentityMatrix(n,i));
						Temp1.set(m, n, i, j, IdentityMatrix(m,i) * invFnpalpha(j,n));
						ddnaFnp1.set(m, n, i, j, 1./(2.*deltaT) * (Temp1.get(m, n, i, j) + Temp2.get(m, n, i, j) + Temp3.get(m, n, i, j) + Temp4.get(m, n, i, j) + Temp5.get(m, n, i, j) + Temp6.get(m, n, i, j)));
					}
				}
			}
		}

		// Term 2.2:
		Matrix3d temp22_1 = Matrix3d::Zero();
		Matrix3d temp22_2 = Matrix3d::Zero();
		Matrix3d temp22_3 = Matrix3d::Zero();
		Matrix3d tempterm22 = Matrix3d::Zero();
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				temp22_1(i,j) = invFnpalpha(j,i);
				for (int g = 0; g < 3; g++) {
					for (int k = 0; k < 3; k++) {
						temp22_2(i,j) = temp22_2(i,j) - _alpha * Fnp1(g,k) * invFnpalpha(k,i) * invFnpalpha(j,g);
						temp22_3(i,j) = temp22_3(i,j) + _alpha * Fn(g,k) * invFnpalpha(k,i) * invFnpalpha(j,g);
					}
				}
			}
		}
		tempterm22 = temp22_1 + temp22_2 + temp22_3;

		FourthOrderTensor term22(3,3,3,3);
		FourthOrderTensor ddevFnp1(3,3,3,3);
		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						term22.set(m, n, i, j, -1./(3. * deltaT) * tempterm22(i,j) * IdentityMatrix(m,n));
						ddevFnp1.set(m,n,i,j, ddnaFnp1.get(m,n,i,j) + term22.get(m,n,i,j));
					}
				}
			}
		}
		// End: This section of the code has been copied from the stress calculation
		FourthOrderTensor term1(3,3,3,3);
		FourthOrderTensor term2(3,3,3,3);
		FourthOrderTensor d2trdnpalphad2Fnp1(3,3,3,3);

		for (int i = 0; i < 3; i++) {
			for (int J = 0; J < 3; J++) {
				for (int m = 0; m < 3; m++) {
					for (int N = 0; N < 3; N++) {
						term1.set(i,J,m,N, -2 * _alpha * invFnpalpha(J,m) * invFnpalpha(N,i));
						for (int k = 0; k < 3; k++) {
							for (int p = 0; p < 3; p++) {
								term2.set(i,J,m,N, term2.get(i,J,m,N) + (Fnp1(k,p) - Fn(k,p)) * _alpha * _alpha * (invFnpalpha(p,m) * invFnpalpha(N,i) * invFnpalpha(J,k) + invFnpalpha(p,i) * invFnpalpha(J,m) * invFnpalpha(N,k)));
							}
						}
						d2trdnpalphad2Fnp1.set(i,J,m,N, 1/deltaT * (term1.get(i,J,m,N) + term2.get(i,J,m,N)));
					}
				}
			}
		}
		
		FourthOrderTensor temp1(3,3,3,3);
		FourthOrderTensor temp2(3,3,3,3);
		FourthOrderTensor temp3(3,3,3,3);
		FourthOrderTensor temp4(3,3,3,3);
		FourthOrderTensor temp5(3,3,3,3);


		for (int i = 0; i < 3; i++) {
			for (int J = 0; J < 3; J++) {
				for (int m = 0; m < 3; m++) {
					for (int N = 0; N < 3; N++) {
						for (int k = 0; k < 3; k++) {
							for (int L = 0; L < 3; L++) {
								temp1.set(i,J,m,N, temp1.get(i,J,m,N) + d2JdF2.get(i,J,m,N) * ddev(k,L) * ddev(k,L));
								temp2.set(i,J,m,N, temp2.get(i,J,m,N) + 2 * dJdF(i,J) * ddevFnp1.get(k,L,m,N) * ddev(k,L));
								temp3.set(i,J,m,N, temp3.get(i,J,m,N) + 2 * dJdF(m,N) * ddev(k,L) * ddevFnp1.get(k,L,i,J));
								temp4.set(i,J,m,N, temp4.get(i,J,m,N) + 2 * Jnpalpha * ddevFnp1.get(k,L,m,N) * ddevFnp1.get(k,L,i,J));
								temp5.set(i,J,m,N, temp5.get(i,J,m,N) + 1./deltaT * Jnpalpha * ddev(k,L) * (-_alpha * IdentityMatrix(L,m) * invFnpalpha(J,k) * invFnpalpha(N,i) - _alpha * IdentityMatrix(k,m) * invFnpalpha(J,L) * invFnpalpha(N,i) - _alpha * IdentityMatrix(L,i) * invFnpalpha(J,m) * invFnpalpha(N,k) - _alpha * IdentityMatrix(k,i) * invFnpalpha(J,m) * invFnpalpha(N,L)) - 2./3. * Jnpalpha * ddev(k,L) * d2trdnpalphad2Fnp1.get(i,J,m,N) * IdentityMatrix(k,L));
								for (int p = 0; p < 3; p++) {
									temp5.set(i,J,m,N, temp5.get(i,J,m,N) + 1./deltaT * Jnpalpha * ddev(k,L) * (-1. * _alpha * _alpha * Fn(L,p) * invFnpalpha(J,m) * invFnpalpha(N,k) * invFnpalpha(p,i) + _alpha * _alpha * Fnp1(L,p) * invFnpalpha(J,m) * invFnpalpha(N,k) * invFnpalpha(p,i) - _alpha * _alpha * Fn(k,p) * invFnpalpha(J,m) * invFnpalpha(N,L) * invFnpalpha(p,i) + _alpha * _alpha * Fnp1(k,p) * invFnpalpha(J,m) * invFnpalpha(N,L) * invFnpalpha(p,i) - _alpha * _alpha * Fn(L,p) * invFnpalpha(J,k) * invFnpalpha(N,i) * invFnpalpha(p,m) + _alpha * _alpha * Fnp1(L,p) * invFnpalpha(J,k) * invFnpalpha(N,i) * invFnpalpha(p,m) - _alpha * _alpha * Fn(k,p) * invFnpalpha(J,L) * invFnpalpha(N,i) * invFnpalpha(p,m) + _alpha * _alpha * Fnp1(k,p) * invFnpalpha(J,L) * invFnpalpha(N,i) * invFnpalpha(p,m)));
								}
							}
						}
						d2phidF2.set(i,J,m,N, _eta* (temp1.get(i,J,m,N) + temp2.get(i,J,m,N) + temp3.get(i,J,m,N) + temp4.get(i,J,m,N) + temp5.get(i,J,m,N)));
					}
				}
			}
		}



		// Return Statement
		return d2phidF2;
	}
}
