#include "PlasticMaterial.h"

namespace voom
{
	void PlasticMaterial::optimizeInternalVariables()
	{
		_Qnp1 = _Q;
		
		int iter = 0;
		for (iter = 0; iter < _maxIter; iter++)
		{
			Vector3d dWdQ = computedWdQ();
			Matrix3d d2WdQ2 = computed2WdQ2();

			Vector3d dQ = -1. * d2WdQ2.inverse() * dWdQ;
			_Qnp1 = _Qnp1 + dQ;

			// cout << "****" << endl << dWdQ << endl << "****" << endl << d2WdQ2 << endl;

			double error = sqrt(square(dQ(0)) + square(dQ(1)) + square(dQ(2)));

			if(error < _hardOptTOL)
				break;
		}
		if (iter == _maxIter)
			cout << "****** Maximum Newton iterations reached for hardening variable optimization. ******" << endl;

	}

	Vector3d PlasticMaterial::computedWdQ()
	{
		Vector3d dWdQ(0,0,0);
		// Compute \mathbf{F}^a_{n+1} according to Flow Rule
		Vector3d deltaQ = _Qnp1 - _Q;
		Matrix3d Fanp1 = Matrix3d::Zero();
		Matrix3d expFanp1 = Matrix3d::Zero();

		Matrix3d A = Matrix3d::Zero();
		for (int i = 0; i < 3; i++){
			A += deltaQ(i) * _M[i];
		}

		Fanp1 = A.exp() * _Fa;

		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = Fanp1.inverse();

		Matrix3d Fenp1 = _Fnp1 * invFanp1;

		FKresults ActiveResults;
		ActiveResults.request = 2;

		_ActiveMaterial->compute(ActiveResults, Fenp1, &_dirVec[0]);
		
		vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
		for (int p = 0; p < 3; p++)
			dFadQ[p] = A.exp() * _M[p] * _Fa;
		
		Vector3d dDsdQnp1(0,0,0);
		for (int a = 0; a < 3; a++)
			for (int k = 0; k < 3; k++)
				for (int L = 0; L < 3; L++)
					for (int i = 0; i < 3; i++)
						for (int J = 0; J < 3; J++)
							for (int p = 0; p < 3; p++)
								dDsdQnp1(a) = dDsdQnp1(a) + ActiveResults.P(i,J) * _Fnp1(i,p) * invFanp1(p,k) * invFanp1(L,J) * dFadQ[a](k,L);
		dDsdQnp1 = dDsdQnp1 * -1;

		// Compute \dpsi^{*}dQnp1 from the Force-velocity potential
		Vector3d dpsidQnp1 = _activation * _KineticPotential->DPsiDQ(deltaQ, _deltaT);
		// TODO: Add in deltaT * dpsidQnp1 to dWdQ	// DONE

		dWdQ = dDsdQnp1 + dpsidQnp1 * _deltaT;

		// cout << dDsdQnp1 << endl << endl << dpsidQnp1 << endl; exit(1);
		// cout << ActiveResults.P << endl; exit(1);
		// cout << Fenp1 << endl; exit(1);

		return dWdQ;
	}

	Matrix3d PlasticMaterial::computed2WdQ2()
	{
		Matrix3d d2WdQ2 = Matrix3d::Zero();

		Matrix3d IdentityTensor = Matrix3d::Identity();

		// Compute \mathbf{F}^a_{n+1} according to Flow Rule
		Vector3d deltaQ = _Qnp1 - _Q;
		Matrix3d Fanp1 = Matrix3d::Zero();

		Matrix3d A = Matrix3d::Zero();
		for (int i = 0; i < 3; i++){
			A += deltaQ(i) * _M[i];
		}

		Fanp1 = A.exp() * _Fa;

		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = Fanp1.inverse();

		Matrix3d Fenp1 = _Fnp1 * invFanp1;

		FKresults ActiveResults;
		ActiveResults.request = 6;

		_ActiveMaterial->compute(ActiveResults, Fenp1, &_dirVec[0]);
		
		vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
		for (int p = 0; p < 3; p++)
			dFadQ[p] = A.exp() * _M[p] * _Fa;

		Matrix3d expA = A.exp();

		FourthOrderTensor d2FadQ2(3,3,3,3);
		FourthOrderTensor dFedFanp1(3,3,3,3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Matrix3d temp = expA * _M[i] * _M[j] * _Fa;
				for (int p = 0; p < 3; p++) {
					for (int q = 0; q < 3; q++) {
						d2FadQ2.set(i,j,p,q, temp(p,q));
						for (int r = 0; r < 3; r++) {
							dFedFanp1.set(i,j,p,q, dFedFanp1.get(i,j,p,q) + _Fnp1(i,r) * invFanp1(r,p) * invFanp1(q,j));
						}
						dFedFanp1.set(i,j,p,q, -1. * dFedFanp1.get(i,j,p,q));
					}
				}
			}
		}

		Matrix3d term1 = Matrix3d::Zero();
		Matrix3d term2 = Matrix3d::Zero();
		Matrix3d term3 = Matrix3d::Zero();

		for (int p = 0; p < 3; p++) {
			for (int q = 0; q < 3; q++) {
				for (int e = 0; e < 3; e++) {
					for (int f = 0; f < 3; f++) {
						for (int c = 0; c < 3; c++) {
							for (int d = 0; d < 3; d++) {
								for (int h = 0; h < 3; h++) {
									for (int i = 0; i < 3; i++) {
										for (int l = 0; l < 3; l++) {
											for (int m = 0; m < 3; m++) {
		                                        					term1(p,q) = term1(p,q) + ActiveResults.K.get(e,f,h,i) * dFedFanp1.get(e,f,c,d) * dFedFanp1.get(h,i,l,m) * dFadQ[p](c,d) * dFadQ[q](l,m);
											}
										}
										term2(p,q) = term2(p,q) + ActiveResults.P(e,f) * -1 * (dFedFanp1.get(e,c,h,i) * invFanp1(d,f) + dFedFanp1.get(e,h,c,d) * invFanp1(i,f)) * dFadQ[p](c,d) * dFadQ[q](h,i);
									}
								}
								term3(p,q) = term3(p,q) + ActiveResults.P(e,f) * dFedFanp1.get(e,f,c,d) * d2FadQ2.get(c,d,p,q);
							}
						}
					}
				}
			}
		}
		
		// DELETE:
		// cout << term1 << endl << term2 << endl << term3 << endl;

		Matrix3d d2DdQ2np1 = term1 + term2 + term3;

		// Compute \dpsi^{*}dQnp1 from the Force-velocity potential
		Matrix3d d2psidQnp12 = _activation * _KineticPotential->D2PsiDQDQ(deltaQ, _deltaT);
		// TODO: Add in deltaT * d2psidQnp12 to d2DdQ2np1	// DONE
		
		d2WdQ2 = d2DdQ2np1 + _deltaT * d2psidQnp12;

		return d2WdQ2;
	}
	// TODO: Change the compute function to take vector<Vector3d> for all three directions.
	void PlasticMaterial::compute(FKresults & R, const Matrix3d & Fnp1, Vector3d * fiber)
	{
		//cout << R.request << endl;
		_Fnp1 = Fnp1;
		computeKinematicParameters();

		// Optimize for the hardening variables
		// Find \mathbf{Q}_{n+1}
		optimizeInternalVariables();
		
		// Compute \mathbf{F}^a_{n+1} according to Flow Rule
		Vector3d deltaQ = _Qnp1 - _Q;
		Matrix3d Fanp1 = Matrix3d::Zero();

		Matrix3d A = Matrix3d::Zero();
		for (int i = 0; i < 3; i++){
			A += deltaQ(i) * _M[i];
		}

		Fanp1 = A.exp() * _Fa;

		// TODO: CHANGE ALL Fanp1 to _Fanp1 but for now:
		_Fanp1 = Fanp1;

		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = Fanp1.inverse();

		Matrix3d Fen = _Fn * invFa;
		Matrix3d Fenp1 = Fnp1 * invFanp1;

		FKresults PassiveResults; FKresults ActiveResults;
		PassiveResults.request = R.request; ActiveResults.request = R.request;
		if (R.request & STIFFNESS)
		{
			// STIFFNESS requires computation of the FORCE of Active and Passive materials
			if (!(R.request & FORCE))
				R.request += 2;
		}
		_PassiveMaterial->compute(PassiveResults, Fnp1, &_dirVec[0]);
		_ActiveMaterial->compute(ActiveResults, Fenp1, &_dirVec[0]);
		
		if (R.request & ENERGY)
		{
			// Compute Viscous Potential Stress
			double phi = _ViscousPotential->phi(_Fn, Fnp1, _deltaT);

			// At the Old Timestep (Only needed for computing Energy)
			// TODO: CAUTION: The dirVec can change between timesteps (How do we deal with this)
			FKresults PassiveResults_OldTimestep = R; FKresults ActiveResults_OldTimestep = R;
			_PassiveMaterial->compute(PassiveResults_OldTimestep, _Fn, &_dirVec[0]);
			_ActiveMaterial->compute(ActiveResults_OldTimestep, Fen, &_dirVec[0]);
			
			// Compute \psi^{*} from the Force-velocity potential
			double psistar = _activation * _KineticPotential->Psi(deltaQ, _deltaT);
			// TODO: Add in deltaT * psistar to the Energy		// DONE
			
			R.W = _deltaT * phi + PassiveResults.W + ActiveResults.W - PassiveResults_OldTimestep.W - ActiveResults_OldTimestep.W + _deltaT * psistar;
		}

		if (R.request & FORCE)
		{
			// Compute Viscous Potential Stress
			Matrix3d dphidF = _ViscousPotential->dphidF(_Fn, Fnp1, _deltaT);

			Matrix3d ActiveTerm = Matrix3d::Zero();
			for (int d = 0; d < 3; d++)
				for (int e = 0; e < 3; e++)
					for (int b = 0; b < 3; b++)
						ActiveTerm(d,e) = ActiveTerm(d,e) + ActiveResults.P(d,b) * invFanp1(e,b);
			
			R.P = dphidF * _deltaT + PassiveResults.P + ActiveTerm;
		}
		
		if (R.request & STIFFNESS)
		{
			// Compute Viscous Potential Stiffness
			FourthOrderTensor d2phidF2 = _ViscousPotential->d2phidF2(_Fn, Fnp1, _deltaT);


			vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
			for (int p = 0; p < 3; p++)
				dFadQ[p] = A.exp() * _M[p] * _Fa;

			FourthOrderTensor dFedFanp1(3,3,3,3);
			for (int i = 0; i < 3; i++){
				for (int J = 0; J < 3; J++){
					for (int K = 0; K < 3; K++){
						for (int L = 0; L < 3; L++){
							for (int R = 0; R < 3; R++){
								dFedFanp1.set(i,J,K,L, dFedFanp1.get(i,J,K,L) + Fnp1(i,R) * invFanp1(R,K) * invFanp1(L,J));
							}
							dFedFanp1.set(i,J,K,L, -1. * dFedFanp1.get(i,J,K,L));
						}
					}
				}
			}

			FourthOrderTensor Term1(3,3,3,3);

			for (int i = 0; i < 3; i++) {
				for (int J = 0; J < 3; J++) {
					for (int k = 0; k < 3; k++) {
						for (int L = 0; L < 3; L++) {
							for (int N = 0; N < 3; N++) {
								for (int T = 0; T < 3; T++) {
									Term1.set(i,J,k,L, Term1.get(i,J,k,L) + ActiveResults.K.get(i,N,k,T) * invFanp1(J,N) * invFanp1(L,T));
								}
							}
							Term1.set(i,J,k,L, Term1.get(i,J,k,L) + PassiveResults.K.get(i,J,k,L));
						}
					}
				}
			}
			
			vector <Matrix3d> Term2(3, Matrix3d::Zero());
			vector <Matrix3d> Term2_2(3, Matrix3d::Zero());
			for (int i = 0; i < 3; i++) {
				for (int J = 0; J < 3; J++) {
					for (int a = 0; a < 3; a++) {
						for (int m = 0; m < 3; m++) {
							for (int N = 0; N < 3; N++) {
								for (int B = 0; B < 3; B++) {
									for (int s = 0; s < 3; s++) {
										for (int T = 0; T < 3; T++) {
											Term2[a](i,J) = Term2[a](i,J) + ActiveResults.K.get(i,B,s,T) * invFanp1(J,B) * dFedFanp1.get(s,T,m,N) * dFadQ[a](m,N);
										}
									}
									Term2_2[a](i,J) = Term2_2[a](i,J) + ActiveResults.P(i,m) * invFanp1(J,N) * invFanp1(B,m) * dFadQ[a](N,B);
								}
							}
						}
						// Term2_2[a](i,J) = -1 * Term2_2[a](i,J);
						Term2[a](i,J) = Term2[a](i,J) - Term2_2[a](i,J);
					}
				}
			}
			Matrix3d Term3 = computed2WdQ2().inverse();

			FourthOrderTensor Term234(3,3,3,3);
			for (int i = 0; i < 3; i++) {
				for (int J = 0; J < 3; J++) {
					for (int k = 0; k < 3; k++) {
						for (int L = 0; L < 3; L++) {
							for (int a = 0; a < 3; a++) {
								for (int b = 0; b < 3; b++) {
									Term234.set(i,J,k,L, Term234.get(i,J,k,L) + Term2[a](i,J) * Term3(a,b) * Term2[b](k,L));
								}
							}
							R.K.set(i,J,k,L, _deltaT * d2phidF2.get(i,J,k,L) + Term1.get(i,J,k,L) - Term234.get(i,J,k,L));
						}
					}
				}
			}
			// cout << _Fanp1 << endl;
			
		}
	}

	void PlasticMaterial::computeKinematicParameters()
	{
		for (int i = 0; i < _dirVec.size(); i++)
		{
			_M[i] = _dirVec[i] * _dirVec[i].transpose();
		}
	}
} // namespace voom
