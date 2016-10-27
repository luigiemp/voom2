#include "PlasticMaterial.h"

namespace voom
{
	void PlasticMaterial::optimizeInternalVariables()
	{
		_Qnp1 = _Q;
		
		int iter = 0;

		// cout << "Qn+1: " << _Qnp1(0) << " " << _Qnp1(1) << " " << _Qnp1(2) << endl;
		// cout << "| \t Iter \t | \t Norm(dW/dQ) \t | \t dQ \t |" << endl;
		for (iter = 0; iter < _maxIter; iter++)
		{
			Vector3d dWdQ = computedWdQ();
			Matrix3d d2WdQ2 = computed2WdQ2();

			Vector3d dQ = -1. * d2WdQ2.inverse() * dWdQ;
			_Qnp1 = _Qnp1 + dQ;
			
			// cout << "\t" << iter << "\t\t" << dWdQ.norm() << "\t\t" << dQ.norm() << endl;

			double error = sqrt(square(dQ(0)) + square(dQ(1)) + square(dQ(2)));

			if(error < _hardOptTOL)
			{
			        // cout << "Qn+1: " << _Qnp1(0) << " " << _Qnp1(1) << " " << _Qnp1(2) << endl;
				break;
			}
		}
		// cout << "Internal Variable Optimization Iterations: \t" << iter << endl;
		if (iter == _maxIter)
		{
			cout << "****** Maximum Newton iterations reached for hardening variable optimization. ******" << endl;
			exit(1);
		}
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
		        A += exp(deltaQ(i)) * _M[i];
		}

		Fanp1 = (A) * _Fa;

		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = Fanp1.inverse();
		Matrix3d Fenp1 = _Fnp1 * invFanp1;

		// Compute Active Stress
		FKresults ActiveResults;
		ActiveResults.request = 2;
		_ActiveMaterial->compute(ActiveResults, Fenp1);
		
		// Method 1: BEGIN
		vector<Matrix3d> Gammap(3, Matrix3d::Zero());
		vector<Matrix3d> Lambdap(3, Matrix3d::Zero());
		for (int p = 0; p < 3; p++) {
		  Gammap[p] = exp(deltaQ[p]) * _M[p] * _Fa;
		  Lambdap[p] = Fenp1 * Gammap[p] * invFanp1;
		}

		Vector3d dDsdQnp1(0.0,0.0,0.0);
		for (int p = 0; p < 3; p++)
		{
		  // Matrix3d PLambda = ActiveResults.P.transpose() * Lambdap[p];
		  dDsdQnp1(p) = -1.0 * (ActiveResults.P.transpose() * Lambdap[p]).trace();
		}
		// Method 1: END

		/* // Method 2: BEGIN
		vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
		for (int p = 0; p < 3; p++)
		{
		    // dFadQ[p] = A.exp() * _M[p] * _Fa;
		    // Modified: 2.23.2016
		    dFadQ[p] = exp(deltaQ[p]) * _M[p] * _Fa;
		}
		
		Vector3d dDsdQnp1(0,0,0);
		for (int a = 0; a < 3; a++)
			for (int k = 0; k < 3; k++)
				for (int L = 0; L < 3; L++)
					for (int i = 0; i < 3; i++)
						for (int J = 0; J < 3; J++)
							for (int p = 0; p < 3; p++)
								dDsdQnp1(a) = dDsdQnp1(a) + ActiveResults.P(i,J) * _Fnp1(i,p) * invFanp1(p,k) * invFanp1(L,J) * dFadQ[a](k,L);
		dDsdQnp1 = dDsdQnp1 * -1;
		// Method 2: END */

		// Compute \dpsi^{*}dQnp1 from the Force-velocity potential
		Vector3d dpsidQnp1 = _activation * _KineticPotential->DPsiDQ(deltaQ, _deltaT);
		// TODO: Add in deltaT * dpsidQnp1 to dWdQ	// DONE

		dWdQ = dDsdQnp1 + dpsidQnp1 * _deltaT;

	        // cout << dDsdQnp1 << endl << endl << dpsidQnp1 << endl; exit(1);
		// cout << ActiveResults.P << endl; exit(1);
	        // cout << Fenp1 << endl; exit(1);
		// cout << dWdQ << endl; exit(1);

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
			A += exp(deltaQ(i)) * _M[i];
		}

		Fanp1 = A * _Fa;

		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = Fanp1.inverse();

		Matrix3d Fenp1 = _Fnp1 * invFanp1;

		FKresults ActiveResults;
		ActiveResults.request = 6;

		_ActiveMaterial->compute(ActiveResults, Fenp1);


		// Method 1: BEGIN
		vector<Matrix3d> Gammap(3, Matrix3d::Zero());
		vector<Matrix3d> Lambdap(3, Matrix3d::Zero());
		vector<vector<Matrix3d> > Thetapq(3, vector<Matrix3d>(3, Matrix3d::Zero()));
		vector<vector<Matrix3d> > Psipq(3, vector<Matrix3d>(3, Matrix3d::Zero()));

		for (int p = 0; p < 3; p++) {
		  Gammap[p] = exp(deltaQ[p]) * _M[p] * _Fa;
		  Lambdap[p] = Fenp1 * Gammap[p] * invFanp1;
		}

		for (int p = 0; p < 3; p++){
		  for (int q = 0; q < 3; q++){
		    Thetapq[p][q] = Lambdap[q] * Gammap[p] * invFanp1 + Lambdap[p] * Gammap[q]  * invFanp1;
		  }
		  Psipq[p][p] = Lambdap[p];
		}

		Matrix3d d2DdQ2np1 = Matrix3d::Zero();
		Matrix3d actstiffnessterm1 = Matrix3d::Zero();
		for (int p = 0; p < 3; p++) {
		  for (int q = 0; q < 3; q++) {
		    for (int e = 0; e < 3; e++) {
		      for (int b = 0; b < 3; b++) {
			for (int g = 0; g < 3; g++) {
			  for (int j = 0; j < 3; j++) {
			    actstiffnessterm1(p,q) += ActiveResults.K(e,b,g,j) * Lambdap[p](e,b) * Lambdap[q](g,j);
			  }
			}
		      }
		    }
		    d2DdQ2np1(p,q) = actstiffnessterm1(p,q) + (ActiveResults.P.transpose() * Thetapq[p][q]).trace() - (ActiveResults.P.transpose() * Psipq[p][q]).trace();
		  }
		}
		// Method 1: END

		/* // Method 2: BEGIN
		vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
		for (int p = 0; p < 3; p++)
		{
		        // dFadQ[p] = A * _M[p] * _Fa;
		        dFadQ[p] = exp(deltaQ[p]) * _M[p] * _Fa;
		}
		Matrix3d expA = A;

		// FourthOrderTensor d2FadQ2(3,3,3,3);
		vector < vector < Matrix3d> > d2FadQ2new(3, vector<Matrix3d>(3, Matrix3d::Zero()));
		for (int p = 0; p < 3; p++)
		  d2FadQ2new[p][p] = exp(deltaQ[p]) * _M[p] * _Fa;

		FourthOrderTensor dFedFanp1(3,3,3,3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
			        // Matrix3d temp = expA * _M[i] * _M[j] * _Fa;
				for (int p = 0; p < 3; p++) {
					for (int q = 0; q < 3; q++) {
					        // d2FadQ2(i,j,p,q) = temp(p,q);
						for (int r = 0; r < 3; r++) {
						  dFedFanp1(i,j,p,q) += _Fnp1(i,r) * invFanp1(r,p) * invFanp1(q,j);
						}
						dFedFanp1(i,j,p,q) = -1. * dFedFanp1(i,j,p,q);
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
		                                        					term1(p,q) = term1(p,q) + ActiveResults.K(e,f,h,i) * dFedFanp1(e,f,c,d) * dFedFanp1(h,i,l,m) * dFadQ[p](c,d) * dFadQ[q](l,m);
											}
										}
										term2(p,q) = term2(p,q) + ActiveResults.P(e,f) * -1 * (dFedFanp1(e,c,h,i) * invFanp1(d,f) + dFedFanp1(e,h,c,d) * invFanp1(i,f)) * dFadQ[p](c,d) * dFadQ[q](h,i);
									}
								}
								// term3(p,q) = term3(p,q) + ActiveResults.P(e,f) * dFedFanp1(e,f,c,d) * d2FadQ2(c,d,p,q);
								term3(p,q) = term3(p,q) + ActiveResults.P(e,f) * dFedFanp1(e,f,c,d) * d2FadQ2new[p][q](c,d);
							}
						}
					}
				}
			}
		}
		
		// DELETE:
		// cout << term1 << endl << term2 << endl << term3 << endl;

		Matrix3d d2DdQ2np1 = term1 + term2 + term3;

		// Method 2: End */

		// Compute \dpsi^{*}dQnp1 from the Force-velocity potential
		Matrix3d d2psidQnp12 = _activation * _KineticPotential->D2PsiDQDQ(deltaQ, _deltaT);
		// TODO: Add in deltaT * d2psidQnp12 to d2DdQ2np1	// DONE
		
		d2WdQ2 = d2DdQ2np1 + _deltaT * d2psidQnp12;

		return d2WdQ2;
	}
  
        // This function computes the Fanp1 and Qnp1 prior to doing the Compute
        void PlasticMaterial::preComputeHelper(const Matrix3d & Fnp1)
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
		
		// cout << "Qn+1: " << _Qnp1 << endl;

		Matrix3d A = Matrix3d::Zero();
		for (int i = 0; i < 3; i++){
			A += deltaQ(i) * _M[i];
		}

		Fanp1 = A.exp() * _Fa;
		_Fanp1 = Fanp1;
	}

	void PlasticMaterial::compute(FKresults & R, const Matrix3d & Fnp1)
	{		
	        preComputeHelper(Fnp1);
		// cout << "Internal Variables: " << _Qnp1[0] << ", " << _Qnp1[1] << ", " << _Qnp1[2] << endl;
		// cout << "Active Def Gradient: " << endl << _Fanp1 << endl;
		Vector3d deltaQ = _Qnp1 - _Q;

		Matrix3d A = Matrix3d::Zero();
		for (int i = 0; i < 3; i++){
			A += deltaQ(i) * _M[i];
		}


		// Compute Elastic Deformation Gradient
		Matrix3d invFa = _Fa.inverse();
		Matrix3d invFanp1 = _Fanp1.inverse();

		Matrix3d Fen = _Fn * invFa;
		Matrix3d Fenp1 = Fnp1 * invFanp1;

		FKresults PassiveResults; FKresults ActiveResults;
		PassiveResults.request = R.request; ActiveResults.request = R.request;
		if (R.request & STIFFNESS)
		{
			// STIFFNESS requires computation of the FORCE of Active and Passive materials
		        if (!(R.request & FORCE))
			  {
				PassiveResults.request = PassiveResults.request | FORCE;
				ActiveResults.request = ActiveResults.request | FORCE;
			  }
		}
		_PassiveMaterial->compute(PassiveResults, Fnp1);
		_ActiveMaterial->compute(ActiveResults, Fenp1);
		
		if (R.request & ENERGY)
		{
			// Compute Viscous Potential Stress
			double phi = _ViscousPotential->phi(_Fn, Fnp1, _deltaT);

			// At the Old Timestep (Only needed for computing Energy)
			// TODO: CAUTION: The dirVec can change between timesteps (How do we deal with this)
			FKresults PassiveResults_OldTimestep = R; FKresults ActiveResults_OldTimestep = R;
			_PassiveMaterial->compute(PassiveResults_OldTimestep, _Fn);
			_ActiveMaterial->compute(ActiveResults_OldTimestep, Fen);
			
			// Compute \psi^{*} from the Force-velocity potential
			double psistar = _activation * _KineticPotential->Psi(deltaQ, _deltaT);
			// TODO: Add in deltaT * psistar to the Energy		// DONE
			
			// cout << "psistar: " << psistar << endl;

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
			_elasticStress = PassiveResults.P + ActiveTerm;
		}
		
		if (R.request & STIFFNESS)
		{
			// Compute Viscous Potential Stiffness
			FourthOrderTensor d2phidF2 = _ViscousPotential->d2phidF2(_Fn, Fnp1, _deltaT);


			vector <Matrix3d> dFadQ(3, Matrix3d::Zero());
			for (int p = 0; p < 3; p++)
			{
				// dFadQ[p] = A.exp() * _M[p] * _Fa;
				// Added Next two lines 2.23.2016
				// Matrix3d P;
				// P << _dirVec[0](0), _dirVec[1](0), _dirVec[2](0),
				//      _dirVec[0](1), _dirVec[1](1), _dirVec[2](1),
			 	//      _dirVec[0](2), _dirVec[1](2), _dirVec[2](2);

				// Matrix3d Dtemp = Matrix3d::Zero();
				// Dtemp(p,p) = exp(deltaQ[p]);
				// dFadQ[p] = P * Dtemp * P.transpose() * _Fa;
			        dFadQ[p] = exp(deltaQ[p]) * _M[p] * _Fa;
			}

			FourthOrderTensor dFedFanp1(3,3,3,3);
			for (int i = 0; i < 3; i++){
				for (int J = 0; J < 3; J++){
					for (int K = 0; K < 3; K++){
						for (int L = 0; L < 3; L++){
							for (int R = 0; R < 3; R++){
							  dFedFanp1(i,J,K,L) += Fnp1(i,R) * invFanp1(R,K) * invFanp1(L,J);
							}
							dFedFanp1(i,J,K,L) = -1. * dFedFanp1(i,J,K,L);
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
								  Term1(i,J,k,L) += ActiveResults.K(i,N,k,T) * invFanp1(J,N) * invFanp1(L,T);
								}
							}
							Term1(i,J,k,L) += PassiveResults.K(i,J,k,L);
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
											Term2[a](i,J) = Term2[a](i,J) + ActiveResults.K(i,B,s,T) * invFanp1(J,B) * dFedFanp1(s,T,m,N) * dFadQ[a](m,N);
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
								  Term234(i,J,k,L) += Term2[a](i,J) * Term3(a,b) * Term2[b](k,L);
								}
							}
							R.K(i,J,k,L) =  _deltaT * d2phidF2(i,J,k,L) + Term1(i,J,k,L) - Term234(i,J,k,L);
						}
					}
				}
			}
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
