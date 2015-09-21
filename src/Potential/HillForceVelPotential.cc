#include "HillForceVelPotential.h"

namespace voom
{
	double HillForceVelPotential::Psi(Vector3d deltaQ, double deltaT)
	{
		double psi = - _b * ( _a + _S0) * log(-_b + deltaQ(0)/deltaT) - _a * (deltaQ(0)/deltaT);
		return psi;
	}

	Vector3d HillForceVelPotential::DPsiDQ(Vector3d deltaQ, double deltaT)
	{
		Vector3d dpsidQ = Vector3d::Zero();
    
		dpsidQ(0) = - _b * ( _a + _S0) * (1/deltaT)/(-_b+deltaQ(0)/deltaT) - _a * (1/deltaT);
		return dpsidQ;
	}

	Vector3d HillForceVelPotential::DPsiDQdot(Vector3d deltaQ, double deltaT){return Vector3d::Zero();}

	Matrix3d HillForceVelPotential::D2PsiDQDQ(Vector3d deltaQ, double deltaT)
	{
		Matrix3d D2PsiDQDQ = Matrix3d::Zero();
    
		D2PsiDQDQ(0,0) = _b * (_a + _S0)/pow(deltaT,2) * 1/pow((-_b+deltaQ(0)/deltaT),2);
		return D2PsiDQDQ;
	}

	Matrix3d HillForceVelPotential::D2PsiDQDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}

	Matrix3d HillForceVelPotential::D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}
}
