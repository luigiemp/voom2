#include "APForceVelPotential.h"

namespace voom
{
	double APForceVelPotential::Psi(Vector3d deltaQ, double deltaT)
	{
		double vel = (deltaQ(0)/deltaT);
		double psi = deltaT * _S0/_a * exp(_a * vel);
		return psi;
	}

	Vector3d APForceVelPotential::DPsiDQ(Vector3d deltaQ, double deltaT)
	{
		double vel = (deltaQ(0)/deltaT);
		Vector3d dpsidQ = Vector3d::Zero();
    
		dpsidQ(0) = _S0 * exp(_a * vel);
		return dpsidQ;
	}

	Vector3d APForceVelPotential::DPsiDQdot(Vector3d deltaQ, double deltaT){return Vector3d::Zero();}

	Matrix3d APForceVelPotential::D2PsiDQDQ(Vector3d deltaQ, double deltaT)
	{
		double vel = (deltaQ(0)/deltaT);
		Matrix3d D2PsiDQDQ = Matrix3d::Zero();
    
		D2PsiDQDQ(0,0) = _a * _S0/deltaT * exp(_a * vel);
		return D2PsiDQDQ;
	}

	Matrix3d APForceVelPotential::D2PsiDQDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}

	Matrix3d APForceVelPotential::D2PsiDQdotDQdot(Vector3d deltaQ, double deltaT){return Matrix3d::Zero();}
}
