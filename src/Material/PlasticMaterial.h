//-*-C++-*-
#ifndef __Plastic_Material_h__
#define __Plastic_Material_h__

#include "VoomMath.h"
#include "MechanicsMaterial.h"
#include "Potential.h"
#include "ViscousPotential.h"

namespace voom
{
  class PlasticMaterial: public MechanicsMaterial
  {
  public:
    // PlasticMaterial(int MatID, MechanicsMaterial _ElasticMaterials, Potential _KineticPotential, Potenial _ViscousPotential)
    PlasticMaterial(int MatID, MechanicsMaterial* ActiveMaterial, MechanicsMaterial* PassiveMaterial, Potential* KineticPotential, ViscousPotential* ViscPotential): MechanicsMaterial(MatID), _ActiveMaterial(ActiveMaterial), _PassiveMaterial(PassiveMaterial), _KineticPotential(KineticPotential), _ViscousPotential(ViscPotential){
      for (int i = 0; i < 3; i++) _M.push_back(Matrix3d::Zero(3,3));
      _maxIter = 100;
      _hardOptTOL = 1.0E-10;
      _activation = 0;
      _deltaT = 0.05;

      _elasticStress = Matrix3d::Zero(3,3);
    }

    //! Destructor
    ~PlasticMaterial(){;}

    //! Compute function
    // TODO: Change Compute Function to take three vectors.
    void compute(FKresults & R, const Matrix3d & F, Vector3d * fiber = NULL);

    //! Helper function prior to doing the compute
    void preComputeHelper(const Matrix3d & F);

    //! Functions for Internal Variable Optimization
    Vector3d computedWdQ();
    Matrix3d computed2WdQ2();

    //! Optimize Internal Variables
    void optimizeInternalVariables();

    //! Update Variables from n+1->n

    //! Set Direction Vectors
    void setDirectionVectors(vector<Vector3d> dirVec) {_dirVec = dirVec;}

    //! Get active deformation gradient at previous timestep
    Matrix3d getActiveDeformationGradient() {return _Fa;}

    //! Get current n+1 active deformation gradient
    Matrix3d getCurrentActiveDeformationGradient() {return _Fanp1;}

    //! Set active deformation gradient at previous timestep
    void setActiveDeformationGradient(Matrix3d Fa) {_Fa = Fa;}

    //! Get total deformation gradient at previous timestep
    Matrix3d getTotalDeformationGradient(){return _Fn;}

    //! Set total deformation gradient at previous timestep
    void setTotalDeformationGradient(Matrix3d F){_Fn = F;}

    //! Get Hardening Parameters
    Vector3d getHardeningParameters() {return _Q;}

    //! Get current n+1 Hardening Parameters
    Vector3d getCurrentHardeningParameters() {return _Qnp1;}

    //! Set Hardening Parameters
    void setHardeningParameters(Vector3d Q){_Q = Q;}

    //! Compute Kinematic Parameters
    void computeKinematicParameters();

    //! Set Timestep
    void setTimestep(double deltaT){_deltaT = deltaT;}

    //! Set Activation Multiplier for Kinetic Potential
    void setActivationMultiplier(double activation){_activation = activation;}

    //! Update State from n+1 -> n
    void updateStateVariables()
    {
      setHardeningParameters(_Qnp1);
      setTotalDeformationGradient(_Fnp1);
      setActiveDeformationGradient(_Fanp1);
    }

    Matrix3d getElasticStressOnly() {return _elasticStress;}

    // FUNCTIONS THAT MUST BE OVERRIDDEN
    //! Clone
    virtual MechanicsMaterial* clone() const {return NULL;}

    //! SetMaterialParameters function
    virtual void setMaterialParameters(const vector<Real > &){;}
    virtual void setInternalParameters(const vector<Real > &){;};
    virtual void setRegularizationParameters(const vector<Real > &){;}

     //! GetMaterialParameters function
    virtual vector<Real > getMaterialParameters(){vector<double> A(2,0); return A;}
    virtual vector<Real > getInternalParameters(){vector <Real> intParam(1,_Q(0)); return intParam;}
    virtual vector<Real > getRegularizationParameters(){vector<double> A(2,0); return A;}

    virtual bool HasHistoryVariables(){return true;}


  protected:
    int _matID;
    
  private:
    //! Active Material
    MechanicsMaterial* _ActiveMaterial;

    //! Passive Material
    MechanicsMaterial* _PassiveMaterial;

    //! Kinetic Potential
    Potential* _KineticPotential;

    //! Viscous Potential
    ViscousPotential* _ViscousPotential;

    //! New Active Deformation Gradient
    Matrix3d _Fanp1;

    //! Active Deformation Gradient
    Matrix3d _Fa;

    //! Old Total Deformation Gradient
    Matrix3d _Fn;

    //! New Total Deformation Gradient (TODO: Find a better way to implement this. Needed for internal variable optimization)
    Matrix3d _Fnp1;

    //! Old Internal Hardening Parameters
    Vector3d _Q;

    //! Internal Hardening Parameters (TODO: Find a better way to implement this. Needed for internal variable optimization)
    Vector3d _Qnp1;

    //! Returns the Elastic Stress (Meaning it doesn't include Viscous stress)
    Matrix3d _elasticStress;
    
    //! Fiber directions
    vector <Vector3d> _dirVec;

    //! Kinematic Parameters
    vector <Matrix3d> _M;

    //! Newton-Raphson Parameters for Optimizing Hardening Variables
    int _maxIter;
    double _hardOptTOL;

    //! Set Timestep
    double _deltaT;

    //! Activation Multiplier
    double _activation;

  }; // class PlasticMaterial

} // namespace voom

#endif
