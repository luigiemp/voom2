//-*-C++-*-
#ifndef __MechanicsModel_h__
#define __MechanicsModel_h__

#include "Model.h"
#include "MechanicsMaterial.h"
#include "EigenResult.h"

// Include files for Writing Output:
#include <boost/lexical_cast.hpp>
#include <vtkVersion.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTensor.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkDoubleArray.h>

namespace voom{

  enum BCPotentialType {
    QUADRATIC, QUARTIC
  };

  // Model Results
  class MechanicsModel: public Model {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs.
     */
    MechanicsModel(Mesh* aMesh, vector<MechanicsMaterial * > Materials,
		   const uint NodeDoF,
		   int PressureFlag = 0, Mesh* SurfaceMesh = NULL,
		   int NodalForcesFlag = 0,
		   int _resetFlag = 1,
		   int _springBCflag = 0);

		   // const vector<string > & ElMatType,
		   // const map<string, MechanicsMaterial* > & ElMaterials);

    //! Input-file-based Constructor
    // MechanicsModel(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Destructor
    ~MechanicsModel() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	  UNIQUEmaterials.insert(_materials[i]);

      for (set<MechanicsMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++)
		delete (*it);
    };



    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();

      for (uint i = 0; i < numNodes; i++)
	for (uint j = 0; j < dim; j++)
	  _field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      // const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      const int nLocalDof = _field.size(); // Doing it this way for flexible membrane. Bad! Go back to the other way.
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }

    void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    void getField(vector<Real > & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    void setPrevField(vector<Real> & prevField) {
      _prevField = prevField;
    }

    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < _nodeDoF; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    void writeField(string OutputFile, int step) {
      // Create outputFile name
      stringstream FileNameStream;
      FileNameStream << OutputFile << step << ".dat";
      ofstream out;
      out.open( (FileNameStream.str()).c_str() );

      out << _field.size() << endl;
      for (uint i = 0; i < _field.size(); i++) {
	out << setprecision(15) << _field[i] << endl;
      }
      out.close();
    }

    uint getNumMat() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);

      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);
      // Assume all materials are of the same type

      return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<MechanicsMaterial * > getMaterials() {
      return _materials;
    }

    void setResetFlag(int ResetFlag) {
      _resetFlag = ResetFlag;
    }

    void setPressureFlag(int PressureFlag) {
      _pressureFlag = PressureFlag;
    }

    void setNodalForcesFlag(int NodalForcesFlag) {
      _nodalForcesFlag = NodalForcesFlag;
    }

    //! Write output
    void writeOutputVTK(const string OutputFile, int step);

    //! GET THIS OUT OF HERE:
    VTKCellType determineVTKCellType(int dim, int NodePerEl);

    //! Write VTK output for normals of pressure
    void writePressurePolyData(string OutputFile, int step);

    //! Write VTK output for linear spring
    void writeLinearSpringPolyData(string OutputFile, int step);

    //! Write VTK output for torsional spring
    void writeTorsionalSpringPolyData(string OutputFile, int step);

    //! Write VTK output for LJ Boundary Condition
    void writeLJBCPolyData(string OutputFile, int step);
    void writeFlexibleMembraneBCPolyData(string OutputFile, int step);
 
    //! Solve the system
    void compute(Result * R);

    //! Finalize Compute (Optional method which includes computations done at the end of a solve step)
    void finalizeCompute();

    // Apply pressure
    void applyPressure(Result * R);

    // Update pressure
    void updatePressure(Real Pressure) {
      _pressure = Pressure;
    }

    // Update nodal forces
    void updateNodalForces(vector<int > * ForcesID, vector<Real > * Forces) {
      _forcesID = ForcesID; _forces = Forces;
    }

    // Check consistency of gradg and Hg
    void checkDmat(EigenResult * R, Real perturbationFactor, Real h, Real tol);

    // Functions for applying spring BC
    // Initialize _springNodes (nodes at which spring BC are applied) and _springElements (elements connected to spring nodes)
    void initSpringBC(const string SpNodes, Mesh* SpMesh, Real SpringK);
    void computeSpringNormals();
    vector<Triplet<Real > > applySpringBC(Result & R);
    void computeAnchorPoints();		// This is for \bar{x} formulation

    // Functions for applying Torsional Spring BC
    void initTorsionalSpringBC(const string torsionalSpringNodes, Real torsionalSpringK);
    //! Computes centroid (x,y) of all the nodes but right now assumes long axis is z.
    void computeCentroid();
    void computeTangents();
    vector<Triplet<Real> > applyTorsionalSpringBC(Result & R);

    Real computeRefVolume();
    Real computeCurrentVolume();

    // Functions for a Lennard-Jones type boundary condition
    void initializeLennardJonesBC(BCPotentialType BCType, const string bodyPotentialBoundaryNodesFile, Mesh* rigidPotentialBoundaryMesh, Mesh* bodyPotentialBoundaryMesh, Real searchRadius, Real depthPotentialWell, Real minDistance);
    void initializeMembraneStiffness(double membraneStiffness);
    vector<Triplet<Real> > imposeLennardJones(Result& R);
    vector<Triplet<Real> > imposeFlexibleMembrane(Result & R);
    void computeLJNormals();
    void findNearestRigidNeighbors();
    void setLJStiffness(double stiffness){_depthPotentialWell = stiffness;};
    void recomputeAverageMinDistance();
    void toggleConstantMinimumDistanceFlag(bool constminDistFlag) {_useConstantMinimumDistance = constminDistFlag;};
    void toggleNumberNeighborFlag(bool useNumberNeighborFlag, int numberNearestNeighborsToUse=1) {
      _useNumberNeighborFlag = useNumberNeighborFlag;
      _numberNearestNeighborsToUse = numberNearestNeighborsToUse;
      this->findNearestRigidNeighbors();
    };
    void initiatePreTensionInMembrane(double preTensionFactor) {_preTensionFactor = preTensionFactor;};

    //! Lagrange multiplier method
    void initializeLagrangeMultiplierMethod(Mesh* EndocardialSurfMesh, Mesh* surfaceCapMesh, vector<int> endoBaseRingNodeSet, double vd);
    void setLagrangeMultiplier(double lambda) {
      int dimMainMesh = _myMesh->getDimension();
      int lagrangeMultiplierIndex = _myMesh->getNumberOfNodes() * dimMainMesh;
      if (_makeFlexible) lagrangeMultiplierIndex += _rigidPotentialBoundaryMesh->getNumberOfNodes() * _rigidPotentialBoundaryMesh->getDimension();
      _field[lagrangeMultiplierIndex] = lambda;
    };
    double getLagrangeMultiplier() {
      int dimMainMesh = _myMesh->getDimension();
      int lagrangeMultiplierIndex = _myMesh->getNumberOfNodes() * dimMainMesh;
      if (_makeFlexible) lagrangeMultiplierIndex += _rigidPotentialBoundaryMesh->getNumberOfNodes() * _rigidPotentialBoundaryMesh->getDimension();
      return _field[lagrangeMultiplierIndex];
    };
    vector<Triplet<Real> > imposeLagrangeMultiplier(Result & R);
    double computeCavityVolume();
    void setTargetVolume(double vd){_vd = vd;};

  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl);

    //! Compute Green Lagrangian Strain Tensor
    void computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl);

    //! List of Material data at each QP in the model
    vector<MechanicsMaterial * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    vector<Real > _field;

    // It should not be done here - maye we should have bodies and forms models from bodies
    int _pressureFlag;
    Real _pressure;
    Mesh* _surfaceMesh;

    int _nodalForcesFlag;
    vector<int > * _forcesID;
    vector<Real > * _forces;

    vector<Real > _prevField;
    int _resetFlag;

    // Spring BC
    int _springBCflag;
    vector<int > _spNodes;
    Mesh* _spMesh;
    Real _springK;
    vector<vector<int > > _spNodesToEle;
    vector<Vector3d > _spNormals;

    // Torsional Spring BC
    int _torsionalSpringBCflag;
    vector<int> _torsionalSpringNodes;
    Vector3d _centroidLocation;
    Real _torsionalSpringK;
    vector<Vector3d> _spTangents;

    // Members for imposing Lennard-Jones type boundary condition
    int _lennardJonesBCFlag;
    BCPotentialType _potentialType;			// Quadratic, Quartic, LJ
    Mesh* _rigidPotentialBoundaryMesh;
    vector<int> _bodyPotentialBoundaryNodes;
    Mesh* _bodyPotentialBoundaryMesh;
    Real _searchRadius;
    Real _depthPotentialWell;
    Real _minDistance;
    vector<vector<int> > _LJNodesToEle; 	// Stores the elements which are connected to a specific node
    vector<Vector3d> _rigidSurfaceNormals;
    vector<vector<int> > _rigidNeighbors;	// Stores each node and its rigid neighbors within a search radius
    bool _useNumberNeighborFlag;		// Flag if true, uses the nearest single neighbor rather than a set of neighbors
    int _numberNearestNeighborsToUse;		// Number of Neighbors to use
    vector<double> _minDistancePerRigidNode;	// For a new formulation where the minimum distance is different for each rigid node
    bool _useConstantMinimumDistance;		// Flag if true, uses the same minimum distance for every node
    bool _makeFlexible;				// Makes the boundary condition flexible
    double _membraneStiffnessCoefficient;		// Used for surface in plane membrane stiffness
    vector<Real> _membraneField;		// Field vector of just the membrane nodes
    double _preTensionFactor;			// preTensionFactor

    // Members for the Lagrange multiplier method and the Windkessel
    // NOT AT ALL THE BEST WAY TO DO THIS:
    bool _lagrangeMultiplierFlag;
    Mesh* _EndocardialSurfMesh;
    Mesh* _surfaceCapMesh;
    vector<int> _endoBaseRingNodeSet;
    vector<pair<int,int> > _ringAndMidsideNodePairs;	// This maps a ring node (first) to it's corresponding midside node (second) in the surfaceCapMesh
    
    double _vd;	// End-diastolic volume
  };

} // namespace voom

#endif
