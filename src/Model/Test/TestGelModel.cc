// Test for Model

#include "GelModel.h"
#include "GelMesh.h"
#include "EigenResult.h"
#include "Spring.h"

using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF GEL MODEL " << endl << endl;
    
    
    GelMesh TestGelMesh("../../Mesh/Test/Filament2.node", "../../Mesh/Test/Filament2.ele");
    
    // Initialize Model
    uint NodeDoF = 3;

    uint NumFil = TestGelMesh.getNumberOfFilaments();
    
    std::cout << "Number of elements is " << NumFil << std::endl;
    
    vector<FilamentMaterial * > materials;
    materials.reserve(NumFil);
    
    Real k = 1.0;
    Vector3d d0;
    d0 << 1.0, 0.0, 0.0;
    
    for (int nFil = 0; nFil < NumFil; nFil++) {        
      Spring* Mat = new Spring(nFil,k,d0);
      materials.push_back(Mat);
    }
    
    GelModel myModel(&TestGelMesh, materials, NodeDoF,0,1);
    
     
    myModel.initializeField(1.0);    
    const vector<GeomFilament* > filElement = TestGelMesh.getFilaments();    

    for (int nFil = 0 ; nFil < NumFil ; nFil++){

      const vector<int > & NodesID = filElement[nFil]->getNodesID();
      const uint nodeNum = NodesID.size();
      vector<Vector3d > dlist(nodeNum-1,Vector3d::Zero());
      
      myModel.computeDeformation(dlist,filElement[nFil]);
      
      FilamentMaterial::Filresults Rf;
      Rf.request = 7;
      for(int n = 0; n<nodeNum-1; n++)
	{
	  cout << dlist[n] << endl;
	  vector<int> bond;
	  bond.push_back(NodesID[n]);
	  bond.push_back(NodesID[n+1]);
	  Vector3d d0;
	  cout << "OK" << endl;
	  d0[0] = TestGelMesh.getX(bond[0],0)-TestGelMesh.getX(bond[1],0);
	  d0[1] = TestGelMesh.getX(bond[0],1)-TestGelMesh.getX(bond[1],1);
	  d0[2] = TestGelMesh.getX(bond[0],2)-TestGelMesh.getX(bond[1],2);
	  cout << "d0" << endl;
	  cout << d0 << endl;
	  cout << dlist[n] << endl;
	  
	  materials[nFil]->compute(Rf,d0);
	  cout << "Energy     = " << Rf.W << endl;
	  cout << "Force  = " << Rf.f << endl;   
	  cout << "Stiffness = " << Rf.k << endl;
	}
    }    
    
    // Run consistency test                                           
    uint PbDoF = (TestGelMesh.getNumberOfNodes())*myModel.getDoFperNode();
    
    cout << "pbdof : " << PbDoF << endl;

    /*
    MechanicsModel myModel(&myFEmesh, materials, NodeDoF, PressureFlag, Pressure, &surfMesh);
    
    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);

    Real perturbationFactor = 0.1;
    uint myRequest = 6; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  

    */

    cout << endl << " END OF TEST OF MECHANICS MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;
  }


}
