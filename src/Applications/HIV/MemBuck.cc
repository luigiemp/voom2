#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "GaugeLipid.h"
#include "GaugeModel.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv) {

    cout << " ------------------------------- " << endl;
    LoopShellMesh disc_mesh("circle_832_nodes.dat","circle_832_conn.dat");
    //LoopShellMesh disc_mesh("sphere_642_nodes.dat","sphere_642_conn.dat");
    
    uint NumMat = disc_mesh.getNumberOfElements();
    uint NodeDoF = 3;
    vector<GaugeLipid *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new GaugeLipid(1.0,0,0));

    GaugeModel model( &disc_mesh, materials, NodeDoF);

    // Set Requests
    uint PbDoF = (disc_mesh.getNumberOfNodes())*model.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
    
    //model.checkConsistency(myResults, 0, FORCE, 1.0e-6, 1.0e-5); 
    //return 0;

    // myResults.setRequest(ENERGY);
    // model.compute( myResults );
    // cout << "Energy : " << myResults.getEnergy()<<endl;
    // return 0;
    
    // myResults.setRequest(FORCE|ENERGY);
    // model.compute( myResults );
        
    // cout << "Energy : " << myResults.getEnergy()<<endl;
    // for (int i=0; i< PbDoF; i++){
    //   cout << (*(myResults._residual))(i) << endl;
    // 	//cout << "Detected at " << i << endl;
    // }
    //return 0;
    double Rbc = 1.35;
    vector<double> LowerBound, UpperBound;
    vector<int> BoundType;
    int numNodes = (model.getMesh())->getNumberOfNodes();

    int NumDOF = numNodes*3;
    LowerBound.resize(NumDOF);
    UpperBound.resize(NumDOF);
    BoundType.resize(NumDOF);
    for (int i = 0; i < numNodes; i++) { 
      const VectorXd X=(model.getMesh())->getX(i);
      BoundType[i*3] = 0; 
      BoundType[i*3+1] = 0;  
      BoundType[i*3+2] = 0;
      LowerBound[i*3] = 0;
      LowerBound[i*3+1] = 0;
      LowerBound[i*3+2] = 0;
      UpperBound[i*3] = 0;
      UpperBound[i*3+1] = 0;
      UpperBound[i*3+2] = 0;
      if ( (X(0)*X(0) + X(1)*X(1)) > Rbc*Rbc ){
	BoundType[i*3] = 2;  // Both lower and upper bound set
    	BoundType[i*3+1] = 2;  // Both lower and upper bound set
    	BoundType[i*3+2] = 2;  // Both lower and upper bound set
    	LowerBound[i*3] = X(0);
    	LowerBound[i*3+1] = X(1);
	LowerBound[i*3+2] = X(2);
    	UpperBound[i*3] = X(0);
    	UpperBound[i*3+1] = X(1);
	UpperBound[i*3+2] = X(2);
      }
      // if (X(0)*X(0)+X(1)*X(1)< 0.2*0.2 ){
      // 	BoundType[i*3+2] = 1; 
      // 	LowerBound[i*3+2] = 0.01;
      // }

    }

    LBFGSB mySolver( & model, & myResults, 5, 0, 1.0e-6, 98, 10000);
    mySolver.setBounds(BoundType, LowerBound, UpperBound);
    mySolver.solve();
    model.writeOutputVTK("ToDelete",0);
    cout << " ------------------------------ " << endl << endl;


  }
