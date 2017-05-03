#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "GaugeLipid.h"
#include "GaugeModel.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv) {

    cout << " ------------------------------- " << endl;
    //LoopShellMesh icosa_mesh("sphere_212_nodes.dat","sphere_212_conn.dat");
    LoopShellMesh icosa_mesh("sphere_642_nodes.dat","sphere_642_conn.dat");
    //LoopShellMesh icosa_mesh("sphere_2562_nodes.dat","sphere_2562_conn.dat");
    //LoopShellMesh icosa_mesh("sphere_10242_nodes.dat","sphere_10242_conn.dat");
    //LoopShellMesh icosa_mesh("paraview_sphere_nodes.dat","paraview_sphere_conn.dat");

    uint NumMat = icosa_mesh.getNumberOfElements();
    uint NodeDoF = 3;
    vector<GaugeLipid *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new GaugeLipid(1.0,0,0));

    GaugeModel model( &icosa_mesh, materials, NodeDoF);

    // Set Requests
    uint PbDoF = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
    
    //model.checkConsistency(myResults, 1e-1, FORCE, 1.0e-5, 1.0e-6); 
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
    LBFGSB mySolver( & model, & myResults, 5, 0, 1.0e-6, 98, 20000);
    mySolver.solve();
    model.writeOutputVTK("ToDelete",0);
    cout << " ------------------------------ " << endl << endl;


  }
