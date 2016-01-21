#include "EigenEllipticResult.h"

#include "LoopShellMesh.h"
#include "SCElastic.h"
#include "LoopShellModel.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv) {

    cout << " ------------------------------- " << endl;
    //LoopShellMesh icosa_mesh("sphere_nodes_1SD.dat","sphere_conn_1SD.dat");
    //LoopShellMesh icosa_mesh("sphere_loop_nodes.dat","sphere_loop_conn.dat");
    //LoopShellMesh icosa_mesh("T4sphere_nodes.dat","T4sphere_conn.dat");
    LoopShellMesh icosa_mesh("T3sphere_nodes.dat","T3sphere_conn.dat");
    //LoopShellMesh icosa_mesh("nonicosa_sphere_nodes.dat","nonicosa_sphere_conn.dat");

    uint NumMat = icosa_mesh.getNumberOfElements();
    uint NodeDoF = 3;
    vector<SCElastic *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new SCElastic(1,0,0));

    LoopShellModel model( &icosa_mesh, materials, NodeDoF);

    // Set Requests
    uint PbDoF = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
    //ComputeRequest myRequest = ENERGY;
    ///myResults.setRequest(myRequest);
    //model.compute( myResults );
    //cout << "Energy : " << myResults.getEnergy()<<endl;
    //model.checkConsistency(myResults, 0.0, FORCE, 1.0e-6, 1.0e-8); 
    //return 0;

    LBFGSB mySolver( & model, & myResults, 1, 1, 1.0e-5, 1, 100000);
    mySolver.solve();
    cout << " ------------------------------ " << endl << endl;


  }
