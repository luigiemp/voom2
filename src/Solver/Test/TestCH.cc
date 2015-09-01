#include "EigenEllipticResult.h"

#include "LoopShellMesh.h"
#include "CahnHilliard.h"
#include "PhaseModel.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv) {

    cout << " ------------------------------- " << endl;
    //LoopShellMesh icosa_mesh("sphere_nodes_1SD.dat","sphere_conn_1SD.dat");
    //LoopShellMesh icosa_mesh("sphere_loop_nodes.dat","sphere_loop_conn.dat");
    LoopShellMesh icosa_mesh("T5sphere_nodes.dat","T5sphere_conn.dat");
    //LoopShellMesh icosa_mesh("nonicosa_sphere_nodes.dat","nonicosa_sphere_conn.dat");

    uint NumMat = icosa_mesh.getNumberOfElements();
    uint NodeDoF = 1;
    vector<CahnHilliard *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new CahnHilliard(1.0));

    PhaseModel model( &icosa_mesh, materials, 1);

    // Set Requests
    uint PbDoF = icosa_mesh.getNumberOfNodes();
    int TotNumMatProp = NumMat*1;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
    
    // ComputeRequest myRequest = ENERGY;
    // myResults.setRequest(myRequest);
    // model.compute( myResults );
    // cout << "Energy : " << myResults.getEnergy()<<endl;
    
    LBFGSB mySolver( & model, & myResults, 1, 1, 1.0e-6, 1, 50);
    mySolver.solve();
    cout << "Energy : " << myResults.getEnergy() << endl;
    model.printField();
    
    
    cout << " ------------------------------ " << endl << endl;


  }
