// Test for Mechanics Body

#include "MechanicsBody.h"
#include "FEMesh.h"
#include "EigenResult.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"

using namespace voom;

int main(int argc, char** argv) {

  //Test Mechanics model
  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST OF MECHANICS BODY " << endl << endl;
  
    cout << " TEST OF 1st MECHANICS BODY " << endl << endl;

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh myFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);   
    uint NumQP = 1;

    

    int NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    for (int k = 0; k < NumMat; k++) {
      vector<Vector3d > Fibers;
      Vector3d N; N << 1.0, 0.0, 0.0;
      Fibers.push_back(N);
      // PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, Fibers);

      CompNeoHookean* Mat = new CompNeoHookean(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX);
      materials.push_back(Mat);
    }
    
    MechanicsBody myBody(&myFEmesh, &myState, materials);
    int PbDoF = myState.getDOFcount();
    cout << "PbDoF = " << PbDoF << endl;
    // myState.printPhi();
    myBody.writeOutputVTK("myBodyTestOutput", 1);



    EigenResult myResult(PbDoF, NumMat*(materials[0]->getMaterialParameters()).size() );
    cout << endl << "Material parameters per element = " << (materials[0]->getMaterialParameters()).size() << endl;
    myResult.initializeResults(PbDoF*4);



    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    myBody.compute(&myResult);
    cout << endl << "Body energy is  = " << myResult.getEnergy() << endl;
    // Change field and recompute energy
    for (int i=0; i<PbDoF; i++) {
      myState.linearizedUpdate(i, double(0.1*rand())/RAND_MAX);
    }
    myResult.resetResults(ENERGY);
    myBody.compute(&myResult);
    cout << endl << "Body energy is  = " << myResult.getEnergy() << endl;

    myBody.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);
    myBody.checkDmat(&myResult, perturbationFactor, myH, myTol);
    

    
    cout << endl << " END OF TEST OF 1st MECHANICS BODY " << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }






//   ////////////////////////////////////////////////







//   //Test Mechanics model
//   {
//     cout << " ------------------------------- " << endl;
//     cout << " TEST OF 2nd MECHANICS MODEL " << endl << endl;

//     FEMesh myFEmesh("/u/home/l/luigiemp/project-cardio/CardiacModels/Contraction/SmallHeart/Small_A_Cavity.node", 
//     		    "/u/home/l/luigiemp/project-cardio/CardiacModels/Contraction/SmallHeart/Small_A_Cavity.ele");
    
//     // Initialize Model
//     uint NodeDoF = 3;
//     uint NumQP = 4;

//     uint NumMat = myFEmesh.getNumberOfElements()*NumQP;
//     vector<MechanicsMaterial * > materials;
//     materials.reserve(NumMat);
    
//     for (int k = 0; k < NumMat; k++) {
//       Jacobian* Mat = new Jacobian(k);
//       materials.push_back(Mat);
//     }
   
//     MechanicsModel myModel(&myFEmesh, materials, NodeDoF);

 
//     // Run consistency test
//     uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
//     int TotNumMatProp = NumMat;
//     EigenEllipticResult myResults(PbDoF, TotNumMatProp);
 
//     Real perturbationFactor = 0.1;
//     uint myRequest = 2; // Check Forces
//     Real myH = 1e-6;
//     Real myTol = 1e-7;

//     myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

//     // Check model volume
//     cout << endl << "Model reference volume is       = " << myModel.computeRefVolume() << endl;
//     cout << endl << "Model current volume is         = " << myModel.computeCurrentVolume() << endl;
//     myResults.setRequest(1);
//     myModel.compute(myResults);
//     cout << endl << "Model energy = model volume is  = " << myResults.getEnergy() << endl;

//     // Change field and recompute volume
//     vector<Real > x(PbDoF, 0.0);
//     myModel.getField(x);
//     for (int i=0; i<PbDoF; i++) {
//       myModel.setField(i, x[i]*0.5);
//     }
//     cout << endl << "Model current volume is   = " << myModel.computeCurrentVolume() << endl;
//     myModel.compute(myResults);
//     cout << endl << "Model energy = model volume is  = " << myResults.getEnergy() << endl;
    
    
// <<<<<<< HEAD
//     cout << " ------------------------------- " << endl;
//     cout << " Testing Loop Shell Model " << endl << endl;
//     //LoopShellMesh icosa_mesh("sphere_nodes_1SD.dat","sphere_conn_1SD.dat");
//     //LoopShellMesh icosa_mesh("sphere_loop_nodes.dat","sphere_loop_conn.dat");
//     LoopShellMesh icosa_mesh("T5sphere_nodes.dat","T5sphere_conn.dat");
//     //LoopShellMesh icosa_mesh("nonicosa_sphere_nodes.dat","nonicosa_sphere_conn.dat");
//     uint NumMat = icosa_mesh.getNumberOfElements();
//     uint NodeDoF = 3;
//     vector<SCElastic *> materials;
//     materials.reserve(NumMat);
//     for(int k = 0; k < NumMat; k++)
//       materials.push_back(new SCElastic(1,0,0));

//     LoopShellModel model( &icosa_mesh, materials, NodeDoF);

//     // Run consistency test
//     uint PbDoF = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode();
//     int TotNumMatProp = NumMat*2;
//     EigenResult myResults(PbDoF, TotNumMatProp);
//     ComputeRequest myRequest = FORCE;
//     myResults.setRequest(myRequest);
//     //model.printField();
//     model.compute( myResults);
//     //cout << "Energy : " << myResults.getEnergy()<<endl;
//     //model.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
//     cout << std::scientific;
//     //cout << (*(myResults._residual)).cwiseAbs().colwise().maxCoeff() << endl;
//     cout << (*(myResults._residual)) << endl;
//     int counter = 0;
//     for (int i=0; i< PbDoF; i++ ){
//       if (abs((*(myResults._residual))(i)) > .0001){
// 	counter++;
// 	//cout << "Detected at " << i << endl;
//       }
//     }
//     //cout << "counter = " << counter << endl;
//     cout << endl << " END OF TEST OF MECHANICS MODEL " << endl;
// =======

//     cout << endl << " END OF TEST OF 2nd MECHANICS MODEL " << endl;
// >>>>>>> 29309bbfa22d9ab99c53776e10f3040558987a65
//     cout << " ------------------------------ " << endl << endl;


//   }

}
