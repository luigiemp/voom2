#include "EigenEllipticResult.h"

#include "LoopShellMesh.h"
#include "LandauBrazovskii.h"
#include "LBModel.h"
#include "LBFGSB.h"
#include <time.h>

using namespace voom;

int main(int argc, char** argv) {

    time_t start, end;
    double dif;
    time(&start);
    const string IC = "ToDelete0.vtk";
    //const string IC = "S5_l24_c0.01_-3.6_1_2.vtk0.vtk";
    cout << " ------------------------------- " << endl;
    //LoopShellMesh icosa_mesh("sphere_nodes_1SD.dat","sphere_conn_1SD.dat");
    //LoopShellMesh icosa_mesh("sphere_loop_nodes.dat","sphere_loop_conn.dat");
    //LoopShellMesh icosa_mesh("T3sphere_nodes.dat","T3sphere_conn.dat");
    //LoopShellMesh icosa_mesh("T4sphere_nodes.dat","T4sphere_conn.dat");
    LoopShellMesh icosa_mesh("T5sphere_nodes.dat","T5sphere_conn.dat");
    //LoopShellMesh icosa_mesh("exact_sphere_nodes.dat","exact_sphere_conn.dat");
    //LoopShellMesh icosa_mesh("T2psphere_nodes.dat","T2psphere_conn.dat");
    //LoopShellMesh icosa_mesh("nonicosa_sphere_nodes.dat","nonicosa_sphere_conn.dat");

    //-------------- Smoothened mesh --------------------------------
    //LoopShellMesh icosa_mesh("Und400_1622_Sm_nodes.dat","Und400_1622_Sm_conn.dat");
    //LoopShellMesh icosa_mesh("Und400_6482_Sm_nodes.dat","Und400_6482_Sm_conn.dat");
    //LoopShellMesh icosa_mesh("Und400_25922_Sm_nodes.dat","Und400_25922_Sm_conn.dat");
     //-------------- Triangle with zig-zag connections ---------------------
    //LoopShellMesh icosa_mesh("Und400_1622_nodes.dat","Und400_1622_conn.dat");
    //LoopShellMesh icosa_mesh("Und400_6482_nodes.dat","Und400_6482_conn.dat");
    //--------------  Triangulations with unidirectional connection -----------
    //LoopShellMesh icosa_mesh("Und400_2722_nodes.dat","Und400_2722_conn.dat");
    //LoopShellMesh icosa_mesh("Und400_1154_nodes.dat","Und400_1154_conn.dat");

    
    uint NumMat = icosa_mesh.getNumberOfElements();
    vector<LandauBrazovskii *> materials;
    materials.reserve(NumMat);
    double l = 15.5;
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new LandauBrazovskii(.0005,sqrt(l*(l+1)),-1,1.5,1)); //for l32 c=0.001

   
    LBModel model( &icosa_mesh, materials, 1);
    // Set Requests
    uint PbDoF = icosa_mesh.getNumberOfNodes();
    int TotNumMatProp = NumMat*1;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);

    //Load Initial configuration
    if (false){
    ifstream icID;
    double phi[PbDoF];
    try{
      icID.open(IC.c_str(), ios::in); //open initial condition file
      string line;
      //Header 1
      for (int cnt =0; cnt <5;cnt++) getline(icID,line);
      while(getline(icID,line)){
	//this while looks removes the nodal information
	double x,y,z;
	istringstream iss(line);
	if(!(iss >> x >> y >> z)) {break;}
      }
       while(getline(icID,line)){
	 //this while looks removes the connectivity info
	 double n1,n2,n3, n4;
	 istringstream iss(line);
	 if(!(iss >> n1 >> n2 >> n3 >>n4)) {break;}
       }
       for (int cnt =0; cnt <3;cnt++) getline(icID,line);
       for (int cnt=0; cnt < PbDoF; cnt++){
	 getline(icID,line);
	 istringstream iss(line);
	 iss >> phi[cnt];
	 cout << phi[cnt] << endl;
       }
       model.initializeField(phi);
       icID.close();
    }
    catch (ifstream::failure e){
      cout<< "Unknown input file"<<endl;
      cout<< "Proceeding with random initial condition ... "<<endl;
    }
    }
    //model.checkConsistency(myResults, 0.0, FORCE, 1.0e-5, 1.0e-8); 
    //return 0;

    int choice = 1;
    if (choice == 1){
      LBFGSB mySolver( & model, & myResults, 5, 0, 1e-8, 98, 4000);
      //LBFGSB mySolver( & model, & myResults, 5, 0.5, 1e-8, 98, 10000);
      mySolver.solve();
      //model.writeOutputVTK("U25922_l24_c0.001_-0.1_2_2.vtk",0);
      //model.writeOutputVTK("S5_l24_c0.01_-3.5_1_2.vtk",0);
      model.writeOutputVTK("ToDelete",0);
      //model.writeOutputVTK("l_7",0);
    }
    if (choice==2){
      // Dynamically evolve using explicit time stepping
      Real dt = 0.00001;
      Real T = 10;
      int printCtr = 1;
      ComputeRequest myRequest = FORCE;
      myResults.setRequest(myRequest);
      int counter = 0;
      model.writeOutputVTK("Udynamics_",counter);
      for (Real t=0; t<=T; t=t+dt){
	model.compute( myResults );
	for (int i=0; i< PbDoF; i++){
	  model.linearizedUpdate(i, -dt*myResults.getResidual(i));
	  //cout << " update: " <<  myResults.getResidual(i) <<endl;
	}
	counter++;
	if (counter % printCtr == 0){
	  cout << "counter = "<< counter << endl;
	  model.printField();
	  model.writeOutputVTK("Udynamics_",counter);
	}
	//if (counter > 2) exit(0);
      }
    }
      if (choice==3){
      // Dynamically evolve using implicit time stepping
	model.setImplicitDynamicsFlag(true);
	Real dt = 0.3;
	Real T = 10;
	int printCtr = 1;
	ComputeRequest myRequest = FORCE;
	myResults.setRequest(myRequest);
	int counter = 0;
	model.writeOutputVTK("dynamics_",counter);
	LBFGSB mySolver( & model, & myResults, 5, 0, 1e-6, 1, 50000);
	model.setPrevField();
      	for (Real t=0; t<=T; t=t+dt){
	  mySolver.solve();
	  model.setPrevField();
	  counter++;
	  if (counter % printCtr == 0){
	    cout << "counter = "<< counter << endl;
	    model.printField();
	    model.writeOutputVTK("dynamics_",counter);
	  }
	  //if (counter > 2) exit(0);
	}
      }
    time (&end);
    dif = difftime (end,start);
    cout << endl << "All done :) in " << dif  << " s" << endl;
    
}
