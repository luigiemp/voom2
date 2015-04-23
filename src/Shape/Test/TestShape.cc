#include "HexShape.h"
#include "LinTetShape.h"
#include "QuadTetShape.h"
#include "LinTriShape.h"
#include "QuadTriShape.h"
#include "BarShape.h"
#include "LMEShape.h"
#include "MRKPMShape.h"
#include "LoopShellShape.h"

using namespace voom;

int main()
{
  cout << endl << "Testing voom shape classes... " << endl;
  cout << ".............................." << endl << endl;

  // HexShape testing
  {
    cout << "Testing HexShape" << endl;
    Vector3d Point3D = Vector3d::Zero();
    srand( time(NULL) );
    Point3D(0) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    Point3D(1) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    Point3D(2) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    
    HexShape HexShp(Point3D);
    HexShp.checkConsistency(Point3D);
    HexShp.checkPartitionUnity(Point3D);
  }

  // TetShape testing
  {
    Vector3d Point3D = Vector3d::Zero();
    srand( time(NULL) );
    Point3D(0) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    Point3D(1) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    Point3D(2) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    
    cout << "Testing LinTetShape" << endl;
    LinTetShape LinTetShp(Point3D);
    LinTetShp.checkConsistency(Point3D);
    LinTetShp.checkPartitionUnity(Point3D);

    cout << "Testing QuadTetShape" << endl;
    QuadTetShape QuadTetShp(Point3D);
    QuadTetShp.checkConsistency(Point3D);
    QuadTetShp.checkPartitionUnity(Point3D);
  }

  // TriShape testing
  {
    Vector2d Point2D = Vector2d::Zero();
    srand( time(NULL) );
    Point2D(0) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    Point2D(1) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    
    cout << "Testing LinTriShape" << endl;
    LinTriShape LinTriShp(Point2D);
    LinTriShp.checkConsistency(Point2D);
    LinTriShp.checkPartitionUnity(Point2D);

    cout << "Testing QuadTriShape" << endl;
    QuadTriShape QuadTriShp(Point2D);
    QuadTriShp.checkConsistency(Point2D);
    QuadTriShp.checkPartitionUnity(Point2D);
  }

  // BarShape testing
  {
    cout << "Testing BarShape" << endl;
    VectorXd Point1D = VectorXd::Zero(1);
    srand( time(NULL) );
    Point1D(0) = 2.0*(Real(rand())/RAND_MAX - 0.5);
    
    BarShape BarShp(Point1D);
    BarShp.checkConsistency(Point1D);
    BarShp.checkPartitionUnity(Point1D);
  }


  // LMEshape testing
  {
    cout << "Testing LMEshape (cube domain)" << endl;
    vector<VectorXd > Nodes(8, Vector3d::Zero(3));
    Nodes[0] << 0.0, 0.0, 0.0;
    Nodes[1] << 1.0, 0.0, 0.0;
    Nodes[2] << 1.0, 1.0, 0.0;
    Nodes[3] << 0.0, 1.0, 0.0;
    Nodes[4] << 0.0, 0.0, 1.0;
    Nodes[5] << 1.0, 0.0, 1.0;
    Nodes[6] << 1.0, 1.0, 1.0;
    Nodes[7] << 0.0, 1.0, 1.0;
    VectorXd Point(3);
    Point << 0.4, 0.5, 0.6;
    Real beta = 0.8, tol = 1.0e-12;
    uint maxIter = 20;
    LMEShape LMEshp(Nodes, Point, beta, tol, maxIter);
    LMEshp.update(Point);
  
    cout << "Print LME at [0.4, 0.5, 0.6]" << endl;
    for (uint i=0; i<LMEshp.getShapeFunctionNum(); i++)
      {
	cout << "Node = " << i << " N = " << LMEshp.getN(i) << " DN = ";
	for (uint j=0; j<Point.size(); j++) {
	  cout << LMEshp.getDN(i,j) << " ";
	};
	cout << endl;
      };
    
    // Check consistency
    Point(0) = (Real(rand())/RAND_MAX );
    Point(1) = (Real(rand())/RAND_MAX );
    Point(2) = (Real(rand())/RAND_MAX );
    cout << "Check consistency at Point " << Point.transpose() << endl;
    LMEshp.checkConsistency(Point, 1.0e-8, 1.0e-6);
    LMEshp.checkPartitionUnity(Point);
    LMEshp.checkReproducingCondition(Point, 1.0e-6);
  }


  {
    cout << "Testing LMEshape (tet domain)" << endl;
    vector<VectorXd > Nodes(4, Vector3d::Zero(3));
    Nodes[0] << 0.0, 0.0, 0.0;
    Nodes[1] << 1.0, 0.0, 0.0;
    Nodes[2] << 0.0, 1.0, 0.0;
    Nodes[3] << 0.0, 0.0, 1.0;
    VectorXd Point(3);
    Point << 0.2, 0.1, 0.3;
    Real beta = 0.8, tol = 1.0e-12;
    uint maxIter = 20;
    LMEShape LMEshp(Nodes, Point, beta, tol, maxIter);
    LMEshp.update(Point);
    
    cout << "Print LME at [0.2, 0.1, 0.3]" << endl;
    for (uint i=0; i<LMEshp.getShapeFunctionNum(); i++)
      {
	cout << "Node = " << i << " N = " << LMEshp.getN(i) << " DN = ";
	for (uint j=0; j<Point.size(); j++) {
	  cout << LMEshp.getDN(i,j) << " ";
	};
	cout << endl;
      };

    // Check consistency
    Point(0) = (Real(rand())/RAND_MAX )*0.25;
    Point(1) = (Real(rand())/RAND_MAX )*0.25;
    Point(2) = (Real(rand())/RAND_MAX )*0.25;
    cout << "Check consistency at Point " << Point.transpose() << endl;
    LMEshp.checkConsistency(Point, 1.0e-8, 1.0e-6);
    LMEshp.checkPartitionUnity(Point);
    LMEshp.checkReproducingCondition(Point, 1.0e-6);
  }



  // MRKPMshape testing
  {
    // using same points as LME
    const Real radius = 0.001, support = 1.414, supportHat = 1.0;
    VectorXd Point2D(2);
    Point2D(0) = Point2D(1) = 0.5;
    
    vector<VectorXd > Nodes2D(4, Vector2d::Zero(2));
    Nodes2D[0] << 0.0, 0.0;
    Nodes2D[1] << 1.0, 0.0;
    Nodes2D[2] << 1.0, 1.0;
    Nodes2D[3] << 0.0, 1.0;
    MRKPMShape MRKPMshp(Nodes2D, Point2D, support, radius, supportHat);
    cout << "Print MRKPM at [0.5, 0.5]" << endl;
    for (uint i=0; i<MRKPMshp.getShapeFunctionNum(); i++)
      {
	cout << "Node = " << i << " N = " << MRKPMshp.getN(i) << " DN = ";
	for (uint j=0; j<Point2D.size(); j++) {
	  cout << MRKPMshp.getDN(i,j) << " ";
	};
	cout << endl;
      };

    Point2D(0) = (Real(rand())/RAND_MAX );
    Point2D(1) = (Real(rand())/RAND_MAX );
    cout << "Check consistency at Point " << Point2D.transpose() << endl;
    MRKPMshp.checkConsistency(Point2D, 1.0e-8, 1.0e-6);
    MRKPMshp.checkPartitionUnity(Point2D);
    MRKPMshp.checkReproducingCondition(Point2D, 1.0e-7);

  }

   // LoopShell testing
  {
    Vector2d Point2D = Vector2d::Zero();
    srand( time(NULL) );
    Point2D(0) = (Real(rand())/RAND_MAX );
    Point2D(1) = (Real(rand())/RAND_MAX );
    Vector3i Valences(6,6,6);
    
    cout << "Testing LoopShellShape" << endl;
    LoopShellShape LpShlShp(12, Valences, Point2D);
    cout << "Check consistency at Point " << Point2D.transpose() << endl;
    LpShlShp.checkConsistency(Point2D);
    LpShlShp.checkPartitionUnity(Point2D);

  }
  
}
