bin_PROGRAMS =  TestLB
INCLUDES =		-I./ -I./../ -I./../../	-I./../../Mesh -I./../Element		\
	 		-I./../../VoomMath/ -I./../../Shape -I./../../Quadrature	\
			-I./../../Model -I./../../Material -I./../../Element		\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3\
			-I./../../Geometry -I./../../HalfEdgeMesh
AM_LDFLAGS = -L./ -L./../ -L./../../VoomMath/ -L./../../Model           	\
             -L./../../Element -L./../../Shape -L./../../Quadrature     	\
	     -L./../../VoomMath -L./../../Material -L./../../Mesh -L./../../Geometry\
             -L./../../HalfEdgeMesh  
LDADD = -lSolver -lModel -lMesh -lElement              \
	-lShape -lQuadrature -lVoomMath                \
	-lMaterials -lGeometry	-lHEMesh	       \
	-lgfortran

TestEigen_SOURCES = TestEigen.cc
TestBVP_SOURCES = TestBVP.cc
TestLV_SOURCES = TestLV.cc
TestPressure_SOURCES = TestPressure.cc
TestLB_SOURCES = TestLB.cc