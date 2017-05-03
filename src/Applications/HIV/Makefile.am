bin_PROGRAMS = MemBuck
INCLUDES =		-I./ -I./../ -I./../../	-I./../../Mesh -I./../Element		\
	 		-I./../../VoomMath/ -I./../../Shape -I./../../Quadrature	\
			-I./../../Model -I./../../Material -I./../../Element		\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3\
			-I./../../Geometry -I./../../HalfEdgeMesh	\
			-I./../../Solver -I./continuation 
AM_LDFLAGS = -L./ -L./../ -L./../../VoomMath/ -L./../../Model           	\
             -L./../../Element -L./../../Shape -L./../../Quadrature     	\
	     -L./../../VoomMath -L./../../Material -L./../../Mesh -L./../../Geometry\
             -L./../../HalfEdgeMesh -L./../../Solver
LDADD = -lSolver -lModel -lMesh -lElement              \
	-lShape -lQuadrature -lVoomMath                \
	-lMaterials -lGeometry	-lHEMesh	       \
	-lgfortran  -lgsl -lgslcblas

MemBuck_SOURCES = MemBuck.cc

