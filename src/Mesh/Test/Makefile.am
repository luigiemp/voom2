bin_PROGRAMS = TestMesh
INCLUDES = -I./../						\
	   -I./../../						\
	   -I./../../Solver/Results				\
	   -I./../../Quadrature					\
	   -I./../../VoomMath 					\
	   -I./../../Shape					\
	   -I./../../Element 					\
           -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3 \
	   -I./../../Geometry 					\
	   -I./../../HalfEdgeMesh				
AM_LDFLAGS = -L./../ 						\
	     -L./../../Quadrature 				\
	     -L./../../VoomMath 				\
	     -L./../../Shape					\
	     -L./../../Element					\
	     -L./../../Geometry					\
	     -L./../../HalfEdgeMesh 
LDADD      = -lMesh -lQuadrature -lVoomMath -lShape -lElement -lGeometry -lHEMesh 
TestMesh_SOURCES = TestMesh.cc

