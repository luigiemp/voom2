bin_PROGRAMS 	= TestGelMesh
INCLUDES = -I./../					\
	   -I./../../					\
	   -I./../../VoomMath/ 				\
	   -I./../../Element -I./../../Quadrature	\
	   -I./../../Shape				\
           -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS = -L./../ -L./../../VoomMath -L./../../Shape	\
	     -L./../../Quadrature -L ./../../Element
LDADD      = -lMesh -lElement -lShape -lQuadrature -lVoomMath
TestGelMesh_SOURCES =  TestGelMesh.cc
