bin_PROGRAMS 	= TestMesh
INCLUDES = -I./../					\
	   -I./../../					\
	   -I./../../VoomMath/ 				\
	   -I./../../Element -I./../../Quadrature	\
	   -I./../../Shape				\
	   -I/u/local/apps/eigen/current/include/eigen3
AM_LDFLAGS = -L./../ -L./../../VoomMath -L./../../Shape	\
	     -L./../../Quadrature -L ./../../Element
LDADD      = -lMesh -lElement -lShape -lQuadrature -lVoomMath
TestMesh_SOURCES = TestMesh.cc
