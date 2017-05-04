bin_PROGRAMS = TestMechanicsBody
INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../Body 						\
	   -I./../../VoomMath					\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Material/MechanicsMaterial			\
	   -I./../../Element					\
	   -I./../../Mesh					\
	   -I./../../Solver/Results				\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
	   -I./../../Geometry -I./../../HalfEdgeMesh
AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
	     	-L./../../Body					 		\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Material/MechanicsMaterial				\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../VoomMath 						\
		-L./../../Geometry						\
		-L./../../HalfEdgeMesh
LDADD   = 	-lBody -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  \
		-lGeometry -lHEMesh -lMechanicsMaterial
TestMechanicsBody_SOURCES = TestMechanicsBody.cc
