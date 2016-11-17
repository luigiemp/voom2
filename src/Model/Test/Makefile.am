bin_PROGRAMS = TestFoundationModel # TestSpringModel # TestModel
INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../Model 					\
	   -I./../../VoomMath					\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Material/MechanicsMaterial			\
	   -I./../../Element					\
	   -I./../../Mesh					\
	   -I./../../Solver					\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
	   -I./../../Geometry -I./../../HalfEdgeMesh
AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
	     	-L./../../Model					 		\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Material/MechanicsMaterial				\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../Solver						\
		-L./../../VoomMath 						\
		-L./../../Geometry						\
		-L./../../HalfEdgeMesh
LDADD   = 	-lModel -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  \
		-lGeometry -lHEMesh -lSolver -lMechanicsMaterial
TestFoundationModel_SOURCES = TestFoundationModel.cc
# TestModel_SOURCES = TestModel.cc
# TestSpringModel_SOURCES = TestSpringModel.cc
