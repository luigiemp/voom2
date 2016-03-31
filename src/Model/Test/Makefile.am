bin_PROGRAMS = TestGelModel
INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../Model 					\
	   -I./../../VoomMath					\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Element					\
	   -I./../../Mesh					\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
	     	-L./../../Model					 		\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../VoomMath 						
LDADD   = -lModel -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath
TestGelModel_SOURCES = TestGelModel.cc
