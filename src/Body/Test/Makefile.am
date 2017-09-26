bin_PROGRAMS = TestMechanicsBody TestPressureBody TestPotentialBody TestTorsionalBody TestPreferredNormalBody
INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../Body 						\
	   -I./../../VoomMath					\
	   -I./../../Solver/Results				\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Material/MechanicsMaterial			\
           -I./../../Material/PairMaterial			\
	   -I./../../Element					\
	   -I./../../Mesh					\
           -I./../../Geometry 					\
	   -I./../../HalfEdgeMesh				\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
	   -I/u/local/apps/vtk/5.8.0/include/vtk-5.8	

AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
	     	-L./../../Body					 		\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Material/MechanicsMaterial				\
		-L./../../Material/PairMaterial				\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../VoomMath 						\
		-L./../../Geometry						\
		-L./../../HalfEdgeMesh						\
		-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8
LDADD   = 	-lBody -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  	\
		-lGeometry -lHEMesh -lMechanicsMaterial -lPairMaterial					\
	 	-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon 			\
		-lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio 		\
		-lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict	
TestMechanicsBody_SOURCES = TestMechanicsBody.cc
TestPressureBody_SOURCES  = TestPressureBody.cc
TestPotentialBody_SOURCES  = TestPotentialBody.cc
TestTorsionalBody_SOURCES  = TestTorsionalBody.cc
TestPreferredNormalBody_SOURCES  = TestPreferredNormalBody.cc

