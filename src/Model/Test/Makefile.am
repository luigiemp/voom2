bin_PROGRAMS = TestModel
INCLUDES = -I./							\
	   -I./../ 						\
           -I./../../						\
	   -I./../../VoomMath					\
	   -I./../../Solver/Results				\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Material/MechanicsMaterial			\
	   -I./../../Element					\
	   -I./../../Mesh					\
	   -I./../../Body					\
	   -I./../../Geometry					\
	   -I./../../HalfEdgeMesh				\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
	   -I/u/local/apps/vtk/5.8.0/include/vtk-5.8	
AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
		-L./../../VoomMath					 	\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../Material 						\
		-L./../../Material/MechanicsMaterial				\
		-L./../../Element 						\
		-L./../../Mesh							\
		-L./../../Body							\
		-L./../../Geometry						\
		-L./../../HalfEdgeMesh						\
		-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8
LDADD   = 	-lModel	-lBody -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  	\
		-lGeometry -lHEMesh -lMechanicsMaterial					\
	 	-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon 			\
		-lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio 		\
		-lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict
TestModel_SOURCES = TestModel.cc
