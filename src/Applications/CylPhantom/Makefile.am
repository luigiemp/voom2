bin_PROGRAMS = CylPhantom

INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../../Body 					\
	   -I./../../VoomMath					\
	   -I./../../Solver/Results				\
	   -I./../../Solver 					\
	   -I./../../Model					\
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
		-L./../../Model							\
	     	-L./../../Body					 		\
                -L./../../VoomMath						\
	        -L./../../Solver/Results					\
	        -L./../../Solver 						\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Material/MechanicsMaterial				\
		-L./../../Material/PairMaterial					\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../VoomMath 						\
		-L./../../Geometry						\
		-L./../../HalfEdgeMesh						\
		-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8

LDADD   = 	-lSolver -lModel -lBody -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  	\
		-lGeometry -lHEMesh -lMechanicsMaterial	-lPairMaterial -lgfortran			\
	 	-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon 			\
		-lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio 		\
		-lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict

CylPhantom_SOURCES = CylPhantom.cc
