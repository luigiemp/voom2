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
           -I./../../Geometry -I./../../HalfEdgeMesh		\
	   -I./../../Solver/Results				\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
	    -I/u/local/apps/vtk/5.8.0/include/vtk-5.8	

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
		-L./../../HalfEdgeMesh						\
		-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8
LDADD   = 	-lBody -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath  	\
		-lGeometry -lHEMesh -lMechanicsMaterial					\
	 	-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon 			\
		-lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio 		\
		-lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict	
TestMechanicsBody_SOURCES = TestMechanicsBody.cc
