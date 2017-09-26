## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 							\
		-I./../							\
	 	-I./../VoomMath	 					\
		-I./../Quadrature					\
		-I./../Shape						\
		-I./../Material						\
		-I./../Material/MechanicsMaterial			\
		-I./../Material/PairMaterial				\
		-I./../Element						\
		-I./../Mesh 						\
		-I./../Solver/Results					\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
		-I./../Geometry 					\
		-I./../HalfEdgeMesh					\
		-I/u/local/apps/vtk/5.8.0/include/vtk-5.8

AM_LDFLAGS  =	-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8
LDADD 	    =	-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon -lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio -lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict

lib_LIBRARIES = libBody.a
libBody_a_SOURCES = Body.cc MechanicsBody.cc PressureBody.cc PotentialBody.cc TorsionalBody.cc PreferredNormalBody.cc
