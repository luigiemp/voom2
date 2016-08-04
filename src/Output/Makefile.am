## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 						\
		-I./../						\
	 	-I./../VoomMath	 					\
		-I./../Quadrature					\
		-I./../Shape						\
		-I./../Material						\
		-I./../Element						\
		-I./../Mesh 						\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
		-I./../Geometry -I./../HalfEdgeMesh			\
		-I./../Model						\
		-I/u/local/apps/vtk/5.8.0/include/vtk-5.8

LDADD	     =  -lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon -lvtksys -ldl -lpthread -lvtkzlib -lvtkexpat -lvtkCommon -lvtkHybrid -lvtkImaging -lvtkFiltering -lvtkRendering                      

lib_LIBRARIES = libOutput.a
libOutput_a_SOURCES = VTKOutput.cc 
