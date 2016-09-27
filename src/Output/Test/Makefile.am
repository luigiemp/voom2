## Makefile.am -- Process this file with automake to produce Makefile.in

bin_PROGRAMS = TestVTKOutput

INCLUDES     = 	-I./ 							\
		-I./../							\
		-I./../../						\
	 	-I./../../VoomMath	 				\
		-I./../../Quadrature					\
		-I./../../Shape						\
		-I./../../Material					\
		-I./../../Element					\
		-I./../../Mesh 						\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
		-I./../../Geometry -I./../HalfEdgeMesh			\
		-I./../../Model						\
		-I/u/local/apps/vtk/5.8.0/include/vtk-5.8		\
		-I./../Output

AM_LDFLAGS   =  -L/u/local/apps/vtk/5.8.0/lib/vtk-5.8			\
		-L./                                                            \
                -L./../                                                         \
		-L./../../							\
                -L./../../Model                                                 \
                -L./../../Mesh                                                  \
                -L./../../Element                                               \
                -L./../../Material                                              \
                -L./../../Quadrature                                            \
                -L./../../Shape                                                 \
                -L./../../VoomMath                                              \
                -L./../../Geometry                                              \
                -L./../../HalfEdgeMesh						\
		-L./../Output

LDADD	     =  -lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon -lvtksys -ldl -lpthread -lvtkzlib -lvtkexpat -lvtkCommon -lvtkHybrid -lvtkImaging -lvtkFiltering -lvtkRendering -lOutput

TestVTKOutput_SOURCES = TestVTKOutput.cc
