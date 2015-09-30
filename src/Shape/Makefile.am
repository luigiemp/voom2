## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=-I ./ 						\
	    -I ./../						\
	    -I ./../VoomMath -I./../Quadrature			\
	    -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3\
	    -I./../Geometry					\
	    -I./../HalfEdgeMesh

lib_LIBRARIES=libShape.a
libShape_a_SOURCES=Shape.cc HexShape.cc BarShape.cc MFShape.cc LMEShape.cc  \
		   MRKPMShape.cc LinTetShape.cc QuadTetShape.cc \
		   LinTriShape.cc QuadTriShape.cc LoopShellShape.cc
