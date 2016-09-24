## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 						\
		-I./../						\
	 	-I./../VoomMath/ -I./../Shape			\
		-I./../Element/ -I./../Quadrature		\
		-I./../Mesh -I./../Material -I./../Model	\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
lib_LIBRARIES = libSolver.a
libSolver_a_SOURCES = EigenNRsolver.cc EigenNR_G_solver.cc LBFGSB.cc lbfgsb-routines.f # EpetraNRsolver.cc
