## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 						\
		-I./../						\
	 	-I./../VoomMath/ -I./../Shape			\
		-I./../Element/ -I./../Quadrature		\
		-I./../Mesh -I./../Material -I./../Model	\
		-I/u/local/apps/trilinos/current/include	\
		-I/u/local/apps/eigen/current/include/eigen3      

lib_LIBRARIES = libSolver.a
libSolver_a_SOURCES = EpetraSparseSolver.cc EpetraNRsolver.cc
