## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 						\
		-I./../						\
	 	-I./../VoomMath	 				\
		-I./../Quadrature				\
		-I./../Shape					\
		-I./../Material					\
		-I./../Element					\
		-I./../Mesh 					\
		-I/u/local/apps/trilinos/current/include	\
		-I/u/local/apps/eigen/current/include/eigen3      

lib_LIBRARIES = libModel.a
libModel_a_SOURCES = Model.cc EllipticModel.cc PoissonModel.cc MechanicsModel.cc
