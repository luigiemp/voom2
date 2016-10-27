## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 							\
		-I./../							\
	 	-I./../VoomMath	 					\
		-I./../Quadrature					\
		-I./../Shape						\
		-I./../Material						\
		-I./../Material/MechanicsMaterial			\
		-I./../Element						\
		-I./../Mesh 						\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
		-I./../Geometry 					\
		-I./../HalfEdgeMesh

lib_LIBRARIES = libModel.a
libModel_a_SOURCES = Model.cc MechanicsModel.cc LoopShellModel.cc LBModel.cc PoissonModel.cc
