## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 					\
		-I./../					\
	 	-I./../VoomMath/ -I./../Shape		\
		-I./../Element/ -I./../Quadrature	\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3\
		-I./../Geometry	\
		-I./../HalfEdgeMesh

lib_LIBRARIES=libMesh.a
libMesh_a_SOURCES = Mesh.cc FEMesh.cc LoopShellMesh.cc
