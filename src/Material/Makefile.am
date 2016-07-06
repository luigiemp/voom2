## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 					\
		-I./../					\
	 	-I./../VoomMath/ 			\
		-I./../Potential/			\
		-I./../Geometry/			\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
# -I/u/home/l/luigiemp/project-cardio/Program/eigenInst/include/eigen3

lib_LIBRARIES=libMaterials.a
libMaterials_a_SOURCES = MechanicsMaterial.cc CompNeoHookean.cc PassMyoA.cc ShellMaterial.cc SCElastic.cc LandauBrazovskii.cc PlasticMaterial.cc Jacobian.cc Humphrey.cc \
LinYinActive.cc LinYinActive_Compressible.cc Humphrey_Compressible.cc
