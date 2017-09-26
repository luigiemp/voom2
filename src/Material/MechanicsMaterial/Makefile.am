## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 					\
		-I./../					\
		-I./../../				\
	 	-I./../../VoomMath/ 			\
		-I./../../Potential/			\
		-I./../ViscousMaterial			\
		-I./../../Geometry/			\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
# -I/u/home/l/luigiemp/project-cardio/Program/eigenInst/include/eigen3

lib_LIBRARIES=libMechanicsMaterial.a
libMechanicsMaterial_a_SOURCES = MechanicsMaterial.cc CompNeoHookean.cc Jacobian.cc PassMyoA.cc Holzapfel.cc Guccione.cc PlasticMaterial.cc Humphrey.cc \
LinYinActive.cc LinYinActive_Compressible.cc Humphrey_Compressible.cc ImposedKinematics.cc
