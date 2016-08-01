## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     = 	-I./ 					\
		-I./../					\
	 	-I./../VoomMath/ 			\
		-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
# -I/u/home/l/luigiemp/project-cardio/Program/eigenInst/include/eigen3

lib_LIBRARIES=libPotentials.a libViscousPotentials.a
libPotentials_a_SOURCES = HillForceVelPotential.cc APForceVelPotential.cc APForceVelPotential_Test.cc
libViscousPotentials_a_SOURCES = NewtonianViscousPotential.cc
