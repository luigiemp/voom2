bin_PROGRAMS = test

INCLUDES = -I./	-I./../	-I./../../			\
	   -I./../../VoomMath -I./../../Quadrature	\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
test_LDFLAGS = -L./../ -L./../../VoomMath/ -L./../../Quadrature
test_LDADD   = -lHEMesh -lVoomMath -lQuadrature
test_SOURCES = test.cc
