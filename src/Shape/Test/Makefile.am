bin_PROGRAMS = TestShape

INCLUDES = -I./	-I./../	-I./../../			\
	   -I./../../VoomMath -I./../../Quadrature	\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
TestShape_LDFLAGS = -L./../ -L./../../VoomMath/ -L./../../Quadrature
TestShape_LDADD   = -lShape -lVoomMath -lQuadrature
TestShape_SOURCES = TestShape.cc
