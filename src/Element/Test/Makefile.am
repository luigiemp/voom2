bin_PROGRAMS = TestFEgeomElement TestLMEgeomElement

INCLUDES = -I./				\
	   -I./../			\
	   -I./../../			\
	   -I./../../VoomMath		\
	   -I./../../Shape		\
           -I./../../Quadrature		\
	   -I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS = -L./../			\
	     -L./../../VoomMath/	\
	     -L./../../Shape 		\
	     -L./../../Quadrature	\
	     -L./../Element		
LDADD   =  -lElement -lShape -lQuadrature -lVoomMath
TestFEgeomElement_SOURCES = TestFEgeomElement.cc
TestFEgeomElement1D_SOURCES = TestFEgeomElement1D.cc
#TestLMEgeomElement_SOURCES = TestLMEgeomElement.cc
