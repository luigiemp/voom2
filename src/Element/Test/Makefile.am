bin_PROGRAMS = TestFEgeomElement TestFEgeomElement1D TestLMEgeomElement

INCLUDES = -I./				\
	   -I./../			\
	   -I./../../			\
	   -I./../../VoomMath		\
	   -I./../../Shape		\
           -I./../../Quadrature		\
	   -I/u/local/apps/eigen/current/include/eigen3
AM_LDFLAGS = -L./../			\
	     -L./../../VoomMath/	\
	     -L./../../Shape 		\
	     -L./../../Quadrature	\
	     -L./../Element		
LDADD   =  -lElement -lShape -lQuadrature -lVoomMath
TestFEgeomElement_SOURCES = TestFEgeomElement.cc
TestFEgeomElement1D_SOURCES = TestFEgeomElement1D.cc
TestLMEgeomElement_SOURCES = TestLMEgeomElement.cc
