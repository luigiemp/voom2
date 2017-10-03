## Makefile.am -- Process this file with automake to produce Makefile.in
bin_PROGRAMS 	= TestPairMaterial
INCLUDES =		-I./../					\
			-I./../../				\
			-I./../../../				\
	 		-I./../../../VoomMath/ 			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS 	= 	-L./../	                                \
	         	-L./../../../VoomMath/ 			
LDADD   = -lPairMaterial -lVoomMath
TestPairMaterial_SOURCES = TestPairMaterial.cc

