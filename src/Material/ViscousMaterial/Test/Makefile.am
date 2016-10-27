bin_PROGRAMS 	= TestViscousMaterial
INCLUDES =		-I./../					\
			-I./../../				\
			-I./../../../				\
	 		-I./../../../VoomMath/ 			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS 	= 	-L./../	                                \
			-L./../../				\
	         	-L./../../../VoomMath/ 			\
			-L./../ViscousMaterial
LDADD   = -lVoomMath -lViscousMaterial  
TestViscousMaterial_SOURCES = TestViscousMaterial.cc
