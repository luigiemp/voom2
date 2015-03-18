bin_PROGRAMS 	= TestQuadrature
INCLUDES =	 	-I./ 					\
			-I./../					\
			-I./../../				\
	 		-I./../../VoomMath/ 			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
TestQuadrature_LDFLAGS = -L./../	                        \
	                 -L./../../VoomMath/ 
TestQuadrature_LDADD   = -lQuadrature -lVoomMath    
TestQuadrature_SOURCES = TestQuadrature.cc
