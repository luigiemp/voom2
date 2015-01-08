bin_PROGRAMS 	= TestMaterial
INCLUDES =		-I./../					\
			-I./../../				\
	 		-I./../../VoomMath/ 			\
			-I/u/local/apps/eigen/current/include/eigen3
TestMaterial_LDFLAGS = 	-L./../	                                \
	         	-L./../../VoomMath/ 
TestMaterial_LDADD   = -lMaterials -lVoomMath    
TestMaterial_SOURCES = TestMaterial.cc
