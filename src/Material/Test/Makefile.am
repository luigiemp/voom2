bin_PROGRAMS 	= TestMaterial
INCLUDES =		-I./../					\
			-I./../../				\
	 		-I./../../VoomMath/ 			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3\
			-I./../../Geometry
TestMaterial_LDFLAGS = 	-L./../	                                \
	         	-L./../../VoomMath/ 			\
			-L./../../Geometry
TestMaterial_LDADD   = -lMaterials -lVoomMath  -lGeometry
TestMaterial_SOURCES = TestMaterial.cc
