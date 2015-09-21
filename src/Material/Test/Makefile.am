bin_PROGRAMS 	= TestMaterial TestPlasticMaterial
INCLUDES =		-I./../					\
			-I./../../				\
	 		-I./../../VoomMath/ 			\
			-I./../../Potential/			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS 	= 	-L./../	                                \
			-L./../../Potential/			\
	         	-L./../../VoomMath/ 			\
			-L./../../ViscousPotential
LDADD   = -lMaterials -lPotentials -lVoomMath -lViscousPotentials   
TestMaterial_SOURCES = TestMaterial.cc
TestPlasticMaterial_SOURCES = TestPlasticMaterial.cc
