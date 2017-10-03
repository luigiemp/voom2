bin_PROGRAMS 	= TestMaterial
INCLUDES =		-I./../					\
			-I./../../				\
			-I./../../../				\
			-I./../../ViscousMaterial		\	
	 		-I./../../../VoomMath/ 			\
			-I./../../../Potential/			\
			-I./../../../Geometry/			\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS 	= 	-L./../../                              \
			-L./../../../Potential/			\
	         	-L./../../../VoomMath/ 			\
			-L./../					\
			-L./../../ViscousMaterial		\
			-L./../../../Geometry
LDADD   = -lMaterials -lMechanicsMaterial -lViscousMaterial -lPotentials -lVoomMath -lGeometry  
TestMaterial_SOURCES = TestMaterial.cc
