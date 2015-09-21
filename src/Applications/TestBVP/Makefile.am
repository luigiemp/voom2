bin_PROGRAMS = TestBVP TestPressure TestMultiStepsBVP TestCube TestStrip StripLV
INCLUDES =		-I./ -I./../ -I./../../	-I./../../Mesh -I./../../Solver  	\
	 		-I./../../VoomMath/ -I./../../Shape -I./../../Quadrature	\
			-I./../../Model -I./../../Potential -I./../../Material 		\
			-I./../../Element 						\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3
AM_LDFLAGS = -L./ -L./../ -L./../../VoomMath/ -L./../../Model           	\
             -L./../../Element -L./../../Shape -L./../../Quadrature     	\
	     -L./../../VoomMath -L./../../Material -L./../../Mesh  -L./../../Solver  \
	     -L./../../Potential -L./../../ViscousPotential
LDADD = -lSolver -lModel -lMesh -lElement              \
	-lShape -lQuadrature -lVoomMath                \
	-lPotentials -lViscousPotentials -lMaterials		       \
	-lgfortran

TestBVP_SOURCES = TestBVP.cc
TestPressure_SOURCES = TestPressure.cc
TestMultiStepsBVP_SOURCES = TestMultiStepsBVP.cc
TestCube_SOURCES = TestCube.cc
TestStrip_SOURCES = TestStrip.cc
StripLV_SOURCES = StripLV.cc
