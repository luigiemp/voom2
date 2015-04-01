bin_PROGRAMS = TestEigen TestBVP TestLV TestPressure
INCLUDES =		-I./ -I./../ -I./../../	-I./../../Mesh -I./../Element		\
	 		-I./../../VoomMath/ -I./../../Shape -I./../../Quadrature	\
			-I./../../Model -I./../../Material -I./../../Element		\
			-I/u/home/l/luigiemp/project-cardio/Program/eigenInst/include/eigen3
AM_LDFLAGS = -L./ -L./../ -L./../../VoomMath/ -L./../../Model           	\
             -L./../../Element -L./../../Shape -L./../../Quadrature     	\
	     -L./../../VoomMath -L./../../Material -L./../../Mesh       
LDADD = -lSolver -lModel -lMesh -lElement              \
	-lShape -lQuadrature -lVoomMath                \
	-lMaterials 				       \
	-lgfortran

TestEigen_SOURCES = TestEigen.cc
TestBVP_SOURCES = TestBVP.cc
TestLV_SOURCES = TestLV.cc
TestPressure_SOURCES = TestPressure.cc