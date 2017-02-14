bin_PROGRAMS = TestStrip EllipsoidLV CheckingConsistency MaterialPoint EllipsoidLVSpring StripSpringBC MaterialPoint_PurelyUniaxial RabbitLV Contraction CylindricalStrip

INCLUDES =		-I./							\
			-I./../							\
			-I./../../						\
			-I./../../Mesh						\
			-I./../../Solver  					\
	 		-I./../../VoomMath					\
			-I./../../Shape						\
			-I./../../Quadrature					\
			-I./../../Model						\
			-I./../../Potential					\
			-I./../../Material 					\
			-I./../../Material/MechanicsMaterial			\
			-I./../../Material/ViscousMaterial			\
			-I./../../Element 					\
			-I/u/local/apps/eigen/3.2.4/gcc-4.4.7/include/eigen3	\
		    	-I/u/local/apps/boost/1_59_0/gcc-4.4.7/include		\
			-I/u/local/apps/vtk/5.8.0/include/vtk-5.8

AM_LDFLAGS = 		-L./							\
			-L./../							\
			-L./../../VoomMath					\
			-L./../../Model           				\
             		-L./../../Element					\
			-L./../../Shape						\
			-L./../../Quadrature     				\
	     		-L./../../VoomMath					\
			-L./../../Material					\
			-L./../../Mesh						\
			-L./../../Solver  					\
	     		-L./../../Potential					\
	     		-L./../../Material/MechanicsMaterial			\
			-L./../../Material/ViscousMaterial			\
			-L/u/local/apps/vtk/5.8.0/lib/vtk-5.8

LDADD = 		-lSolver						\
			-lModel							\
			-lMesh							\
			-lElement 				             	\
			-lShape							\
			-lQuadrature						\
			-lVoomMath                				\
			-lPotentials						\
			-lMaterials						\
			-lMechanicsMaterial					\
			-lViscousMaterial					\
			-lgfortran						\
			-lvtkIO -lvtkGraphics -lvtkGenericFiltering -lvtkFiltering -lvtkCommon -lvtksys -ldl -lpthread -lvtkzlib -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio -lvtkNetCDF_cxx -lvtksqlite -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtkverdict

TestStrip_SOURCES = TestStrip.cc
EllipsoidLV_SOURCES = EllipsoidLV.cc
CheckingConsistency_SOURCES = CheckingConsistency.cc
MaterialPoint_SOURCES = MaterialPoint.cc
EllipsoidLVSpring_SOURCES = EllipsoidLVSpring.cc
StripSpringBC_SOURCES = StripSpringBC.cc
MaterialPoint_PurelyUniaxial_SOURCES = MaterialPoint_PurelyUniaxial.cc
RabbitLV_SOURCES = RabbitLV.cc
Contraction_SOURCES = Contraction.cc
CylindricalStrip_SOURCES = CylindricalStrip.cc
