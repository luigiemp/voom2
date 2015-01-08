bin_PROGRAMS = TestModel
INCLUDES = -I./							\
	   -I./../ 						\
	   -I./../../  						\
	   -I./../Model 					\
	   -I./../../VoomMath					\
	   -I./../../Quadrature					\
	   -I./../../Shape 					\
	   -I./../../Material 					\
	   -I./../../Element					\
	   -I./../../Mesh					\
	   -I/u/local/apps/eigen/current/include/eigen3		\
	   -I/u/local/apps/trilinos/current/include
AM_LDFLAGS = 	-L./ 								\
	     	-L./../ 							\
	     	-L./../../Model					 		\
		-L./../../Mesh							\
		-L./../../Element 						\
		-L./../../Material 						\
		-L./../../Quadrature 						\
		-L./../../Shape 						\
		-L./../../VoomMath 						\
	     	-L/u/local/apps/trilinos/current/lib				\
	     	-L/u/local/compilers/intel/current/current/mkl/lib/em64t 	\
	     	-L/u/local/apps/mumps/4.10.0/intel-11.1/openmpi-1.4/lib  	\
	     	-L/u/local/compilers/intel/11.1/080/mkl/lib/em64t  	 	\
	     	-L/u/local/compilers/intel/11.1/080/mkl/lib/em64t
LDADD   = -lModel -lMesh -lElement -lQuadrature -lShape -lMaterials -lVoomMath	\
	  -lm -lepetra -lifpack -lepetraext -ltriutils -lpord			\
	  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lguide          	\
	  -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran      	\
	  -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack    	\
	  -lifcore -limf
TestModel_SOURCES = TestModel.cc
