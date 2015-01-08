bin_PROGRAMS = TestEpetraSparseSolver TestEpetraNRsolver

INCLUDES = -I./	-I./../ -I./../../ -I./../../Mesh -I./../Element	\
	   -I./../../VoomMath -I./../../Shape -I./../../Quadrature	\
	   -I./../../Model -I./../../Material -I./../../Element		\
	   -I/u/local/apps/eigen/current/include/eigen3			\
	   -I/u/local/apps/trilinos/current/include
AM_LDFLAGS = -L./ -L./../ -L./../../VoomMath/ -L./../../Model           	\
             -L./../../Element -L./../../Shape -L./../../Quadrature     	\
	     -L./../../VoomMath -L./../../Material -L./../../Mesh       	\
             -L/u/local/apps/trilinos/current/lib                       	\
             -L/u/local/apps/superlu_dist/3.1/intel-11.1/openmpi-1.4/lib	\
	     -L/u/local/compilers/intel/current/current/mkl/lib/em64t   	\
	     -L/u/local/apps/mumps/4.10.0/intel-11.1/openmpi-1.4/lib    	\
	     -L/u/local/compilers/intel/11.1/080/mkl/lib/em64t          	\
	     -L/u/local/compilers/intel/11.1/080/mkl/lib/em64t
LDADD = -lSolver -lModel -lMesh -lElement                       	\
	-lShape -lQuadrature  -lVoomMath                                \
	-lMaterials -lbelos -lifpack -lamesos  -laztecoo                	\
	-lteuchos  -lepetra  -lepetraext -ltriutils                     \
	-ldmumps -lmumps_common -lpord -lbelosepetra -lsuperlu_dist_3.1 \
	-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lguide            \
	-lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran        \
	-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack      \
	-lifcore -limf

TestEpetraSparseSolver_SOURCES = TestEpetraSparseSolver.cc
TestEpetraNRsolver_SOURCES = TestEpetraNRsolver.cc
