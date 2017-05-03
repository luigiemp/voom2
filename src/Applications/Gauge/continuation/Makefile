ARPACKPP_DIR=/u/home/s/sanjaydm/arpack++
include LAPACKMake.inc
include $(ARPACKPP_DIR)/Makefile.inc

CBLAS= -lgsl -lgslcblas

DBG = #-g #Debugger Option
OPT = -O0 -pg -fPIC
# The libraries and include directories below are produced by pkg-config --cflags --libs plplotd-c++

DENSEMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym
MY_FLAGS = -g -Dlinux -I/u/home/s/sanjaydm/arpack++/include -I -Wall -ansi 

uvw: uvw_fem.cpp uvw_model2_FV.o uvw_model2_SV.o continuation.o kbhit.o myfunctions.o arnoldi_eig.o
	g++ -c $(MY_FLAGS) -I$(DENSEMAT_DIR) uvw_fem.cpp  arnoldi_eig.o
	g++ uvw_fem.o uvw_model2_FV.o uvw_model2_SV.o continuation.o  arnoldi_eig.o kbhit.o myfunctions.o\
         $(CBLAS) $(ALL_LIBS)  -lstdc++ -o uvw_fem.out $(OPT)

continuation.o: continuation.h continuation.cpp
	g++ -c  $(MY_FLAGS) -I$(DENSEMAT_DIR) continuation.h continuation.cpp -lgsl -lgslcblas $(DBG) $(OPT)

continuation: continuation.o kbhit.o
	ar rcs libcontinuation.a continuation.o kbhit.o

arnoldi_eig.o: arnoldi_eig.cpp arnoldi_eig.h
	g++ -c arnoldi_eig.h arnoldi_eig.cpp -I$(ARPACKPP_DIR)/include -I$(DENSEMAT_DIR) -fPIC

kbhit.o: kbhit.h kbhit.c
	g++ -c kbhit.h kbhit.c $(DBG) $(OPT)

uvw_model2_FV: uvw_model2_FV.h uvw_model2_FV.cpp
	g++ -c uvw_model2_FV.h uvw_model2_FV.cpp $(DBG) $(OPT)

uvw_model2_SV: uvw_model2_SV.h uvw_model2_SV.cpp
	g++ -c uvw_model2_SV.h uvw_model2_SV.cpp $(DBG) $(OPT)

fou_plot: fou_plot.cpp
	g++ fou_plot.cpp -Wl,-rpath -Wl,/usr/local/lib  $(PLPLOT_LIBS) -lgsl -lgslcblas  -o fou_plot.out $(DBG)

myfunctions: myfunction.h myfunction.cpp
	g++ -c myfunctions.h myfunctions.cpp $(DBG) $(OPT)

clean:
	rm *.o *.gch *.txt *.out

cleantxt:
	rm *.txt
