CC=g++
MPICC=mpicc
COPT=-O4
CFLAGS=
LIBS=   -llapack -lblas -lm 
LIBMPI= -lmpi
DIRBLACS=/usr/lib64
LIBBLACS=   $(DIRBLACS)/libmpiblacsCinit.a  $(DIRBLACS)/libmpiblacsF77init.a    $(DIRBLACS)/libmpiblacs.a
LIBSCALAPACK= -lscalapack 

principal: principal.cpp
	$(CC) $(COPT) $(CFLAGS) principal.cpp -o principal $(LIBS)

clean:
	rm principal
