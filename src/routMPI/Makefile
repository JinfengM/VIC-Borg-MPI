########################################################################
### rout.f makefile #################################################### 
########################################################################
#
# Routing algorithm written D. Lohmann
#
# This is a slightly modified code (main algotrithms unchanged -IO and
# array dimensions simplified).
# Maintained by G. O'Donnell (tempgd@hydro.washington.edu) and Andy Wood
# 
# $Id: Makefile,v 1.1 2005/04/07 05:07:28 vicadmin Exp $
#

#This program uses the non-standard Fortran argument GETARG
#Different compilers require different flags to link with this function
#Comment out one of the following depending on your compiler

#If compiling on SUN and LINUX use  (remember -O)
FFLAGS = -fPIC -ffixed-line-length-none
#If compiling on HP use
#FFLAGS = -C -O +U77 -ffixed-line-length-none
#for debugging
#FFLAGS = -C -g -lm -ffixed-line-length-none

FC=gfortran

HFILES=    	parameter.h

OBJECTS=	routmpi.o 			\
		make_convolution.o	\
                init_routines.o		\
                read_routines.o		\
                write_routines.o	\
                unit_hyd_routines.o

exe:		$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -shared -o routMPI.so

routmpi.o:			routmpi.f
make_convolution.o:		make_convolution.f
init_routines.o:		init_routines.f
read_routines.o:		read_routines.f
write_routines.o:		write_routines.f
unit_hyd_routines.o:		unit_hyd_routines.f

clean:
	/bin/rm *.o
