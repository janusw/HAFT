
# default compiler
FORT=gfortran
FORT_NOPATH = $(notdir $(FORT))

# flags for different compilers
# GCC
ifeq ($(findstring gfortran,$(FORT_NOPATH)),gfortran)
	FFLAGS=-std=f2003 -Wall -Wextra
endif
# INTEL
ifeq ($(FORT_NOPATH),ifort)
	FFLAGS=-warn
endif
# SUN/ORACLE
ifeq ($(findstring sunf95,$(FORT_NOPATH)),sunf95)
	FFLAGS=-g -w2
endif
# PGI
ifeq ($(FORT_NOPATH),pgf95)
	FFLAGS=-g -Mstandard -Mextend -Minform=inform
endif


# default target: build the HAFT code with a Fortran compiler
HAFT: readHAFT2.o

%.o: %.f90
	$(FORT) $(FFLAGS) -c $<

# target 'f2py': build the HAFT Python bindings via f2py
f2py: HAFT.so

HAFT.so: readHAFT2.f90
	f2py -m HAFT -c $< skip: momspread readhaftmatrix getmatrixval gettableval param smearhades4momentum readhaftpairmatrix

all: HAFT f2py

clean:
	@rm -f *.o *.so *.mod
