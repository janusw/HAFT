
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


%.o: %.f90
	$(FORT) $(FFLAGS) -c $<

HAFT: readHAFT2.o

clean:
	@rm -f *.o *.mod
