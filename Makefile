
FORT=gfortran
FFLAGS=-std=f2003 -Wall -Wextra

%.o:	%.f90
	$(FORT) $(FFLAGS) -c $<

HAFT: readHAFT2.o

clean:
	@rm -f *.o *.mod
