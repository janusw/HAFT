# HAFT
The HADES Acceptance Filter for Theorists (HAFT)

This repository contains a Fortran2003 version of the HAFT code, based on the original F77 code from https://hades-wiki.gsi.de/cgi-bin/view/Public/HadesHAFT. This version is being used in the transport codes GiBUU and SMASH (in the latter via f2py).

In contrast to the original code, the reading of the acceptance matrices is done in a compiler-indepedent fashion. This version of the code has been verified to compile with gfortran, ifort, sunf95 and pgf95 and should work with any compiler that supports the Fortran2003 standard.

There are two modules, 'HAFT_single' and 'HAFT_pair'.

The module 'HAFT_single' provides routines for single-particle acceptance filtering. Usage:

1. set acceptance file name with
        `call setFileName(fname)`
2. sample single acceptance values with calls to
        `acc = getHadesAcceptance(id,mom,theta,phi,mode)`
3. apply detector resolution with calls to
        `call smearhadesmomentum(...)`

The module 'HAFT_pair' provides routines for pair acceptance filtering. Usage:

1. set acceptance file name with
        `call setPairFileName(fname)`
2. sample pair acceptance values with calls to
        `acc = getHadesPairAcceptance(mass,pt,rapidity,mode)`
3. set resolution parameters via
        `call setResolutionParameters(...)`
4. apply detector resolution via
        `call smearHadesPair(...)`
