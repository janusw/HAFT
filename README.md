# HAFT
The HADES Acceptance Filter for Theorists (HAFT)

This repository contains a Fortran2003 version of the HAFT code, based on the original F77 code from https://hades-wiki.gsi.de/cgi-bin/view/Public/HadesHAFT.

In contrast to the original code, the reading of the acceptance matrices is done in a compiler-indepedent fashion. This version of the code has been verified to compile with gfortran, ifort, sunf95 and pgf95.


Usage :

1. set acceptance file name with
       `call setFileName(fname)`
2. sample single acceptance values with calls to
       `acc = getHadesAcceptance(id,mom,theta,phi,mode)`
3. sample pair acceptance values with calls to
       `acc = getHadesPairAcceptance(mass,pt,rapidity,mode)`
