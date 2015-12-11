# Dependencies #

  * Tapenade
  * GMRES library
  * BLAS
  * Fortran 90 compiler

## Tapenade ##

Tapenade is an automatic differentiation tool developed at INRIA, Sophia Antipolis. Instruction on intalling tapenade are available here

http://nuwtun.berlios.de/tapenade

Download the latest version of tapenade from [ftp://ftp-sop.inria.fr/tropics/tapenade](ftp://ftp-sop.inria.fr/tropics/tapenade)

## GMRES library ##

Download the GMRES library from

http://www.cerfacs.fr/algor/Softs/GMRES/index.html

You must get the real double precision version of the library. The file required is dPackgmres.f which must be copied to the directory flo2d/gmres.

# Compilation #

Set some variables like the fortran compiler in flo2d/Makefile.in file. Then go into src-flo directory and type make.