include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS}

TARGETS = flo

#
# program dependencies
#

flo   =   geometric.o \
          common.o \
          main.o \
			 fvresidual.o \
          average.o \
          vaverage.o \
	       killnormvel.o \
	       solid_flux.o \
	       farfield_flux.o \
	       kfvs_flux.o \
	       lax_flux.o \
	       roe_flux.o \
	       hcusp_flux.o \
	       limit.o \
	       result.o \
	       clcd.o \
	       avgfact1.o \
	       avgfact2.o \
	       lusgs.o \
	       smooth.o \
			 gradient.o \
			 visc.o \
			 sutherland.o

#
# programs
#

ALL:  $(TARGETS)

flo:  ${flo}
		${FC} -o flo ${flo}
		@echo "==> Built ${TOLIMIT} version of flo"

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

#Flux subroutines: They must be named as somename_flux.F
#Example: roe_flux.F
%_flux.o:  %_flux.F
		${GCC} -DDEBUG -E -C -P $*_flux.F > $*_flux.f
		${FC} ${FFLAGS} -c $*_flux.f
		${RM} -f $*_flux.f

#Example: roe_flux.F
limit.o:  limit.F ../Makefile.in
		${GCC} -D${TOLIMIT} -E -C -P limit.F > limit.f
		${FC} ${FFLAGS} -c limit.f
		${RM} -f limit.f
##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o $(TARGETS)
