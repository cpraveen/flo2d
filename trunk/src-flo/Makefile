include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS} -free -Tf

#Location of forward differentiated code
SRCDIFF = .src-diff

#Tapenade command
TPND    = ${TAPENADE_HOME}/bin/tapenade -noinclude \
        -ext ${GenLib} -extAD ${ADLib} -d

TARGETS = flo2d

#Some colorful printing
RED=\033[30;43m
GREEN=\033[30;42m
BLACK=\033[0m
PRINT := $(shell \
	[ -z `echo -e` ] && PRINT='echo -e' || PRINT='echo';\
	echo $$PRINT;\
	)

#
# program dependencies
#

HDR   = $(wildcard *.h)

fvm   =   geometric.o \
          common.o \
          read_input.o \
          con2prim.o \
          main.o \
          fvresidual.o \
          average.o \
          vaverage.o \
          solid_flux.o \
          farfield_flux.o \
          kfvs_flux.o \
          lax_flux.o \
          roe_flux.o \
          hcusp_flux.o \
          limit.o \
          result.o \
			 screen.o \
          clcd.o \
          avgfact1.o \
          avgfact2.o \
          lusgs.o \
          smooth.o \
          gradient.o \
          visc.o \
          viscres.o \
          viscres_adiabatic.o \
          sutherland.o \
          read_grid.o \
          pres_match.o

impl  =   flux_jacob.o \
			 jacobian.o \
          sparsekit.o \
			 dPackgmres.o \
			 fvresidual_dq.o \
			 con2prim_dq.o \
          test_resd.o

FVRESIDUAL_DQ = fvresidual.f95 \
					 roe_flux_preproc.f95 \
					 kfvs_flux_preproc.f95 \
					 farfield_flux_preproc.f95 \
					 limit_preproc.f95 \
					 solid_flux.f95 \
					 visc.f95 \
					 viscres.f95 \
					 viscres_adiabatic.f95 \
					 gradient.f95 \
					 average.f95 \
					 vaverage.f95 \
					 sutherland.f95

gmres =   gmres.o

flo2d =   ${fvm} ${impl} ${gmres}
#
# programs
#

ALL:  $(TARGETS)

flo2d:  ${flo2d}
		${FC} -o flo2d ${flo2d} ${LIBS}
		@${RM} -f *.msg *~
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of flo2d"

#This is to prevent these files from being automatically deleted by make
.SECONDARY: ${SRCDIFF}/kfvs_flux_dq.f ${SRCDIFF}/roe_flux_dq.f

##############################################################################
# transformation rules
##############################################################################
#Flux subroutines: They must be named as somename_flux.F
#Example: roe_flux.F
%_flux.o:  %_flux.f95 $(HDR)
		${FC} -DDEBUG -DSECONDORDER -E -C -P $*_flux.f95 > preproc.f95
		${FC} -c -o $*_flux.o ${FFLAGS} preproc.f95
		${RM} -f preproc.f95

#Example: roe_flux.F
limit.o:  limit.f95 ../Makefile.in $(HDR)
		${FC} -D${TOLIMIT} -E -C -P limit.f95 > preproc.f95
		${FC} -c -o limit.o ${FFLAGS} preproc.f95
		${RM} -f preproc.f95

%_preproc.f95:  %.f95
		${FC} -DSECONDORDER -E -C -P $*.f95 > $*_preproc.f95

limit_preproc.f95:  limit.f95 ../Makefile.in $(HDR)
		${FC} -D${TOLIMIT} -E -C -P limit.f95 > limit_preproc.f95
##############################################################################
# Forward derivatives for implicit shcheme
##############################################################################
fvresidual_dq.f95: $(FVRESIDUAL_DQ) $(HDR)
	$(TPND) \
		-head fvresidual    \
		-vars "qc qv qx qy res" \
		-outvars "qv qx qy res"      \
		-difffuncname "_dq" \
		-o fvresidual       \
		$(FVRESIDUAL_DQ)

con2prim_dq.o: con2prim.f95 $(HDR)
	$(TPND) \
		-head con2prim      \
		-vars "con" \
		-outvars "prim"      \
		-difffuncname "_dq" \
		con2prim.f95
	${FC} -c ${FFLAGS} $*.f95
##############################################################################
dPackgmres.o: ../gmres/dPackgmres.f 
	${FC} -c -O3 $<

sparsekit.o: sparsekit.f
	${FC} -c -O3 $<

%.o:  %.f95 $(HDR)
		${FC} -c ${FFLAGS} $*.f95

##############################################################################
##############################################################################
# clean things up
##############################################################################

info:
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of flo2d"

clean:	
	${RM} -f *.o *_dq.f95 *_preproc.f95 preproc.f95 *.msg *~ $(TARGETS)

allclean:	
	${RM} -f *_dq.f95
