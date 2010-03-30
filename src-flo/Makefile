include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS} -free -Tf

#Location of forward differentiated code
SRCDIFF = .src-diff

#Tapenade command
TPN     = ${TAPENADE_HOME}/bin/tapenade -noinclude \
        -ext ${GenLib} -extAD ${ADLib}
TPNDQ   = $(TPN) -d -difffuncname "_dq"
TPNBQ   = $(TPN) -b -difffuncname "_bq"

TARGETS = flo2d adj2d

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
          sutherland.o \
          read_grid.o \
          pres_match.o

impl  =   flux_jacob.o \
			 jacobian.o \
          sparsekit.o \
			 dPackgmres.o \
			 fvresidual_dq.o \
			 con2prim_dq.o \
			 gmres.o \
          test_resd.o

FVRESIDUAL_Q =  fvresidual.f95 \
					 roe_flux_preproc.f95 \
					 kfvs_flux_preproc.f95 \
					 farfield_flux_preproc.f95 \
					 limit_preproc.f95 \
					 solid_flux.f95 \
					 visc.f95 \
					 gradient.f95 \
					 average.f95 \
					 vaverage.f95 \
					 sutherland.f95

COSTFUN_Q = costfun.f95 \
            clcd.f95 \
				average.f95 \
				vaverage.f95 \
				gradient.f95 \

adj = adj2d.o \
		sutherland.o \
		common.o \
		geometric.o \
		limit.o \
		gradient.o \
		average.o \
		avgfact1.o \
		avgfact2.o \
		clcd.o \
		read_grid.o \
		read_input.o \
		vaverage.o \
		write_vtk.o \
      fvresidual_bq.o \
      costfun_bq.o \
      con2prim_bq.o \
		clcd_bq.o \
      adStack.o \
      adBuffer.o

flo2d =   ${fvm} ${impl}

adj2d =   ${adj}
#
# programs
#

ALL:  $(TARGETS)

flo2d:  ${flo2d}
		${FC} -o flo2d ${flo2d} ${LIBS}
		@${RM} -f *.msg *~
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of flo2d"

adj2d:  ${adj2d}
		${FC} -o adj2d ${adj2d} ${LIBS}
		@${RM} -f *.msg *~
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of adj2d"

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
fvresidual_dq.f95: $(FVRESIDUAL_Q) $(HDR)
	$(TPNDQ) \
		-head    fvresidual    \
		-vars    "qc qv qx qy res" \
		-outvars "qv qx qy res"      \
		-o       fvresidual       \
		$(FVRESIDUAL_Q)

con2prim_dq.o: con2prim.f95 $(HDR)
	$(TPNDQ) \
		-head    con2prim      \
		-vars    "con" \
		-outvars "prim"      \
		con2prim.f95
	${FC} -c ${FFLAGS} $*.f95
##############################################################################
fvresidual_bq.f95: $(FVRESIDUAL_Q) $(HDR)
	$(TPNBQ) \
		-head    fvresidual    \
		-vars    "qc qv qx qy res" \
		-outvars "qc qv qx qy res"      \
		-o       fvresidual       \
		$(FVRESIDUAL_Q)
##############################################################################
costfun_bq.f95: $(COSTFUN_Q) $(HDR)
	$(TPNBQ) \
		-head    costfun            \
		-vars    "qc qv qx qy cost" \
		-outvars "qc qv qx qy cost" \
		$(COSTFUN_Q)
##############################################################################
con2prim_bq.o: con2prim.f95 $(HDR)
	$(TPNBQ) \
		-head    con2prim  \
		-vars    "con"     \
		-outvars "prim"    \
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
adBuffer.o: $(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
	$(FC) -c $<

adStack.o: $(TAPENADE_HOME)/ADFirstAidKit/adStack.c
	$(CC) -c $<

##############################################################################
# clean things up
##############################################################################

info:
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of flo2d"

clean:	
	${RM} -f *.o 
	${RM} -f *_dq.f95
	${RM} -f *_bq.f95
	${RM} -f *_preproc.f95
	${RM} -f preproc.f95
	${RM} -f *.msg *~ $(TARGETS)

allclean:	
	${RM} -f *_dq.f95
