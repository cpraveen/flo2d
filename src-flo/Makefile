include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS}

#Location of forward differentiated code
SRCDIFF = .src-diff

#Tapenade command
TPN     = ${TAPENADE_HOME}/bin/tapenade -noinclude \
			 -ext ${GenLib}
TPNDQ   = $(TPN) -d -tgtfuncname "_dq"
TPNBQ   = $(TPN) -b -adjfuncname "_bq"

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

FVRESIDUAL_Q =  fvresidual.f90 \
					 roe_flux_preproc.f90 \
					 kfvs_flux_preproc.f90 \
					 farfield_flux_preproc.f90 \
					 limit_preproc.f90 \
					 solid_flux.f90 \
					 visc.f90 \
					 gradient.f90 \
					 average.f90 \
					 vaverage.f90 \
					 sutherland.f90

COSTFUN_Q = costfun.f90 \
            clcd.f90 \
				average.f90 \
				vaverage.f90 \
				gradient.f90 \

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
      agmres.o \
      jacobian.o \
      flux_jacob.o \
      sparsekit.o \
      dPackgmres.o \
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
%_flux.o:  %_flux.f90 $(HDR)
		${FC} -DDEBUG -DSECONDORDER -E -cpp -C -P $*_flux.f90 > preproc.f90
		${FC} -c -o $*_flux.o ${FFLAGS} preproc.f90
		${RM} -f preproc.f90

#Example: roe_flux.F
limit.o:  limit.f90 ../Makefile.in $(HDR)
		${FC} -D${TOLIMIT} -E -cpp -C -P limit.f90 > preproc.f90
		${FC} -c -o limit.o ${FFLAGS} preproc.f90
		${RM} -f preproc.f90

%_preproc.f90:  %.f90
		${FC} -DSECONDORDER -E -cpp -C -P $*.f90 > $*_preproc.f90

limit_preproc.f90:  limit.f90 ../Makefile.in $(HDR)
		${FC} -D${TOLIMIT} -E -cpp -C -P limit.f90 > limit_preproc.f90
##############################################################################
# Forward derivatives for implicit shcheme
##############################################################################
fvresidual_dq.f90: $(FVRESIDUAL_Q) $(HDR)
	$(TPNDQ) \
		-head    "fvresidual(qc qv qx qy res)/(qv qx qy res)"   \
		-o       fvresidual       \
		$(FVRESIDUAL_Q)

con2prim_dq.o: con2prim.f90 $(HDR)
	$(TPNDQ) \
		-head  "con2prim(con prim)/(con prim)" \
		con2prim.f90
	#mv con2prim_dq.f90 con2prim_dq.f90
	${FC} -c ${FFLAGS} $*.f90
##############################################################################
fvresidual_bq.f90: $(FVRESIDUAL_Q) $(HDR)
	$(TPNBQ) \
		-head    "fvresidual(qc qv qx qy res)/(qc qv qx qy res)" \
		-o       fvresidual       \
		$(FVRESIDUAL_Q)
##############################################################################
costfun_bq.f90: $(COSTFUN_Q) $(HDR)
	$(TPNBQ) \
		-head    "costfun(qc qv qx qy cost)/(qc qv qx qy cost)" \
		$(COSTFUN_Q)
##############################################################################
con2prim_bq.o: con2prim.f90 $(HDR)
	$(TPNBQ) \
		-head    "con2prim(con prim)/(con prim)"    \
		con2prim.f90
	${FC} -c ${FFLAGS} $*.f90
##############################################################################
dPackgmres.o: ../gmres/dPackgmres.f 
	${FC} -c -O3 $<

sparsekit.o: sparsekit.f
	${FC} -c -O3 $<

%.o:  %.f90 $(HDR)
		${FC} -c ${FFLAGS} $*.f90

%_bq.o:  %_bq.f90 $(HDR)
		${FC} -c ${FFLAGS} $*_bq.f90

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
	${RM} -f *_dq.f90
	${RM} -f *_bq.f90
	${RM} -f *_preproc.f90
	${RM} -f preproc.f90
	${RM} -f *.msg *~ $(TARGETS)

allclean:	
	${RM} -f *_dq.f90
