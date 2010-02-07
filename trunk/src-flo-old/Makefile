include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS}

#Location of forward differentiated code
SRCDIFF = .src-diff

#Tapenade command
TPND    = ${TAPENADE_HOME}/bin/tapenade -noinclude \
        -ext ${GenLib} -extAD ${ADLib} -d -O ${SRCDIFF}

TARGETS = flo

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

fvm   =   geometric.o \
          common.o \
          con2prim.o \
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
			 screen.o \
          clcd.o \
          avgfact1.o \
          avgfact2.o \
          lusgs.o \
          smooth.o \
          gradient.o \
          visc.o \
          sutherland.o

impl  =   fvresidual_q.o \
          average_q.o \
          ${SRCDIFF}/vaverage_dq.o \
          ${SRCDIFF}/killnormvel_dq.o \
          ${SRCDIFF}/limit_dq.o \
          ${SRCDIFF}/kfvs_flux_dq.o \
          ${SRCDIFF}/roe_flux_dq.o  \
          ${SRCDIFF}/solid_flux_dq.o \
          ${SRCDIFF}/farfield_flux_dq.o \
          ${SRCDIFF}/con2prim_dq.o \
			 flux_jacob.o \
			 jacobian.o \
          sparsekit.o \
			 dPackgmres.o \
          test_resd.o

gmres =   gmres.o

flo   =   ${fvm} ${impl} ${gmres}
#
# programs
#

ALL:  $(TARGETS)

flo:  ${flo}
		${FC} -o flo ${flo} ${LIBS}
		@${RM} -f *.msg *~
		@$(PRINT) "==> Built $(GREEN)${TOLIMIT}$(BLACK) version of flo"

#This is to prevent these files from being automatically deleted by make
.SECONDARY: ${SRCDIFF}/kfvs_flux_dq.f ${SRCDIFF}/roe_flux_dq.f

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

${SRCDIFF}/%.o: ${SRCDIFF}/%.f
	cd ${SRCDIFF}; ${FC} ${FFLAGS} -c $*.f

##############################################################################
# Print warning if forward diff code is older than original code
##############################################################################
average_q.f: average.f
fvresidual_q.f: fvresidual.f
		@$(PRINT) "$(RED) WARNING: $@ is older than $< $(BLACK)"

#Flux subroutines: They must be named as somename_flux.F
#Example: roe_flux.F
%_flux.o:  %_flux.F
		${FC} -DDEBUG -DSECONDORDER -E -C -P $*_flux.F > preproc.f
		${FC} ${FFLAGS} -c -o $*_flux.o preproc.f
		${RM} -f preproc.f

#Example: roe_flux.F
limit.o:  limit.F ../Makefile.in
		${FC} -D${TOLIMIT} -E -C -P limit.F > preproc.f
		${FC} ${FFLAGS} -c -o limit.o preproc.f
		${RM} -f preproc.f
##############################################################################
# Forward derivatives for implicit shcheme
##############################################################################
${SRCDIFF}/%_flux_dq.f:  %_flux.F
		${FC} -E -C -P -DSECONDORDER $*_flux.F > preproc.f
		${TPND} \
		       -head       $*_flux                     \
		       -vars       "qcl qcr qvl qvr resl resr" \
		       -outvars    "resl resr"                 \
		       -difffuncname "_dq"                     \
		        preproc.f
		${RM} -f preproc.f

#-----------------------------------------------------------------------------
# Flux for solid edges
${SRCDIFF}/solid_flux_dq.f:	solid_flux.f
		${TPND} \
		       -head       solid_flux \
		       -vars       "qc res"   \
		       -outvars    "res"      \
		       -difffuncname "_dq"    \
		        solid_flux.f;

#-----------------------------------------------------------------------------
# Flux for farfield face
${SRCDIFF}/farfield_flux_dq.f: farfield_flux.F
		${FC} -E -C -P farfield_flux.F > preproc.f
		${TPND} \
		       -head       farfield_flux \
		       -vars       "qc res"      \
		       -outvars    "res"         \
		       -difffuncname "_dq"       \
		        preproc.f
		${RM} -f preproc.f

#-----------------------------------------------------------------------------
# Limiter function
${SRCDIFF}/limit_dq.f:	limit.F ../Makefile.in
		${FC} -D${TOLIMIT} -E -C -P limit.F > preproc.f
		${TPND} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "limit"           \
		       -difffuncname "_dq"           \
		        preproc.f
		${RM} -f preproc.f

#-----------------------------------------------------------------------------
# Vertex averaging
${SRCDIFF}/vaverage_dq.f: vaverage.f
		${TPND} \
		       -head       vaverage         \
		       -vars       "qc qv1 qv2 qv3" \
		       -outvars    "qv1 qv2 qv3"    \
		       -difffuncname "_dq"          \
		        vaverage.f;

#-----------------------------------------------------------------------------
${SRCDIFF}/killnormvel_dq.f: killnormvel.f
		${TPND} \
		       -head       killnormvel \
		       -vars       "qv"        \
		       -outvars    "qv"        \
		       -difffuncname "_dq"     \
		        killnormvel.f;

#-----------------------------------------------------------------------------
${SRCDIFF}/con2prim_dq.f: con2prim.f
		${TPND} \
		       -head       con2prim \
		       -vars       "con"    \
		       -outvars    "prim"   \
		       -difffuncname "_dq"  \
		        con2prim.f;

##############################################################################
dPackgmres.o: ../gmres/dPackgmres.f 
	${FC} ${FFLAGS} -c $<

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *.msg *~ $(TARGETS)

allclean:	
	${RM} -f *.o $(TARGETS) $(SRCDIFF)/*
