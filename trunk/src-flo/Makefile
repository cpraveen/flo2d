include ../Makefile.in

#
# fortran flags
#
FFLAGS   = ${CFLAGS}

#Tapenade command
TPND    = ${TAPENADE_HOME}/bin/tapenade \
        -ext ${GenLib} -extAD ${ADLib} -d

TARGETS = flo

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
          vaverage_dq.o \
          killnormvel_dq.o \
          limit_dq.o \
          kfvs_flux_dq.o \
          roe_flux_dq.o  \
          solid_flux_dq.o \
          farfield_flux_dq.o \
          con2prim_dq.o \
          diag_matrix.o \
          precond.o

gmres =   gmres.o

flo   =   ${fvm} ${impl} ${gmres}
#
# programs
#

ALL:  $(TARGETS)

flo:  ${flo}
		${FC} -o flo ${flo} ${LIBS}
		@echo "==> Built ${TOLIMIT} version of flo"

#This is to prevent these files from being automatically deleted by make
.SECONDARY: kfvs_flux_dq.f roe_flux_dq.f

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

#Flux subroutines: They must be named as somename_flux.F
#Example: roe_flux.F
%_flux.o:  %_flux.F
		${GCC} -DDEBUG -DSECONDORDER -E -C -P $*_flux.F > $*_flux.f
		${FC} ${FFLAGS} -c $*_flux.f
		${RM} -f $*_flux.f

#Example: roe_flux.F
limit.o:  limit.F ../Makefile.in
		${GCC} -D${TOLIMIT} -E -C -P limit.F > limit.f
		${FC} ${FFLAGS} -c limit.f
		${RM} -f limit.f
##############################################################################
# Forward derivatives for implicit shcheme
##############################################################################
%_flux_dq.f:  %_flux.F
		${GCC} -E -C -P -DSECONDORDER $*_flux.F > $*_flux.f
		${TPND} \
		       -head       $*_flux                     \
		       -vars       "qcl qcr qvl qvr resl resr" \
		       -outvars    "resl resr"                 \
		       -difffuncname "_dq"                     \
		        $*_flux.f
		${RM} -f $*_flux.f

#-----------------------------------------------------------------------------
# Flux for solid edges
solid_flux_dq.f:	solid_flux.f
		${TPND} \
		       -head       solid_flux \
		       -vars       "qc res"   \
		       -outvars    "res"      \
		       -difffuncname "_dq"    \
		        solid_flux.f;

#-----------------------------------------------------------------------------
# Flux for farfield face
farfield_flux_dq.f: farfield_flux.F
		${GCC} -E -C -P farfield_flux.F > farfield_flux.f
		${TPND} \
		       -head       farfield_flux \
		       -vars       "qc res"      \
		       -outvars    "res"         \
		       -difffuncname "_dq"       \
		        farfield_flux.f;
		${RM} -f farfield_flux.f

#-----------------------------------------------------------------------------
# Limiter function
limit_dq.f:	limit.F ../Makefile.in
		${GCC} -D${TOLIMIT} -E -C -P limit.F > limit.f
		${TPND} \
		       -head       limit             \
		       -vars       "dl dr"           \
		       -outvars    "limit"           \
		       -difffuncname "_dq"           \
		        limit.f;
		${RM} -f limit.f

#-----------------------------------------------------------------------------
# Vertex averaging
vaverage_dq.f: vaverage.f
		${TPND} \
		       -head       vaverage         \
		       -vars       "qc qv1 qv2 qv3" \
		       -outvars    "qv1 qv2 qv3"    \
		       -difffuncname "_dq"          \
		        vaverage.f;

#-----------------------------------------------------------------------------
killnormvel_dq.f: killnormvel.f
		${TPND} \
		       -head       killnormvel \
		       -vars       "qv"        \
		       -outvars    "qv"        \
		       -difffuncname "_dq"     \
		        killnormvel.f;

#-----------------------------------------------------------------------------
con2prim_dq.f: con2prim.f
		${TPND} \
		       -head       con2prim \
		       -vars       "con"    \
		       -outvars    "prim"   \
		       -difffuncname "_dq"  \
		        con2prim.f;

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o *_bq.msg *~ $(TARGETS)

allclean:	
	${RM} -f *.o *_bq.msg *~ *_dq.f $(TARGETS)
