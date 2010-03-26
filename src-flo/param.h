      include 'common.h'
      include 'inf.h'
      include 'size.h'
      include 'visc.h'

      integer istart, scratch, restart
      parameter(scratch=1, restart=2)
      common/starttype/istart

! Range of bounding box
      real(dp) :: xmin, xmax, ymin, ymax
      common/range/xmin, xmax, ymin, ymax

! Type of grid
! any other value implies hybrid grid
      character gridfile*64, inpfile*32
      common/files/gridfile,inpfile

      real(dp) :: CFL, MINRES, dtglobal, gerrtol
      character(len=24) :: timemode
      integer  :: iter, ITERLAST, MAXITER, saveinterval, &
                  gmaxiter, prectype, scrinterval
      common/itparam/CFL,MINRES,dtglobal,gerrtol,iter,ITERLAST, &
                     MAXITER,saveinterval, timemode, gmaxiter, &
                     prectype, scrinterval

      real(dp) :: airk(3), birk(3)
      integer  :: NIRK
      common/timeintegration/airk,birk,NIRK

!     Number of contours
      integer  :: niso
      common/contour/niso

! Range of some variables
! rmin,rmax = density
! pmin,pmax = pressure
! mmin,mmax = mach number
      real(dp) :: rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, &
                  mmin, mmax, emin, emax, nmin, nmax
      common/minmaxprim/rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, &
                        mmin, mmax, emin, emax, nmin, nmax

      real(dp) :: fres, fres_old, fres1, fresi
      integer  :: iresi
      common/resparam/fres,fres_old,fres1,fresi,iresi

      real(dp) :: wd1(nspmax), wdmin, wdmax
      common/wdparam/wd1, wdmin, wdmax

!     Define tags for point types
      integer  :: interior, solid, farfield, outflow, boundary
      parameter(interior=0)
      parameter(solid=3)
      parameter(outflow=4)
      parameter(farfield=5)
      parameter(boundary=-1)

! Small number
      real(dp) :: EPSILON
      parameter(EPSILON=1.0d-16)

! Limiter factor for MUSCL
      integer  :: GRADTYPE, ILIMIT
      real(dp) :: LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22
      common/lim/LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22, &
                 GRADTYPE, ILIMIT

! Size of connectivity list; required by mayavi
      integer  :: lsize
      common/maya/lsize

      character(len=24) :: flux_type
      common/flux/flux_type

      integer  :: inviscid, laminar, turbulent
      parameter(inviscid=1, laminar=2, turbulent=3)

      integer  :: xvstatus, display
      common/xv/xvstatus, display

      real(dp) :: minelarea, maxelarea, mincvarea, maxcvarea
      real(dp) :: minflen, maxflen
      common/minmaxarea/minelarea, maxelarea, mincvarea, maxcvarea, &
                        minflen, maxflen
