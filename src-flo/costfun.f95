!------------------------------------------------------------------------------
! Objective or cost function
!------------------------------------------------------------------------------
      subroutine costfun(elem, edge, tedge, vedge, spts, &
                         coord, qc, qv, qx, qy, af, tarea, varea, &
                         cl, cd, cost)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*), tedge(2,*), vedge(2,*),  &
                  spts(*)
      real(dp) :: coord(2,*), qc(nvar,*), af(3,*), qv(nvar,*),  &
                  tarea(*), varea(*)
      real(dp) :: cl, cd, cost
      real(dp) :: qx(3,*), qy(3,*)

      integer  :: i, j, ie, v1, v2, e1, e2, c1, c2

!     Compute area averaged value at vertices
      call average(spts, elem, edge, coord, tarea, af, qc, qv)

!     Viscous terms
      if(iflow .ne. inviscid)then
         call gradient(elem, edge, spts, coord, varea, qc, qv, &
                       qx, qy)
      endif

      call clcd(edge, tedge, coord, qc, cl, cd)

      cost = cl

      end
