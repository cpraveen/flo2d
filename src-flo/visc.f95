!------------------------------------------------------------------------------
! Add viscous fluxes to the FV residual
!------------------------------------------------------------------------------
      subroutine viscflux(edge, tedge, coord, qv, qx, qy, res)
      implicit none
      include 'param.h'
      integer  :: edge(2,*), tedge(2,*)
      real(dp) :: coord(2,*), qv(nvar,*), qx(3,*), &
                  qy(3,*), res(nvar,*)

      integer  :: i, e1, e2, c1, c2
      real(dp) :: resdummy(nvar)

! Viscous flux for interior edges
      do i=1,nin
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         c2 = tedge(2,i)
         call viscres(coord(1,e1), coord(1,e2), qv(1,e1), qv(1,e2), &
                      qx(1,e1), qx(1,e2), qy(1,e1), qy(1,e2), &
                      res(1,c1), res(1,c2))
      enddo

! Viscous flux for boundary edges
      do i=1,nvar
         resdummy(i) = 0.0d0
      enddo

!     Solid wall edges
      do i=nsw1,nsw2
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call viscres_adiabatic &
                     (coord(1,e1), coord(1,e2), qv(1,e1), qv(1,e2), &
                      qx(1,e1), qx(1,e2), qy(1,e1), qy(1,e2), &
                      res(1,c1), resdummy)
      enddo

!     Farfield edges
      do i=nff1,nff2
         e1 = edge(1,i)
         e2 = edge(2,i)
         c1 = tedge(1,i)
         call viscres(coord(1,e1), coord(1,e2), qv(1,e1), qv(1,e2), &
                      qx(1,e1), qx(1,e2), qy(1,e1), qy(1,e2), &
                      res(1,c1), resdummy)
      enddo

      end
