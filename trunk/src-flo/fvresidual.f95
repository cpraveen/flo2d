!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
      subroutine fvresidual(elem, edge, tedge, vedge, spts, &
                            coord, qc, qv, qx, qy, af, tarea, varea, &
                            cl, cd, res)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*), tedge(2,*), vedge(2,*),  &
                  spts(*)
      real(dp) :: coord(2,*), qc(nvar,*), af(3,*), qv(nvar,*),  &
                  tarea(*), varea(*), res(nvar,*)
      real(dp) :: cl, cd
      real(dp) :: qx(3,*), qy(3,*)

      integer  :: i, j, ie, v1, v2, e1, e2, c1, c2

      do i=1,nt
         do j=1,nvar
            res(j,i) = 0.0d0
         enddo
      enddo

!     Compute area averaged value at vertices
      call average(spts, elem, edge, coord, tarea, af, qc, qv)

!     Compute flux for interior edges
      select case(iflux)

      case(ikfvs)
         do ie=nin1,nin2
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            c1 = tedge(1,ie)
            c2 = tedge(2,ie)
            v1 = vedge(1,ie)
            v2 = vedge(2,ie)
            call kfvs_flux(coord(1,e1), coord(1,e2), &
                           qc(1,c1), qc(1,c2), qv(1,v1), qv(1,v2), &
                           res(1,c1), res(1,c2))
            enddo

      case(iroe)
         do ie=nin1,nin2
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            c1 = tedge(1,ie)
            c2 = tedge(2,ie)
            v1 = vedge(1,ie)
            v2 = vedge(2,ie)
            call roe_flux(coord(1,e1), coord(1,e2), &
                          qc(1,c1), qc(1,c2), qv(1,v1), qv(1,v2),  &
                          res(1,c1), res(1,c2))
         enddo

      case default
         print*,'fvresidual: Unknown flux type ',iflux
         stop

      end select

!     Compute flux for solid wall edges
      do ie=nsw1,nsw2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call solid_flux(coord(1,e1), coord(1,e2), qc(1,c1), res(1,c1))
      enddo

!     Flux for far-field points
      do ie=nff1,nff2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call farfield_flux(coord(1,e1), coord(1,e2), qc(1,c1), &
                            cl, cd, res(1,c1))
      enddo

!     Viscous terms
      if(iflow .ne. inviscid)then
         call gradient(elem, edge, spts, coord, varea, qc, qv, &
                       qx, qy)
         call viscflux(edge, tedge, coord, qv, qx, qy, res)
      endif

      end
