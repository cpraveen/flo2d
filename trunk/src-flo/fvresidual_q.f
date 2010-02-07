C------------------------------------------------------------------------------
C Computes: finite volume residual jacobian * given vector qcd
C------------------------------------------------------------------------------
      subroutine fvresidual_q(elem, edge, tedge, vedge, spts, bdedge,
     +                       coord, qc, qv, qx, qy, af, tarea, cl, cd, 
     +                       qcd1, resd)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*), tedge(2,*), vedge(2,*), spts(*),
     1            bdedge(2,*)
      real(dp) :: coord(2,*), qc(nvar,*), af(3,*), qv(nvar,*), tarea(*),
     1            resd(nvar,*), qcd1(nvar,*)
      real(dp) :: cl, cd
      real(dp) :: qx(3,*), qy(3,*)

      integer  :: i, j, ie, v1, v2, e1, e2, c1, c2
      real(dp) :: resl(nvar), resr(nvar)
      real(dp) :: qcd(nvar,nt), qpdummy(nvar), qvd(nvar,np)
      real(dp) :: qcon(nvar)

C     Initialize to zero
      do i=1,nt
         do j=1,nvar
            resd(j,i) = 0.0d0
         enddo
      enddo

C     Dummy residue variables
      do j=1,nvar
         resl(j) = 0.0d0
         resr(j) = 0.0d0
      enddo

C     Transform conserved qcd1 to primitive variables qcd
      do i=1,nt
         call prim2con(qc(1,i), qcon)
         call con2prim_dq(qcon, qcd1(1,i), qpdummy, qcd(1,i))
      enddo

C     Compute area averaged value at vertices
      call average_q(spts, elem, edge, coord, tarea, af, qc,
     +               qcd, qv, qvd)

C     Compute flux for interior edges
      select case(iflux)

      case(ikfvs)
         do ie=nin1,nin2
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            c1 = tedge(1,ie)
            c2 = tedge(2,ie)
            v1 = vedge(1,ie)
            v2 = vedge(2,ie)
            call kfvs_flux_dq(coord(1,e1), coord(1,e2),
     +                        qc(1,c1), qcd(1,c1), qc(1,c2), qcd(1,c2),
c    +                        qv(1,v1), qv(1,v2),
     +                        qv(1,v1), qvd(1,v1), qv(1,v2), qvd(1,v2),
     +                        resl, resd(1,c1), resr, resd(1,c2))
         enddo

      case(iroe)
         do ie=nin1,nin2
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            c1 = tedge(1,ie)
            c2 = tedge(2,ie)
            v1 = vedge(1,ie)
            v2 = vedge(2,ie)
            call  roe_flux_dq(coord(1,e1), coord(1,e2),
     +                        qc(1,c1), qcd(1,c1), qc(1,c2), qcd(1,c2),
c    +                        qv(1,v1), qv(1,v2),
     +                        qv(1,v1), qvd(1,v1), qv(1,v2), qvd(1,v2),
     +                        resl, resd(1,c1), resr, resd(1,c2))
         enddo

      case default
         print*,'fvresidual: Unknown flux type ',iflux
         stop

      end select

C     Compute flux for solid wall edges
      do ie=nsw1,nsw2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call solid_flux_dq(coord(1,e1), coord(1,e2),
     +                      qc(1,c1), qcd(1,c1), resl, resd(1,c1))
      enddo

C     Flux for far-field points
      do ie=nff1,nff2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call farfield_flux_dq(coord(1,e1), coord(1,e2), qc(1,c1),
     +                         qcd(1,c1), cl, cd, resl, resd(1,c1))
      enddo

C     TO BE CONTINUED
C     Viscous terms
      if(iflow .ne. inviscid)then
c        call gradient(elem, edge, bdedge, spts, coord, qc, qv, 
c    +                 qx, qy)
c        call viscflux(edge, tedge, coord, qv, qx, qy, res)
         print*,'fvresidual_q: viscous not implemented yet'
         stop
      endif

      return
      end

