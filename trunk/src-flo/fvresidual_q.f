C------------------------------------------------------------------------------
C Computes: finite volume residual jacobian * given vector qcd
C------------------------------------------------------------------------------
      subroutine fvresidual_q(elem, edge, tedge, vedge, spts, bdedge,
     +                       coord, qc, qv, qx, qy, af, tarea, cl, cd, 
     +                       qcd1, resd)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 vedge(2,nemax), spts(nspmax), bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), tarea(ntmax), resd(nvar,ntmax),
     +                 qcd1(nvar,ntmax)
      double precision cl, cd
      double precision qx(3,npmax), qy(3,npmax)

      integer          i, j, ie, v1, v2, e1, e2, c1, c2
      double precision resl(nvar), resr(nvar)
      double precision qcd(nvar,ntmax), qpdummy(nvar), qvd(nvar,npmax)
      double precision qcon(nvar)

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

C     Transform conserved to primitive variables
      do i=1,nt
         call prim2con(qc(1,i), qcon)
         call con2prim_dq(qcon, qcd1(1,i), qpdummy, qcd(1,i))
      enddo

C     Compute area averaged value at vertices
      call average_q(spts, elem, edge, bdedge, coord, tarea, af, qc,
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

