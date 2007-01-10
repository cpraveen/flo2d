C------------------------------------------------------------------------------
C Computes finite volume residual
C------------------------------------------------------------------------------
      subroutine fvresidual(elem, edge, tedge, vedge, spts, bdedge,
     +                      coord, qc, qv, qx, qy, af, carea, cl, cd, 
     +                      res)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 vedge(2,nemax), spts(nspmax), bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax), res(nvar,ntmax)
      double precision cl, cd
      double precision qx(3,npmax), qy(3,npmax)

      integer          i, j, ie, v1, v2, e1, e2, c1, c2

      do i=1,nt
         do j=1,nvar
            res(j,i) = 0.0d0
         enddo
      enddo

C     Compute area averaged value at vertices
      call average(spts, elem, edge, bdedge, coord, carea, af, qc,
     +             qv)

C     Compute flux for interior edges
      if(iflux .eq. ikfvs) goto 10
      if(iflux .eq. iroe ) goto 20
      print*,'flux is not defined'
      stop

C     KFVS flux
10    do ie=nin1,nin2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         c2 = tedge(2,ie)
         v1 = vedge(1,ie)
         v2 = vedge(2,ie)
         call kfvs_flux(coord(1,e1), coord(1,e2),
     +                 qc(1,c1), qc(1,c2), qv(1,v1), qv(1,v2), 
     +                 res(1,c1), res(1,c2))
      enddo
      goto 100

C     Roe flux
20    do ie=nin1,nin2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         c2 = tedge(2,ie)
         v1 = vedge(1,ie)
         v2 = vedge(2,ie)
         call roe_flux(coord(1,e1), coord(1,e2),
     +                qc(1,c1), qc(1,c2), qv(1,v1), qv(1,v2), 
     +                res(1,c1), res(1,c2))
      enddo

C Compute flux for solid wall edges
100   do ie=nsw1,nsw2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call solid_flux(coord(1,e1), coord(1,e2),
     +                   qc(1,c1), res(1,c1))       
      enddo

C Flux for far-field points
      do ie=nff1,nff2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         call farfield_flux(coord(1,e1), coord(1,e2), qc(1,c1),
     +                      cl, cd, res(1,c1))
      enddo

C Viscous terms
      if(iflow .ne. inviscid)then
         call gradient(elem, edge, bdedge, spts, coord, qc, qv, 
     +                 qx, qy)
         call viscflux(edge, tedge, coord, qv, qx, qy, res)
      endif

      return
      end
