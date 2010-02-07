C------------------------------------------------------------------------------
C Compute vertex values from cell-center values using my new averaging
C procedure which has linear consistency
C------------------------------------------------------------------------------
      subroutine average(spts, elem, edge, bdedge, coord, tarea, af, qc,
     +                   qv)
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), tarea(ntmax)

      integer          i, j, v1, v2, v3, e1, e2

      do i=1,np
         do j=1,nvar
            qv(j,i) = 0.0d0
         enddo
      enddo
      do i=1,nt
         v1 = elem(1,i)
         v2 = elem(2,i)
         v3 = elem(3,i)
         call vaverage(coord(1,v1), coord(1,v2), coord(1,v3),
     +                 af(1,v1), af(1,v2), af(1,v3), tarea(i), 
     +                 qc(1,i), qv(1,v1), qv(1,v2), qv(1,v3))
      enddo
      do i=1,np
         do j=1,nvar
            qv(j,i) = qv(j,i)/af(3,i)
         enddo
      enddo

      if(iflow .ne. inviscid)then
C        Viscous flow: zero velocity condition
         do i=1,nsp
            j = spts(i)
            qv(2,j) = 0.0d0
            qv(3,j) = 0.0d0
         enddo
      endif

      return
      end

