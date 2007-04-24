      subroutine jacobian(elem, esue, coord, qc, jmat, ja, ia)
      implicit none
      include 'param.h'
      integer          elem(3,*), esue(3,*), ja(*), ia(*)
      double precision coord(2,*), qc(nvar,*), jmat(*)

      integer          it, i, j, p1, p2, p3, t1, t2, t3, nzrow, irow
      double precision j0(nvar,nvar), j1(nvar,nvar), j2(nvar,nvar),
     +                 j3(nvar,nvar)

      count = 0

      do it=1,nt

         do i=1,nvar
            do j=1,nvar
               j0(i,j) = 0.0d0
               j1(i,j) = 0.0d0
               j2(i,j) = 0.0d0
               j3(i,j) = 0.0d0
            enddo
         enddo

c        vertices of the triangle
         p1 = elem(1,it)
         p2 = elem(2,it)
         p3 = elem(3,it)

c        neighbouring triangles
         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         if(t1.gt.0)then
            call flux_jacob(coord(1,p1),coord(1,p2),qc(1,it),j0,+1)
            call flux_jacob(coord(1,p1),coord(1,p2),qc(1,t1),j1,-1)
         elseif(t1.eq.-solid)then
            call solid_flux_jacob(coord(1,p1),coord(1,p2),qc(1,it),j0)
         elseif(t1.eq.-farfield)then
            call flux_jacob(coord(1,p1),coord(1,p2),qc(1,it),j0,+1)
         endif

         if(t2.gt.0)then
            call flux_jacob(coord(1,p2),coord(1,p3),qc(1,it),j0,+1)
            call flux_jacob(coord(1,p2),coord(1,p3),qc(1,t2),j2,-1)
         elseif(t2.eq.-solid)then
            call solid_flux_jacob(coord(1,p2),coord(1,p3),qc(1,it),j0)
         elseif(t2.eq.-farfield)then
            call flux_jacob(coord(1,p2),coord(1,p3),qc(1,it),j0,+1)
         endif

         if(t3.gt.0)then
            call flux_jacob(coord(1,p3),coord(1,p1),qc(1,it),j0,+1)
            call flux_jacob(coord(1,p3),coord(1,p1),qc(1,t3),j3,-1)
         elseif(t3.eq.-solid)then
            call solid_flux_jacob(coord(1,p3),coord(1,p1),qc(1,it),j0)
         elseif(t3.eq.-farfield)then
            call flux_jacob(coord(1,p3),coord(1,p1),qc(1,it),j0,+1)
         endif

c        number of non-zero elements in each row
         nzrow = nvar
         if(t1.gt.0) nzrow = nzrow + nvar
         if(t2.gt.0) nzrow = nzrow + nvar
         if(t3.gt.0) nzrow = nzrow + nvar

c        Put jacobian in CSR format
         do irow=1,nvar

            do icol=1,nvar
               count       = count + 1
               jmat(count) = j0(irow,icol)
               ja(count)   = (it-1)*nvar + icol
            enddo

            if(t1.gt.0)then
               do icol=1,nvar
                  count       = count + 1
                  jmat(count) = j1(irow,icol)
                  ja(count)   = (t1-1)*nvar + icol
               enddo
            endif

            if(t2.gt.0)then
               do icol=1,nvar
                  count       = count + 1
                  jmat(count) = j2(irow,icol)
                  ja(count)   = (t2-1)*nvar + icol
               enddo
            endif

            if(t3.gt.0)then
               do icol=1,nvar
                  count       = count + 1
                  jmat(count) = j3(irow,icol)
                  ja(count)   = (t3-1)*nvar + icol
               enddo
            endif

            ia(rcount) = 1
         enddo

      enddo

      return
      end
