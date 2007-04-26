C------------------------------------------------------------------------------
C Compute and assemble first order jacobian in CSR format
C computes by looping over triangles
C------------------------------------------------------------------------------
      subroutine jacobian(elem, esue, coord, carea, dt, qc)
      implicit none
      include 'param.h'
      include 'data.h'
      integer          elem(3,*), esue(3,*)
      double precision coord(2,*), carea(*), dt(*), qc(nvar,*)

      integer          it, i, j, p1, p2, p3, t1, t2, t3, irow, 
     +                 icol, jcount, rcount, iw(nvar*ntmax), ierr
      double precision j0(nvar,nvar), j1(nvar,nvar), j2(nvar,nvar),
     +                 j3(nvar,nvar), fact

c     Count number of elements in jmat
      jcount = 0

c     Count number of rows
      rcount = 0

c     Loop over triangles
      do it=1,nt

         fact = carea(it)/dt(it)

         do i=1,nvar
            do j=1,nvar
               j0(i,j) = 0.0d0
               j1(i,j) = 0.0d0
               j2(i,j) = 0.0d0
               j3(i,j) = 0.0d0
            enddo
            j0(i,i)    = fact
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

c        Put jacobian in CSR format
         do irow=1,nvar
            rcount     = rcount + 1
            ia(rcount) = jcount + 1

            do icol=1,nvar
               jcount       = jcount + 1
               jmat(jcount) = j0(irow,icol)
               ja(jcount)   = (it-1)*nvar + icol
            enddo

            if(t1.gt.0)then
               do icol=1,nvar
                  jcount       = jcount + 1
                  jmat(jcount) = j1(irow,icol)
                  ja(jcount)   = (t1-1)*nvar + icol
               enddo
            endif

            if(t2.gt.0)then
               do icol=1,nvar
                  jcount       = jcount + 1
                  jmat(jcount) = j2(irow,icol)
                  ja(jcount)   = (t2-1)*nvar + icol
               enddo
            endif

            if(t3.gt.0)then
               do icol=1,nvar
                  jcount       = jcount + 1
                  jmat(jcount) = j3(irow,icol)
                  ja(jcount)   = (t3-1)*nvar + icol
               enddo
            endif

         enddo

      enddo

      ia(rcount+1) = jcount + 1

c     Perform two transpositions to order the columns
      call csrcsc(4*nt, 1, 1, jmat, ja, ia, alu, jlu, ju)
      call csrcsc(4*nt, 1, 1, alu, jlu, ju, jmat, ja, ia)

c     Compute ILU(0)
      call ilu0(nvar*nt, jmat, ja, ia, alu, jlu, ju, iw, ierr)

c     do irow=1,1
c     t1=ia(irow)
c     t2=ia(irow+1)-1
c     print*,'begin,end=',t1,t2
c     do i=t1,t2
c     print*,ja(i),jmat(i)
c     enddo
c     enddo

c     print*,'jcount =',jcount
c     print*,'rcount =',rcount

      return
      end
