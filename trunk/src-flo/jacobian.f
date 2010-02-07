C------------------------------------------------------------------------------
C Compute and assemble first order jacobian in CSR format
C computes by looping over triangles
C------------------------------------------------------------------------------
      subroutine jacobian(elem, esue, coord, tarea, dt, qc)
      implicit none
      include 'param.h'
      include 'gmres.h'
      integer          elem(3,*), esue(3,*)
      double precision coord(2,*), tarea(*), dt(*), qc(nvar,*)

      integer          it, i, j, p1, p2, p3, t1, t2, t3, irow, 
     +                 icol, jcount, rcount, iw(nvar*nt), ierr, iwk
      double precision j0(nvar,nvar), j1(nvar,nvar), j2(nvar,nvar),
     +                 j3(nvar,nvar), fact

c     Count number of elements in jmat
      jcount = 0

c     Count number of rows
      rcount = 0

c     Loop over triangles
      do it=1,nt

         fact = tarea(it)/dt(it)

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
            call flux_jacob_diag(coord(1,p1),coord(1,p2),
     1                           qc(1,it),qc(1,t1),j0)
            call flux_jacob_lax(coord(1,p1),coord(1,p2),
     1                          qc(1,it),qc(1,t1),j1,-1)
         elseif(t1.eq.-solid)then
c           call solid_flux_jacob(coord(1,p1),coord(1,p2),qc(1,it),j0)
            call flux_jacob_diag(coord(1,p1),coord(1,p2),
     1                           qc(1,it),qc(1,it),j0)
         elseif(t1.eq.-farfield)then
            call flux_jacob(coord(1,p1),coord(1,p2),qc(1,it),j0,+1)
         endif

         if(t2.gt.0)then
            call flux_jacob_diag(coord(1,p2),coord(1,p3),
     1                           qc(1,it),qc(1,t2),j0)
            call flux_jacob_lax(coord(1,p2),coord(1,p3),
     1                          qc(1,it),qc(1,t2),j2,-1)
         elseif(t2.eq.-solid)then
c           call solid_flux_jacob(coord(1,p2),coord(1,p3),qc(1,it),j0)
            call flux_jacob_diag(coord(1,p2),coord(1,p3),
     1                           qc(1,it),qc(1,it),j0)
         elseif(t2.eq.-farfield)then
            call flux_jacob(coord(1,p2),coord(1,p3),qc(1,it),j0,+1)
         endif

         if(t3.gt.0)then
            call flux_jacob_diag(coord(1,p3),coord(1,p1),
     1                           qc(1,it),qc(1,t3),j0)
            call flux_jacob_lax(coord(1,p3),coord(1,p1),
     1                          qc(1,it),qc(1,t3),j3,-1)
         elseif(t3.eq.-solid)then
c           call solid_flux_jacob(coord(1,p3),coord(1,p1),qc(1,it),j0)
            call flux_jacob_diag(coord(1,p3),coord(1,p1),
     1                           qc(1,it),qc(1,it),j0)
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
c     call csrcsc(4*nt, 1, 1, jmat, ja, ia, alu, jlu, ju)
c     call csrcsc(4*nt, 1, 1, alu, jlu, ju, jmat, ja, ia)

c     Compute ILU(0)
c     call ilu0(nvar*nt, jmat, ja, ia, alu, jlu, ju, iw, ierr)

      iwk = 2*(4*nvar*nvar*nt) + (nvar*nt)
      call iluk(nvar*nt, jmat, ja, ia, lfil, alu, jlu, ju, levs, iwk,
     1          w, jw, ierr)
      if(ierr.gt.0)then
         print*,'iluk: Zero pivot encountered at step number',ierr
         stop
      elseif(ierr.eq.-1)then
         print*,'iluk: Error. input matrix may be wrong'
         stop
      elseif(ierr.eq.-2)then
         print*,'iluk: The matrix L overflows the array alu'
         stop
      elseif(ierr.eq.-3)then
         print*,'iluk: The matrix U overflows the array alu.'
         stop
      elseif(ierr.eq.-4)then
         print*,'iluk: Illegal value for lfil'
         stop
      elseif(ierr.eq.-5)then
         print*,'iluk: zero row encountered in A or U'
         stop
      endif

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

C------------------------------------------------------------------------------
      subroutine flux_jacob_diag(x1, x2, qc1, qc2, jac)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qc1(*), qc2(*), jac(nvar,nvar)

      integer          i
      double precision r1, u1, v1, p1, a1, r2, u2, v2, p2, a2, 
     &                 ua, va, aa,
     &                 una, ct, st, lent, lam

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)

C     State
      r1 = qc1(1)
      u1 = qc1(2)
      v1 = qc1(3)
      p1 = qc1(4)
      a1 = dsqrt(GAMMA*p1/r1)

      r2 = qc2(1)
      u2 = qc2(2)
      v2 = qc2(3)
      p2 = qc2(4)
      a2 = dsqrt(GAMMA*p2/r2)

      ua = 0.5d0*(u1 + u2)
      va = 0.5d0*(v1 + v2)
      aa = 0.5d0*(a1 + a2)

C     Rotated velocity
      una = ua*ct + va*st

c     Maximum eigenvalue
      lam = dabs(una) + aa * lent

      do i=1,nvar
         jac(i,i) = jac(i,i) + 0.5d0*lam
      enddo

      return
      end

C------------------------------------------------------------------------------
      subroutine flux_jacob_lax(x1, x2, qc1, qc2, jac, mode)
      implicit none
      include 'common.h'
      integer          mode
      double precision x1(2), x2(2), qc1(*), qc2(*), jac(nvar,nvar)

      integer          i, j
      double precision r1, u1, v1, p1, a1, r2, u2, v2, p2, a2, 
     &                 q22, h2, un2, ua, va, aa,
     &                 una, ct, st, lent, lam, amat(nvar,nvar)

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)

C     State
      r1 = qc1(1)
      u1 = qc1(2)
      v1 = qc1(3)
      p1 = qc1(4)
      a1 = dsqrt(GAMMA*p1/r1)

      r2 = qc2(1)
      u2 = qc2(2)
      v2 = qc2(3)
      q22= u2**2 + v2**2
      p2 = qc2(4)
      a2 = dsqrt(GAMMA*p2/r2)
      h2 = a2**2/(GAMMA-1.0d0) + 0.5d0*q22

      ua = 0.5d0*(u1 + u2)
      va = 0.5d0*(v1 + v2)
      aa = 0.5d0*(a1 + a2)

C     Rotated velocity
      una = ua*ct + va*st

c     Maximum eigenvalue
      lam = dabs(una) + aa * lent

      un2 = u2*ct + v2*st

      amat(1,1) = 0.0d0
      amat(2,1) = -u2*un2 + 0.5d0*(GAMMA-1.0d0)*q22*ct
      amat(3,1) = -v2*un2 + 0.5d0*(GAMMA-1.0d0)*q22*st
      amat(4,1) = -un2*(h2 - 0.5d0*(GAMMA-1.0d0)*q22)

      amat(1,2) = ct
      amat(2,2) = (3.0d0-GAMMA)*u2*ct + v2*st
      amat(3,2) = -(GAMMA-1.0d0)*u2*st + v2*ct
      amat(4,2) = h2*ct - (GAMMA-1.0d0)*u2*un2

      amat(1,3) = st
      amat(2,3) = -(GAMMA-1.0d0)*v2*ct + u2*st
      amat(3,3) = (3.0d0-GAMMA)*v2*st + u2*ct
      amat(4,3) = h2*st - (GAMMA-1.0d0)*v2*un2

      amat(1,4) = 0.0d0
      amat(2,4) = (GAMMA-1.0d0)*ct
      amat(3,4) = (GAMMA-1.0d0)*st
      amat(4,4) = GAMMA*un2

      do i=1,nvar
         do j=1,nvar
            jac(i,j) = 0.5d0*amat(i,j)
         enddo
         jac(i,i) = jac(i,i) - 0.5d0*lam
      enddo

      return
      end

