C------------------------------------------------------------------------------
C     Computes roe dissipation matrix
C------------------------------------------------------------------------------
      subroutine flux_jacob(x1, x2, qc, jac, mode)
      implicit none
      include 'common.h'
      integer          mode
      double precision x1(2), x2(2), qc(*), jac(nvar,nvar)

      integer          i, j, k
      double precision 
     &                 ua, va, pa, qa2, aa2, aa, ha, ra,
     &                 una, vna, ct, st, lent,
     &                 m2, t1, t2, t3, t4, t5, l1, l2, l3, l4,
     &                 S(nvar,nvar), R(nvar,nvar)

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)
      ct = ct/lent
      st = st/lent

C     State
      ra = qc(1)
      ua = qc(2)
      va = qc(3)
      pa = qc(4)
      qa2= ua**2 + va**2
      aa2= GAMMA*pa/ra
      aa = dsqrt(aa2)
      ha = aa2/GAMMA1 + 0.5d0*qa2

C     Rotated velocity
      una = ua*ct + va*st
      vna =-ua*st + va*ct

C     Eigenvalues with entropy fix
      if(mode.eq.-1)then
         l1 = dmin1(una, 0.0d0)
         l2 = l1
         l3 = dmin1(una + aa, 0.0d0)
         l4 = dmin1(una - aa, 0.0d0)
      else
         l1 = dmax1(una, 0.0d0)
         l2 = l1
         l3 = dmax1(una + aa, 0.0d0)
         l4 = dmax1(una - aa, 0.0d0)
      endif

c     Right eigenvector matrix
      t1 = 0.5d0*ra/aa
      m2 = qa2/aa2

      R(1,1) = 1.0d0
      R(2,1) = ua
      R(3,1) = va
      R(4,1) = 0.5d0*qa2

      R(1,2) = 0.0d0
      R(2,2) = ra*st
      R(3,2) = -ra*ct
      R(4,2) = -ra*vna

      R(1,3) = t1
      R(2,3) = t1*(ua + aa*ct)
      R(3,3) = t1*(va + aa*st)
      R(4,3) = t1*(ha + aa*una)

      R(1,4) = t1
      R(2,4) = t1*(ua - aa*ct)
      R(3,4) = t1*(va - aa*st)
      R(4,4) = t1*(ha - aa*una)

c     Inverse of right eigenvector matrix
      t1     = 0.5d0*gamma1*m2*aa/ra
      t2     = una/ra
      t3     = gamma1*ua/aa
      t4     = gamma1*va/aa
      t5     = gamma1/(ra*aa)

      S(1,1) = 1.0d0 - 0.5d0*gamma1*m2
      S(2,1) = vna/ra
      S(3,1) = t1 - t2
      S(4,1) = t1 + t2

      S(1,2) = t3/aa
      S(2,2) = st/ra
      S(3,2) = (ct - t3)/ra
      S(4,2) =-(ct + t3)/ra

      S(1,3) = t4/aa
      S(2,3) = -ct/ra
      S(3,3) = (st - t4)/ra
      S(4,3) =-(st + t4)/ra

      S(1,4) = -gamma1/aa2
      S(2,4) = 0.0d0
      S(3,4) = t5
      S(4,4) = t5

c     Multiply lambda * S
      do i=1,nvar
         R(i,1) = l1*R(i,1)
         R(i,2) = l2*R(i,2)
         R(i,3) = l3*R(i,3)
         R(i,4) = l4*R(i,4)
      enddo

c     Now multiply with S
      do i=1,nvar
         do j=1,nvar
            jac(i,j) = 0.0d0
            do k=1,nvar
               jac(i,j) = jac(i,j) + R(i,k)*S(k,j)
            enddo
            jac(i,j) = lent*jac(i,j)
         enddo
      enddo

      return
      end

C------------------------------------------------------------------------------
C     Computes roe dissipation matrix
C------------------------------------------------------------------------------
      subroutine solid_flux_jacob(x1, x2, qc, jac)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qc(*), jac(nvar,nvar)

      double precision u, v, q2, ct, st

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))

C     State
      u  = qc(2)
      v  = qc(3)
      q2 = u**2 + v**2

      jac(2,1) = jac(2,1) + 0.5d0*gamma1*q2*ct
      jac(2,2) = jac(2,2) - gamma1*u*ct
      jac(2,3) = jac(2,3) - gamma1*v*ct
      jac(2,4) = jac(2,4) - gamma1*ct

      jac(3,1) = jac(3,1) + 0.5d0*gamma1*q2*st
      jac(3,2) = jac(3,2) - gamma1*u*st
      jac(3,3) = jac(3,3) - gamma1*v*st
      jac(3,4) = jac(3,4) - gamma1*st

      return
      end
