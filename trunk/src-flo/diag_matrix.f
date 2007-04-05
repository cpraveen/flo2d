C------------------------------------------------------------------------------
C Computes block diagonal elements using roe scheme
C------------------------------------------------------------------------------
      subroutine diag_matrix(edge, tedge, coord, qc, carea, dt)
      implicit none
      include 'param.h'
      include 'data.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision coord(2,npmax), qc(nvar,ntmax),
     +                 carea(ntmax), dt(ntmax)

      integer          i, j, k, ie, e1, e2, c1, c2
      double precision dm(nvar,nvar), f1, f2
      integer          ipvt(nvar), info
      double precision work(nvar), det

C     Initialize to identity
      do i=1,nt
         do j=1,nvar
            do k=1,nvar
               dmat(j,k,i) = 0.0d0
            enddo
            dmat(j,j,i) = 1.0d0
         enddo
      enddo

C     Roe flux
      do ie=nin1,nin2
         e1 = edge(1,ie)
         e2 = edge(2,ie)
         c1 = tedge(1,ie)
         c2 = tedge(2,ie)
         call roe_dmat(coord(1,e1), coord(1,e2), qc(1,c1), qc(1,c2), dm)
         f1 = 0.5d0*dt(c1)/carea(c1)
         f2 = 0.5d0*dt(c2)/carea(c2)
         do i=1,nvar
            do j=1,nvar
               dmat(i,j,c1) = dmat(i,j,c1) + f1*dm(i,j)
               dmat(i,j,c2) = dmat(i,j,c2) + f2*dm(i,j)
            enddo
         enddo
      enddo

C     Compute inverse of block diagonal matrix
      do i=1,nt
c        write(*,'(4e14.4)')((dmat(j,k,i),k=1,4),j=1,4)
         call dgefa(dmat(1,1,i),nvar,nvar,ipvt,info)
         call dgedi(dmat(1,1,i),nvar,nvar,ipvt,det,work,01)
c        write(*,'(4e14.4)')((dmat(j,k,i),k=1,4),j=1,4)
c        stop
      enddo

      return
      end

C------------------------------------------------------------------------------
C     Computes roe dissipation matrix
C------------------------------------------------------------------------------
      subroutine roe_dmat(x1, x2, qcl, qcr, dissm)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), qcl(nvar), qcr(nvar), 
     +                 dissm(nvar,nvar)

      integer          i, j, k
      double precision rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr,
     &                 ua, va, qa2, aa2, aa, ha, ra,
     &                 ql2, qr2, rl12, rr12, rd,
     &                 unl, unr, una, vna, ct, st, lent,
     &                 m2, t1, t2, t3, t4, t5, l1, l2, l3, l4,
     &                 S(nvar,nvar), R(nvar,nvar)

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)
      ct = ct/lent
      st = st/lent

C     Left state
      rl = qcl(1)
      ul = qcl(2)
      vl = qcl(3)
      pl = qcl(4)

C     Right state
      rr = qcr(1)
      ur = qcr(2)
      vr = qcr(3)
      pr = qcr(4)

      ql2= ul**2 + vl**2
      al2= GAMMA*pl/rl
      hl = al2/GAMMA1 + 0.5d0*ql2

      qr2= ur**2 + vr**2
      ar2= GAMMA*pr/rr
      hr = ar2/GAMMA1 + 0.5d0*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

C     Roe average
      rl12 = dsqrt(rl)
      rr12 = dsqrt(rr)
      rd   = 1.0d0/(rl12 + rr12)

      ra   = dsqrt(rl*rr)
      ua   = (ul*rl12 + ur*rr12)*rd
      va   = (vl*rl12 + vr*rr12)*rd
      ha   = (hl*rl12 + hr*rr12)*rd
      qa2  = ua**2 + va**2
      aa2  = GAMMA1*(ha - 0.5d0*qa2)
      aa  = dsqrt(aa2)
      una = ua*ct + va*st
      vna =-ua*st + va*ct

C     Eigenvalues with entropy fix
      l1 = dabs(una)
      l2 = l1
      l3 = dabs(una + aa)
      l4 = dabs(una - aa)

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
            dissm(i,j) = 0.0d0
            do k=1,nvar
               dissm(i,j) = dissm(i,j) + R(i,k)*S(k,j)
            enddo
            dissm(i,j) = lent*dissm(i,j)
         enddo
      enddo

      return
      end
