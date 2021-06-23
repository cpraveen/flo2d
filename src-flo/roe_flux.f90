! Roe flux function
      subroutine roe_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar), &
                  qvr(nvar), resl(nvar), resr(nvar)

      integer  :: i
      real(dp) :: rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr, &
                  ua, va, qa2, aa2, aa, ha, &
                  ql2, qr2, rl12, rr12, rd, &
                  unl, unr, una, vna, ct, st, Fc(4), Fd(4), &
                  m1, m2, a1, a2, a3, a4, l1, l2, l3, l4, &
                  a1l1, a2l2, a3l3, a4l4, aact, aast, &
                  du1, du2, du3, du4, flux, dl, dr, li(4), limit, &
                  lent, e1, e4, del
      real(dp), parameter :: ETOL=0.01d0

      ct =  (x2(2) - x1(2))
      st = -(x2(1) - x1(1))
      lent = sqrt(ct**2 + st**2)
      ct = ct/lent
      st = st/lent

#if defined SECONDORDER
      do i=1,4
         dl    = qcl(i) - qvl(i)
         dr    = qvr(i) - qcr(i)
         li(i) = LIMIT(dl, dr)
      enddo

!     Left state
      rl = qcl(1) + 0.5d0*li(1)
      ul = qcl(2) + 0.5d0*li(2)
      vl = qcl(3) + 0.5d0*li(3)
      pl = qcl(4) + 0.5d0*li(4)

!     Right state
      rr = qcr(1) - 0.5d0*li(1)
      ur = qcr(2) - 0.5d0*li(2)
      vr = qcr(3) - 0.5d0*li(3)
      pr = qcr(4) - 0.5d0*li(4)

#elif defined FIRSTORDER
!     Left state
      rl = qcl(1)
      ul = qcl(2)
      vl = qcl(3)
      pl = qcl(4)

!     Right state
      rr = qcr(1)
      ur = qcr(2)
      vr = qcr(3)
      pr = qcr(4)
#else
      print*,'roe_flux: fatal error, left/right state not defined'
      stop
#endif

      ql2= ul**2 + vl**2
      al2= GAMMA*pl/rl
      hl = al2/GAMMA1 + 0.5d0*ql2

      qr2= ur**2 + vr**2
      ar2= GAMMA*pr/rr
      hr = ar2/GAMMA1 + 0.5d0*qr2

!     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

!     Centered flux
      Fc(1) = rl*unl            + rr*unr
      Fc(2) = pl*ct + rl*ul*unl + pr*ct + rr*ur*unr
      Fc(3) = pl*st + rl*vl*unl + pr*st + rr*vr*unr
      Fc(4) = rl*hl*unl         + rr*hr*unr

!     Roe average
      rl12 = sqrt(rl)
      rr12 = sqrt(rr)
      rd   = 1.0d0/(rl12 + rr12)

      ua   = (ul*rl12 + ur*rr12)*rd
      va   = (vl*rl12 + vr*rr12)*rd
      ha   = (hl*rl12 + hr*rr12)*rd
      qa2  = ua**2 + va**2
      aa2  = GAMMA1*(ha - 0.5d0*qa2)

#ifdef DEBUG
      if(aa2 .le. 0.0d0)then
         print*,'Sonic speed is negative'
         print*,'Left/right cell values'
         print*,qcl(1),qcl(2),qcl(3),qcl(4)
         print*,qcr(1),qcr(2),qcr(3),qcr(4)
         print*,'Left/right vertex values'
         print*,qvl(1),qvl(2),qvl(3),qvl(4)
         print*,qvr(1),qvr(2),qvr(3),qvr(4)
         print*
         print*,rl,ul,vl,pl
         print*,rr,ur,vr,pr
         print*,li
         stop
      endif
#endif
      aa  = sqrt(aa2)
      una = ua*ct + va*st
      vna =-ua*st + va*ct

!     Eigenvalues with entropy fix
      e1 = abs(una - aa)
      l2 = abs(una)
      l3 = l2
      e4 = abs(una + aa)

      del= ETOL*aa
      if(e1 .lt. del)then
         l1 = 0.5d0*(del + e1**2/del)
      else
         l1 = e1
      endif

      if(e4 .lt. del)then
         l4 = 0.5d0*(del + e4**2/del)
      else
         l4 = e4
      endif

!     Difference of conserved variables
      du1 = rr           - rl
      du2 = rr*ur        - rl*ul
      du3 = rr*vr        - rl*vl
      du4 = (rr*hr - pr) - (rl*hl - pl)

!     Amplitudes
      m1 = (ct*du2 + st*du3 - una*du1)/aa
      m2 = GAMMA1*(du4 - ua*du2 - va*du3 + qa2*du1)/aa**2

      a4 = 0.5d0*(m1 + m2)
      a1 = 0.5d0*(m2 - m1)
      a3 = du1 - a1 - a4
      a2 = ( st*du2 - ct*du3 + vna*du1 )/aa

!     Diffusive flux
      a1l1  = a1*l1
      a2l2  = a2*l2
      a3l3  = a3*l3
      a4l4  = a4*l4
      aact  = aa*ct
      aast  = aa*st

      Fd(1) = a1l1               +               a3l3           + a4l4
      Fd(2) = a1l1*(ua - aact)   + a2l2*aa*st  + a3l3*ua        + &
              a4l4*(ua + aact)
      Fd(3) = a1l1*(va - aast)   - a2l2*aa*ct  + a3l3*va        + &
              a4l4*(va + aast)
      Fd(4) = a1l1*(ha - una*aa) + a2l2*aa*vna + a3l3*0.5d0*qa2 + &
              a4l4*(ha + una*aa)

!     Total flux
      do i=1,4
         flux    = 0.5d0*lent*( Fc(i) - Fd(i) )
         resl(i) = resl(i) + flux
         resr(i) = resr(i) - flux
      enddo

      end
