!.....Kinetic split fluxes
      subroutine kfvs_flux(x1, x2, qcl, qcr, qvl, qvr, resl, resr) 
      implicit none
      include 'common.h'
      real(dp) :: x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar), &
                  qvr(nvar), resl(nvar), resr(nvar)

      integer  :: i
      real(dp) :: rl, ul, vl, pl, el, rr, ur, vr, pr, er, &
                  ql2, qr2, unl, unr, vnl, vnr, &
                  lent, ct, st, Al, Bl, Ar, Br, &
                  sl, betal, sr, betar, &
                  Fp(4), Fm(4), Ff(4), dl, dr, li(4), &
                  limit, flux(nvar)

      ct   =  (x2(2) - x1(2))
      st   = -(x2(1) - x1(1))
      lent = sqrt(ct**2 + st**2)
      ct   = ct/lent
      st   = st/lent

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
      print*,'kfvs_flux: fatal error, left/right state not defined'
      stop
#endif

      ql2= ul**2 + vl**2
      el = pl/GAMMA1 + 0.5d0*rl*ql2

      qr2= ur**2 + vr**2
      er = pr/GAMMA1 + 0.5d0*rr*qr2

!     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

      vnl =-ul*st + vl*ct
      vnr =-ur*st + vr*ct

!     Positive flux
      betal = 0.5d0*rl/pl
      sl    = unl*sqrt(betal)
      Al    = 0.5d0*(1.0d0 + DERF(sl))
      Bl    = 0.5d0*exp(-sl**2)/sqrt(M_PI*betal)

      Fp(1) = rl*(unl*Al + Bl)
      Fp(2) = (pl + rl*unl**2)*Al + rl*unl*Bl
      Fp(3) = rl*(unl*Al + Bl)*vnl
      Fp(4) = (el + pl)*unl*Al + (el + 0.5d0*pl)*Bl

!     Negative flux
      betar = 0.5d0*rr/pr
      sr    = unr*sqrt(betar)
      Ar    = 0.5d0*(1.0d0 - DERF(sr))
      Br    = 0.5d0*exp(-sr**2)/sqrt(M_PI*betar)

      Fm(1) = rr*(unr*Ar - Br)
      Fm(2) = (pr + rr*unr**2)*Ar - rr*unr*Br
      Fm(3) = rr*(unr*Ar - Br)*vnr
      Fm(4) = (er + pr)*unr*Ar - (er + 0.5d0*pr)*Br

!     Total flux
      do i=1,4
         Ff(i) = Fp(i) + Fm(i)
      enddo

      flux(1) = Ff(1)*lent
      flux(2) = (ct*Ff(2) - st*Ff(3))*lent
      flux(3) = (st*Ff(2) + ct*Ff(3))*lent
      flux(4) = Ff(4)*lent

!     Total flux
      do i=1,4
         resl(i) = resl(i) + flux(i)
         resr(i) = resr(i) - flux(i)
      enddo

      end
