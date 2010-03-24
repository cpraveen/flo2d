!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------
      subroutine lusgs(elem, esue, edge, tedge, coord, qcold, qc, res, &
                       dt, tarea)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), esue(3,*), edge(2,*), tedge(2,*)
      real(dp) :: coord(2,*), qcold(nvar,*), qc(nvar,*), res(nvar,*), &
                  dt(*), tarea(*)

      integer  :: i, it, iv, p1, p2, p3, t1, t2, t3, e1, e2
      real(dp) :: u, v, a, lam, sx1, sx2, sx3, sy1, sy2, sy3, &
                  ds1, ds2, ds3, D(nt), dqc(nvar,nt), &
                  cres(nvar), flux1(nvar), flux2(nvar), &
                 lam1, lam2, lam3, omega, c1(nvar), c2(nvar)
      real(dp) :: mul1, mul2, mul, sutherland, qa(nvar), r, p, lv
      external :: sutherland

! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
      omega = 1.5d0

! Compute diagonal term
      do it=1,nt
         D(it) = tarea(it)/dt(it)
      enddo

      do i=1,nin
         e1 = edge(1,i)
         e2 = edge(2,i)
         t1 = tedge(1,i)
         t2 = tedge(2,i)

         do iv=1,4
            qa(iv) = 0.5d0*(qc(iv,t1) + qc(iv,t2))
         enddo
         r   = qa(1)
         u   = qa(2)
         v   = qa(3)
         p   = qa(4)
         a   = sqrt(GAMMA*p/r)

         sx1 =  coord(2,e2) - coord(2,e1)
         sy1 = -(coord(1,e2) - coord(1,e1))
         ds1 =  sqrt(sx1*sx1 + sy1*sy1)
         lam =  abs(u*sx1 + v*sy1) + a*ds1
         lam = 0.5d0*omega*lam
         D(t1) = D(t1) + lam
         D(t2) = D(t2) + lam

         if(iflow .ne. inviscid)then
            mul1 = sutherland(qc(1,t1))
            mul2 = sutherland(qc(1,t2))
            mul  = 0.5d0*(mul1 + mul2)
            lv   = (GAMMA/r)*(mul/prandtl)*ds1**2
            D(t1)= D(t1) + lv/tarea(t1)
            D(t2)= D(t2) + lv/tarea(t2)
         endif
      enddo

      do i=nin+1,ne
         e1 = edge(1,i)
         e2 = edge(2,i)
         t1 = tedge(1,i)

         r   = qc(1,t1)
         u   = qc(2,t1)
         v   = qc(3,t1)
         p   = qc(4,t1)
         a   = sqrt(GAMMA*p/r)

         sx1 =  coord(2,e2) - coord(2,e1)
         sy1 = -(coord(1,e2) - coord(1,e1))
         ds1 =  sqrt(sx1*sx1 + sy1*sy1)
         lam =  abs(u*sx1 + v*sy1) + a*ds1
         lam = 0.5d0*omega*lam
         D(t1) = D(t1) + lam

         if(iflow .ne. inviscid)then
            mul  = sutherland(qc(1,t1))
            lv   = (GAMMA/r)*(mul/prandtl)*ds1**2
            D(t1)= D(t1) + lv/tarea(t1)
         endif
      enddo

! Forward loop
      do it=1,nt
         p1  = elem(1,it)
         p2  = elem(2,it)
         p3  = elem(3,it)

         sx1 =  coord(2,p2) - coord(2,p1)
         sy1 = -(coord(1,p2) - coord(1,p1))
         ds1 =  sqrt(sx1*sx1 + sy1*sy1)

         sx2 =  coord(2,p3) - coord(2,p2)
         sy2 = -( coord(1,p3) - coord(1,p2) )
         ds2 =  sqrt(sx2*sx2 + sy2*sy2)

         sx3 =  coord(2,p1) - coord(2,p3)
         sy3 = -( coord(1,p1) - coord(1,p3) )
         ds3 =  sqrt(sx3*sx3 + sy3*sy3)

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         do iv=1,nvar
            cres(iv) = 0.0d0
         enddo

         if(t1 .lt. it .and. t1 .gt. 0)then
            call normalflux(sx1, sy1, qcold(1,t1),  flux1)
            call normalflux(sx1, sy1, qc(1,t1), flux2)
            call maxeig(sx1, sy1, ds1, qcold(1,t1), lam1)
            lam1 = omega*lam1
            if(iflow .ne. inviscid)then
               call visceig(it, t1, coord, elem, qcold, ds1, lam1)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) -  &
                          lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .lt. it .and. t2 .gt. 0)then
            call normalflux(sx2, sy2, qcold(1,t2),  flux1)
            call normalflux(sx2, sy2, qc(1,t2), flux2)
            call maxeig(sx2, sy2, ds2, qcold(1,t2), lam2)
            lam2 = omega*lam2
            if(iflow .ne. inviscid)then
               call visceig(it, t2, coord, elem, qcold, ds2, lam2)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                          lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .lt. it .and. t3 .gt. 0)then
            call normalflux(sx3, sy3, qcold(1,t3),  flux1)
            call normalflux(sx3, sy3, qc(1,t3), flux2)
            call maxeig(sx3, sy3, ds3, qcold(1,t3), lam3)
            lam3 = omega*lam3
            if(iflow .ne. inviscid)then
               call visceig(it, t3, coord, elem, qcold, ds3, lam3)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                          lam3*dqc(iv,t3)
            enddo
         endif

         call prim2con(qcold(1,it), c1)
         do iv=1,nvar
            dqc(iv,it) = ( -res(iv,it) - 0.5d0*cres(iv) )/D(it)
            c2(iv)     = c1(iv) + dqc(iv,it)
         enddo
         call con2prim(c2, qc(1,it))

      enddo

! Reverse loop
      do it=nt,1,-1
         p1  = elem(1,it)
         p2  = elem(2,it)
         p3  = elem(3,it)

         sx1 =   coord(2,p2) - coord(2,p1)
         sy1 = -(coord(1,p2) - coord(1,p1))
         ds1 =  sqrt(sx1*sx1 + sy1*sy1)

         sx2 =   coord(2,p3) - coord(2,p2)
         sy2 = -(coord(1,p3) - coord(1,p2))
         ds2 =  sqrt(sx2*sx2 + sy2*sy2)

         sx3 =   coord(2,p1) - coord(2,p3)
         sy3 = -(coord(1,p1) - coord(1,p3))
         ds3 =  sqrt(sx3*sx3 + sy3*sy3)

         t1 = esue(1,it)
         t2 = esue(2,it)
         t3 = esue(3,it)

         do iv=1,nvar
            cres(iv) = 0.0d0
         enddo

         if(t1 .gt. it .and. t1 .gt. 0)then
            call normalflux(sx1, sy1, qcold(1,t1),  flux1)
            call normalflux(sx1, sy1, qc(1,t1), flux2)
            call maxeig(sx1, sy1, ds1, qcold(1,t1), lam1)
            lam1 = omega*lam1
            if(iflow .ne. inviscid)then
               call visceig(it, t1, coord, elem, qcold, ds1, lam1)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                          lam1*dqc(iv,t1)
            enddo
         endif

         if(t2 .gt. it .and. t2 .gt. 0)then
            call normalflux(sx2, sy2, qcold(1,t2),  flux1)
            call normalflux(sx2, sy2, qc(1,t2), flux2)
            call maxeig(sx2, sy2, ds2, qcold(1,t2), lam2)
            lam2 = omega*lam2
            if(iflow .ne. inviscid)then
               call visceig(it, t2, coord, elem, qcold, ds2, lam2)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                          lam2*dqc(iv,t2)
            enddo
         endif

         if(t3 .gt. it .and. t3 .gt. 0)then
            call normalflux(sx3, sy3, qcold(1,t3),  flux1)
            call normalflux(sx3, sy3, qc(1,t3), flux2)
            call maxeig(sx3, sy3, ds3, qcold(1,t3), lam3)
            lam3 = omega*lam3
            if(iflow .ne. inviscid)then
               call visceig(it, t3, coord, elem, qcold, ds3, lam3)
            endif
            do iv=1,nvar
               cres(iv) = cres(iv) + (flux2(iv) - flux1(iv)) - &
                          lam3*dqc(iv,t3)
            enddo
         endif

         call prim2con(qcold(1,it), c1)
         do iv=1,nvar
            dqc(iv,it) = ( D(it)*dqc(iv,it) - 0.5d0*cres(iv) )/D(it)
            c2(iv)     = c1(iv) + dqc(iv,it)
         enddo
         call con2prim(c2, qc(1,it))

      enddo

      end
!-----------------------------------------------------------------------------
! Computes flux along (sx,sy)
!-----------------------------------------------------------------------------
      subroutine normalflux(sx, sy, qc, flux)
      implicit none
      include 'param.h'
      real(dp) :: sx, sy, qc(nvar), flux(nvar)

      real(dp) :: d, u, v, p, un, e

      d       = qc(1)
      u       = qc(2)
      v       = qc(3)
      p       = qc(4)
      e       = p/GAMMA1 + 0.5d0*d*(u*u + v*v)
      un      = u*sx + v*sy
      flux(1) = d*un
      flux(2) = p*sx + d*u*un
      flux(3) = p*sy + d*v*un
      flux(4) = (e + p)*un

      end
!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (sx, sy)
! and ds = sqrt(sx*sx + sy*sy) is face length
!-----------------------------------------------------------------------------
      subroutine maxeig(sx, sy, ds, qc, lam)
      implicit none
      include 'param.h'
      real(dp) :: sx, sy, ds, qc(nvar), lam

      real(dp) :: d, u, v, p, a

      d   = qc(1)
      u   = qc(2)
      v   = qc(3)
      p   = qc(4)
      a   = sqrt(GAMMA*p/d)
      lam = abs(u*sx + v*sy) + a*ds

      end
!-----------------------------------------------------------------------------
! Adds viscous eigenvalue
!-----------------------------------------------------------------------------
      subroutine visceig(it, t1, coord, elem, qc, ds, lam)
      implicit none
      include 'param.h'
      integer  :: it, t1, elem(3,*) 
      real(dp) :: coord(2,*), qc(nvar,*), ds, lam

      integer  :: p1, p2, p3
      real(dp) :: xt, yt, x1, y1, dr, mul, mup, lv, sutherland
      external :: sutherland

      p1 = elem(1,it)
      p2 = elem(2,it)
      p3 = elem(3,it)
      xt = (coord(1,p1) + coord(1,p2) + coord(1,p3))/3.0d0
      yt = (coord(2,p1) + coord(2,p2) + coord(2,p3))/3.0d0

      p1 = elem(1,t1)
      p2 = elem(2,t1)
      p3 = elem(3,t1)
      x1 = (coord(1,p1) + coord(1,p2) + coord(1,p3))/3.0d0
      y1 = (coord(2,p1) + coord(2,p2) + coord(2,p3))/3.0d0

      dr = sqrt( (xt-x1)**2 + (yt-y1)**2 )

      mul= sutherland(qc(1,t1))
      mup= mul/prandtl
      lv = (ds/dr)*(GAMMA/qc(1,t1))*mup
      lam= lam + lv

      end
