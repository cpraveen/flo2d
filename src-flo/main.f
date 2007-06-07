C------------------------------------------------------------------------------
C Main program of vertex-centroid scheme similar to Jameson, Frink
C------------------------------------------------------------------------------
      program main
      implicit none
      include 'param.h'
      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 esue(3,ntmax), vedge(2,nemax), spts(nspmax),
     +                 bdedge(2,nbpmax), esubp(mesubp,nbpmax),
     +                 ptype(npmax)
      double precision coord(2,npmax), qc(nvar,ntmax), 
     +                 qcold(nvar,ntmax), dt(ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax),
     +                 drmin(ntmax), res(nvar,ntmax), c1(nvar),
     +                 c2(nvar), c3(nvar)
      double precision cl, cd
      double precision qx(3,npmax), qy(3,npmax), qcd(nvar,ntmax)
      real etime, elapsed(2), totaltime

      integer          i, j, irk

      call math
      call read_input
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, carea, af)

C Set initial condition
      call initialize(qc, cl, cd)

C For testing AD gradients
c     call test_resd(elem, edge, tedge, vedge, spts, bdedge,
c    +                 coord, qc, qv, qx, qy, af, carea, dt, cl, cd,
c    +                 res, qcd)
c     call write_result(coord, elem, edge, qc, qv, cl, cd)
c     stop

c     call time_step2(edge, tedge, carea, coord, qc, dt)
c     call jacobian(elem, esue, coord, carea, dt, qc)
c     stop

      iter = 0
      fres = 1.0d20
      call system('rm -f FLO.RES')
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
c        call time_step(drmin, qc, dt)
         call time_step2(edge, tedge, carea, coord, qc, dt)
         call save_old(qc, qcold)

         do irk=1,nirk

C           Compute finite volume residual
            call fvresidual(elem, edge, tedge, vedge, spts, bdedge,
     +                      coord, qc, qv, qx, qy, af, carea, cl, cd, 
     +                      res)

C           Update the solution
            if(timemode .eq. 1)then
               do i=1,nt
                  call prim2con(qcold(1,i), c1)
                  call prim2con(qc(1,i),    c2)
                  do j=1,nvar
                     c3(j) = airk(irk)*c1(j) + 
     +                       birk(irk)*(c2(j)-dt(i)*res(j,i)/carea(i))
                  enddo
                  call con2prim(c3, qc(1,i))
               enddo
            elseif(timemode .eq. 2)then
               call lusgs(elem, esue, edge, tedge, coord, qcold, qc, 
     +                    res, dt, carea)
            elseif(timemode .eq. 3)then
               call gmres(elem, esue, edge, tedge, vedge, spts, bdedge,
     +                    coord, qc, qv, qx, qy, af, carea, dt, cl, cd,
     +                    res, qcd)
               do i=1,nt
                  call prim2con(qc(1,i),    c2)
                  do j=1,nvar
                     c3(j) = c2(j) + qcd(j,i)
                  enddo
                  call con2prim(c3, qc(1,i))
               enddo
            endif

         enddo

         iter = iter + 1
         call residue(res)
         call clcd(edge, tedge, coord, qc, cl, cd)
         open(unit=99, file='FLO.RES', access='append')
         write(99,'(i6,4e16.6)') iter, fres, fresi, cl, cd
         close(99)
         if(mod(iter,saveinterval) .eq. 0)then
            call write_result(coord, elem, edge, qc, qv, cl, cd)
         endif

         if(timemode.eq.3)then
            cfl = dmax1(1.0d0, 1.0d0/fres)
c           cfl = -2.0d0 + 3.0d0*iter
c           cfl = 1.0d0*iter
            cfl = dmin1(cfl,20.0d0)
         endif

      enddo

      call write_result(coord, elem, edge, qc, qv, cl, cd)
      call write_sol(iter, fres, cl, cd)

      totaltime = etime(elapsed)
      totaltime = totaltime/60.0
      elapsed(1)= elapsed(1)/60.0
      elapsed(2)= elapsed(1)/60.0
      print *, 'Time: total=', totaltime, ' user=', elapsed(1),
     +         ' system=', elapsed(2)

      stop
      end
