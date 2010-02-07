C------------------------------------------------------------------------------
C Main program of vertex-centroid scheme similar to Jameson, Frink
C------------------------------------------------------------------------------
      program main
      implicit none
      include 'param.h'
      include 'gmres.h'
      integer, allocatable, dimension(:)   :: ptype
      integer, allocatable, dimension(:,:) :: elem, esue, edge, tedge,
     1                                        vedge
      real(dp),allocatable, dimension(:)   :: dt, drmin, tarea
      real(dp),allocatable, dimension(:,:) :: coord, qc, qcold, af,
     1                                        qv, res, qx, qy, qcd

      integer  ::  spts(nspmax), bdedge(2,nbpmax)
      real(dp) :: c1(nvar), c2(nvar), c3(nvar)
      real(dp) :: cl, cd
      real     :: etime, elapsed(2), totaltime

      integer  :: i, j, irk

c     Set some constants
      call math

c     Read input parameters from file, also sets np, nt
      call read_input

c     Allocate memory for variables
      include 'allocate.h'

c     Geometry preprocessor
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, coord, drmin, tarea, af)

C Set initial condition
      call initialize(qc, cl, cd)

C For testing AD gradients
c     call test_resd(elem, edge, tedge, vedge, spts, bdedge,
c    +                 coord, qc, qv, qx, qy, af, tarea, dt, cl, cd,
c    +                 res, qcd)
c     call write_result(coord, elem, edge, qc, qv, cl, cd)
c     stop

c     call time_step2(edge, tedge, tarea, coord, qc, dt)
c     call jacobian(elem, esue, coord, tarea, dt, qc)
c     stop

      iter = 0
      fres = 1.0d20
      call system('rm -f FLO.RES')
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
c        call time_step(drmin, qc, dt)
         call time_step2(edge, tedge, tarea, coord, qc, dt)
         call save_old(qc, qcold)

         do irk=1,nirk

C           Compute finite volume residual
            call fvresidual(elem, edge, tedge, vedge, spts, bdedge,
     +                      coord, qc, qv, qx, qy, af, tarea, cl, cd, 
     +                      res)

C           Update the solution
            if(timemode .eq. 1)then
               do i=1,nt
                  call prim2con(qcold(1,i), c1)
                  call prim2con(qc(1,i),    c2)
                  do j=1,nvar
                     c3(j) = airk(irk)*c1(j) + 
     +                       birk(irk)*(c2(j)-dt(i)*res(j,i)/tarea(i))
                  enddo
                  call con2prim(c3, qc(1,i))
               enddo
            elseif(timemode .eq. 2)then
               call lusgs(elem, esue, edge, tedge, coord, qcold, qc, 
     +                    res, dt, tarea)
            elseif(timemode .eq. 3)then
               call gmres(elem, esue, edge, tedge, vedge, spts, bdedge,
     +                    coord, qc, qv, qx, qy, af, tarea, dt, cl, cd,
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
            call average(spts, elem, edge, coord, tarea, af, qc, qv)
            call write_result(coord, elem, edge, qc, qv, cl, cd)
         endif
         if(mod(iter,scrinterval) .eq. 0)then
            call screen(qc, cl, cd)
         endif

         if(timemode.eq.3)then
            cfl = 1.0e20                ! for subsonic
c           if(iter.gt.20) cfl = 1.0e20 ! for transonic
c           cfl = dmax1(1.0d0, 10.0d0/fres)
c           cfl = -2.0d0 + 3.0d0*iter
         endif

      enddo

c     Save final solution
      call average(spts, elem, edge, coord, tarea, af, qc, qv)
      call write_result(coord, elem, edge, qc, qv, cl, cd)
      call write_sol(iter, fres, cl, cd)

c     Pressure matching objective function
      open(25,file='p.dat')
      write(25,*)(qv(4,i),i=1,np)
      close(25)
c     call pres_match(edge, coord, qv)

      totaltime = etime(elapsed)
      totaltime = totaltime/60.0
      elapsed(1)= elapsed(1)/60.0
      elapsed(2)= elapsed(1)/60.0
      print *, 'Time: total=', totaltime, ' user=', elapsed(1),
     +         ' system=', elapsed(2)

      stop
      end
