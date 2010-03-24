!------------------------------------------------------------------------------
! Main program of vertex-centroid scheme similar to Jameson, Frink
!------------------------------------------------------------------------------
      program main
      implicit none
      include 'param.h'
      include 'gmres.h'
      integer, allocatable, dimension(:)   :: ptype
      integer, allocatable, dimension(:,:) :: elem, esue, edge, tedge, &
                                              vedge
      real(dp),allocatable, dimension(:)   :: dt, drmin, tarea, varea
      real(dp),allocatable, dimension(:,:) :: coord, qc, qcold, af, &
                                              qv, res, qx, qy, qcd, &
                                              qvd, qxd, qyd
      ! Arrays for gmres                                        
      integer,target,allocatable,dimension(:)  :: a_ja, a_ia, a_jlu, &
                                                  a_ju, a_levs, a_jw
      real(dp),target,allocatable,dimension(:) :: a_jmat, a_alu, a_w

      integer  ::  spts(nspmax)
      real(dp) :: c1(nvar), c2(nvar), c3(nvar)
      real(dp) :: cl, cd
      real     :: etime, elapsed(2), totaltime

      integer  :: i, j, irk

!     Set some constants
      call math

!     Read input parameters from file, also sets np, nt
      call read_input

!     Allocate memory for variables
      include 'allocate.h'

!     Geometry preprocessor
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype, &
                     coord, drmin, tarea, varea, af)

! Set initial condition
      call initialize(qc, cl, cd)

! For testing AD gradients
      !call test_resd(elem, edge, tedge, vedge, spts, coord, &
      !                 qc, qcd, qv, qvd, qx, qxd, qy, qyd, af, &
      !                 tarea, varea, dt, cl, cd, res)
      !call write_result(coord, elem, edge, qc, qv, cl, cd)
      !stop

!     call time_step2(edge, tedge, tarea, coord, qc, dt)
!     call jacobian(elem, esue, coord, tarea, dt, qc)
!     stop

      iter = 0
      fres = 1.0d0
      call system('rm -f FLO.RES')
      print*,'Beginning of iterations ...'
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
!        call time_step(drmin, qc, dt)
         call time_step2(edge, tedge, tarea, coord, qc, dt)
         call save_old(qc, qcold)

         do irk=1,nirk

!           Compute finite volume residual
            call fvresidual(elem, edge, tedge, vedge, spts, &
                            coord, qc, qv, qx, qy, af, tarea, varea, &
                            cl, cd, res)

!           Update the solution
            if(timemode .eq. 1)then
               do i=1,nt
                  call prim2con(qcold(1,i), c1)
                  call prim2con(qc(1,i),    c2)
                  do j=1,nvar
                     c3(j) = airk(irk)*c1(j) + &
                             birk(irk)*(c2(j)-dt(i)*res(j,i)/tarea(i))
                  enddo
                  call con2prim(c3, qc(1,i))
               enddo
            elseif(timemode .eq. 2)then
               call lusgs(elem, esue, edge, tedge, coord, qcold, qc, &
                          res, dt, tarea)
            elseif(timemode .eq. 3)then
               call gmres(elem, esue, edge, tedge, vedge, spts, &
                          coord, qc, qcd, qv, qvd, qx, qxd, qy, qyd, &
                          af, tarea, varea, dt, &
                          cl, cd, res)
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
            !cfl = 1.0e20                ! for subsonic
            !if(iter.gt.20) cfl = 1.0e20 ! for transonic
            !cfl = max(1.0d0, 10.0d0/fres)
            !cfl = -2.0d0 + 3.0d0*iter
            cfl = 1.1d0*cfl*fres_old/fres
         endif

      enddo

!     Save final solution
      call average(spts, elem, edge, coord, tarea, af, qc, qv)
      call write_result(coord, elem, edge, qc, qv, cl, cd)
      call write_sol(iter, fres, cl, cd)

!     Pressure matching objective function
      open(25,file='p.dat')
      write(25,*)(qv(4,i),i=1,np)
      close(25)
!     call pres_match(edge, coord, qv)

      totaltime = etime(elapsed)
      totaltime = totaltime/60.0
      elapsed(1)= elapsed(1)/60.0
      elapsed(2)= elapsed(1)/60.0
      print *, 'Time: total=', totaltime, ' user=', elapsed(1), &
               ' system=', elapsed(2)

      stop
      end
