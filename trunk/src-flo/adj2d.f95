      program adj2d
      implicit none
      include 'param.h'
      include 'gmres.h'
      include 'adj.h'

      integer, parameter :: mode = 2
      integer, allocatable, dimension(:)   :: ptype
      integer, allocatable, dimension(:,:) :: elem, esue, edge, tedge, &
                                              vedge
      real(dp),allocatable, dimension(:)   :: dt, drmin, tarea, varea
      real(dp),allocatable, dimension(:,:) :: coord, qc, qcold, af, &
                                              qv, res, qx, qy
      real(dp),allocatable, dimension(:,:) :: qcd, qvd, qxd, qyd, resd
      real(dp),allocatable, dimension(:,:) :: con, conb, qcb, qvb, qxb, &
                                              qyb, resb, resbold, dJdu

      ! Arrays for gmres                                        
      integer,target,allocatable,dimension(:)  :: a_ja, a_ia, a_jlu, &
                                                  a_ju, a_levs, a_jw
      real(dp),target,allocatable,dimension(:) :: a_jmat, a_alu, a_w

      integer  ::  spts(nspmax)
      real(dp) :: c1(nvar), c2(nvar), c3(nvar)
      real(dp) :: cl, cd
      real(dp) :: cost, costb
      real     :: etime, elapsed(2), totaltime

      integer  :: i, j, irk

      costtyp = 'CL'

!     Set some constants
      call math

!     Read input parameters from file, also sets np, nt
      call read_input

      ! TEMPORARY STUFF
      timemode= 'rk3'
      nirk    = 3

!     Allocate memory for variables
      include 'allocate.h'
      include 'allocadj.h'

!     Geometry preprocessor
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype, &
                     coord, drmin, tarea, varea, af)

! Set initial condition
      istart = restart
      call initialize(qc, cl, cd)

      call clcd(edge, tedge, coord, qc, cl, cd)

!     call time_step(drmin, qc, dt)
      call time_step2(edge, tedge, tarea, coord, qc, dt)

      qcb   = 0.0d0
      qvb   = 0.0d0
      qxb   = 0.0d0
      qyb   = 0.0d0
      costb = 1.0d0
      call COSTFUN_BQ(elem, edge, tedge, vedge, spts, coord, qc, qcb, &
                      qv, qvb, qx, qxb, qy, qyb, af, tarea, varea, &
                      cl, cd, cost, costb)
      do i=1,nt
         call prim2con(qc(1,i), con(1,i))
         call con2prim_bq(con(1,i), dJdu(1,i), qc(1,i), qcb(1,i))
      enddo

      iter = 0
      fres = 1.0d0
      resb = 0.0d0
      qvb  = 0.0d0
      qxb  = 0.0d0
      qyb  = 0.0d0
      call system('rm -f ADJ.RES')
      print*,'Beginning of iterations ...'
      do while(iter .lt. MAXITER .and. fres .gt. MINRES)
         call save_old(resb, resbold)

         do irk=1,nirk

            ! Compute finite volume residual
            ! Use conb as temp variable. On output, conb=0
            qcb = 0.0d0
            conb = resb
            qvb  = 0.0d0
            qxb  = 0.0d0
            qyb  = 0.0d0
            call fvresidual_bq(elem, edge, tedge, vedge, spts, &
                            coord, qc, qcb, qv, qvb, qx, qxb, qy, qyb, &
                            af, tarea, varea, cl, cd, res, conb)

!           Update the solution
            if(timemode == 'rk3')then
               do i=1,nt
                  call con2prim_bq(con(1,i), conb(1,i), qc(1,i), qcb(1,i))
                  conb(:,i) = conb(:,i) + dJdu(:,i)
                  do j=1,nvar
                     resb(j,i) = airk(irk)*resbold(j,i) + &
                                 birk(irk)*(resb(j,i) - &
                                 dt(i)*conb(j,i)/tarea(i))
                  enddo
               enddo
            else
               print*,'adj2d: Time integration scheme not done'
               stop
            endif

         enddo

         iter = iter + 1
         call residue(conb)
         write(*,'(i6,4e16.6)') iter, fres, fresi, cl, cd
         open(unit=99, file='ADJ.RES', access='append')
         write(99,'(i6,4e16.6)') iter, fres, fresi, cl, cd
         close(99)

         if(mod(iter,saveinterval) == 0)then
            call average(spts, elem, edge, coord, tarea, af, resb, qv)
            call write_vtk(coord, elem, qv)
         endif

      enddo

      totaltime = etime(elapsed)
      totaltime = totaltime/60.0
      elapsed(1)= elapsed(1)/60.0
      elapsed(2)= elapsed(1)/60.0
      print *, 'Time: total=', totaltime, ' user=', elapsed(1), &
               ' system=', elapsed(2)

      stop
      end
