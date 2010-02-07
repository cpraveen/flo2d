C-----------------------------------------------------------------------------
C.....Definition of some constants
C-----------------------------------------------------------------------------
      subroutine math
      implicit none
      include 'param.h'

      GAMMA        = 1.4d0
      GAMMA1       = GAMMA-1.0d0
      GAS_CONST    = 1.0d0
      M_PI         = 4.0d0*datan2(1.0d0, 1.0d0)

      prandtl      = 0.72d0
      prandtl_turb = 0.9d0

      ALBADA11     = 2.0d0/3.0d0
      ALBADA12     = 1.0d0 - ALBADA11

      ALBADA21     = 4.0d0/3.0d0
      ALBADA22     = 1.0d0 - ALBADA21

      xvstatus     = no

c Constants for Spalart-Allmaras model
      Cb1          = 0.1355d0
      Cb2          = 0.622d0
      sigma_sa     = 2.0d0/3.0d0
      kolm         = 0.41d0
      Cw1          = Cb1/kolm**2 + (1.0d0 + Cb2)/sigma_sa
      Cw2          = 0.3d0
      Cw3          = 2.0d0
      Cv1          = 7.1d0
      Cv2          = 5.0d0

      Cv11         = Cv1**3
      Cw31         = 1.0d0 + Cw3**6
      Cw32         = Cw3**6
      kolm2        = kolm**2
      Cb2Sig1      = (1.0d0 + Cb2)/sigma_sa
      Cb2Sig2      = Cb2/sigma_sa

      return
      end

C-----------------------------------------------------------------------------
C.....Read parameters from an input file and set freestream values
C-----------------------------------------------------------------------------
      subroutine read_input
      implicit none
      include 'param.h'
      integer inp, iargc, n, inpstatus
      character sdummy*32
      character*24 flux_type

      n = iargc()
      if(n .eq. 0)then
         print*,'You must specify an input file.'
         stop
      endif

      call getarg(1,inpfile)

      inp = 11
      open(unit=inp, file=inpfile, status='old')
      print*,'Reading parameters from ',inpfile
      read(inp,*)sdummy, istart
      read(inp,*)sdummy, iflow
      read(inp,*)sdummy, mach_inf
      read(inp,*)sdummy, aoa_deg
      read(inp,*)sdummy, Rey
      read(inp,*)sdummy, cfl
      read(inp,*)sdummy, timemode
      read(inp,*)sdummy, gmaxiter, prectype, gerrtol
      read(inp,*)sdummy, iterlast
      read(inp,*)sdummy, maxiter
      read(inp,*)sdummy, minres
      read(inp,*)sdummy, saveinterval
      read(inp,*)sdummy, niso
      read(inp,*)sdummy, flux_type
      read(inp,*)sdummy, ILIMIT
      read(inp,*)sdummy, vortex, xref, yref
      read(inp,*)sdummy, display
      read(inp,*)sdummy, gridfile
      close(inp)

      inpstatus = yes

      if(istart .ne. scratch .and. istart .ne. restart)then
         print*,'Unknown start option',istart
         print*,'Possible values: 1=scratch or 2=restart'
         inpstatus = no
      endif

      if(iflow .ne. inviscid .and. iflow .ne. laminar .and.
     &   iflow .ne. turbulent)then
         print*,'Unknown flow type',iflow
         print*,'Possible values: 1=inviscid, 2=laminar, 3=turbulent'
         inpstatus = no
      endif

      if(flux_type == 'roe')  iflux = iroe
      if(flux_type == 'kfvs') iflux = ikfvs
      if(iflux .ne. iroe .and. iflux .ne. ikfvs)then
         print*,'Unknown flux ', flux_type
         print*,'Possible values: roe, kfvs'
         inpstatus = no
      endif


      if(ilimit .ne. no .and. ilimit .ne. yes)then
         print*,'Unknown limiter option',ilimit
         print*,'Possible values: 0=no, 1=yes'
         inpstatus = no
      endif

      if(vortex .ne. yes .and. vortex .ne. no)then
         print*,'Unknown vortex option',vortex
         print*,'Possible values: 0=no, 1=yes'
         inpstatus = no
      endif

      if(inpstatus .eq. no) stop

      open(10, file=gridfile, status="old")
      rewind(10)
      read(10,*) np, nt
      write(*, '( " Number of points    :", i8)') np
      write(*, '( " Number of triangles :", i8)') nt
      close(10)

      return
      end

C-----------------------------------------------------------------------------
C.....Variables stored are primitive - density, u, v, pressure
C.....Initialize primitive variables to free stream values
C-----------------------------------------------------------------------------
      subroutine initialize(qc, cl, cd)
      implicit none
      include 'param.h'
      double precision qc(nvar, ntmax), cl, cd
      
      integer          j

      q_inf   = 1.0d0
      r_inf   = 1.0d0
      p_inf   = 1.0d0/(GAMMA*mach_inf**2)
      T_inf   = p_inf/(r_inf*GAS_CONST)
      aoa     = aoa_deg*M_PI/180.0d0
      u_inf   = q_inf*dcos(aoa)
      v_inf   = q_inf*dsin(aoa)
      ent_inf = p_inf/r_inf**GAMMA
      a_inf   = dsqrt(GAMMA*p_inf/r_inf)
      H_inf   = a_inf**2/GAMMA1 + 0.5d0*q_inf

c     Required by Sutherland law
      T_infd  = 300.0d0
      SCONST  = 110.4d0*T_inf/T_infd

c     Primitive variables in free-stream
      prim_inf(1) = r_inf
      prim_inf(2) = u_inf
      prim_inf(3) = v_inf
      prim_inf(4) = p_inf

c     Print some useful info
      if(iflow.eq.inviscid) print*,'Euler computation'
      if(iflow.eq.laminar)  print*,'Laminar Navier-Stokes computation'
      if(iflow.eq.turbulent)print*,'Turbulent Navier-Stokes computation'
      print*,'Free-stream values:'
      write(*,'(5x, " Mach number =", f8.4)')mach_inf
      write(*,'(5x, " AOA         =", f8.4)')aoa_deg
      write(*,'(5x, " u velocity  =", f8.4)')u_inf
      write(*,'(5x, " v velocity  =", f8.4)')v_inf
      write(*,'(5x, " Pressure    =", f8.4)')p_inf

      if(vortex .eq. yes)then
            print*,'Using point vortex correction for far-field points'
            write(*,'(" Vortex center =", 2e12.4)') xref, yref
      endif

C     Runge-Kutta time stepping
      NIRK    = 3
      airk(1) = 0.0d0
      airk(2) = 3.0d0/4.0d0
      airk(3) = 1.0d0/3.0d0
      birk(1) = 1.0d0
      birk(2) = 1.0d0/4.0d0
      birk(3) = 2.0d0/3.0d0

C     For implicit scheme, set to single stage RK
      if(timemode .ne. 1) NIRK = 1

      cl = 0.0d0
      cd = 0.0d0

      if(istart .eq. scratch)then
         print*,'Initializing solution to free stream values'
         do j=1,nt
            qc(1,j) = prim_inf(1)
            qc(2,j) = prim_inf(2)
            qc(3,j) = prim_inf(3)
            qc(4,j) = prim_inf(4)
         enddo

      else
         print*,'Initializing solution to old values'
         call read_flow(qc)
      endif

      return
      end

C-----------------------------------------------------------------------------
C Save flow solution into a file. Conserved variables are saved.
C-----------------------------------------------------------------------------
      subroutine save_flow(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='FLO.DAT')
      do i=1,nt
         write(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

C-----------------------------------------------------------------------------
C Read flow solution from file
C-----------------------------------------------------------------------------
      subroutine read_flow(qc)
      implicit none
      include 'common.h'
      include 'size.h'
      double precision qc(nvar,ntmax)

      integer          i, j

      open(unit=50, file='FLO.DAT', status='OLD')
      do i=1,nt
         read(50,'(4e20.10)') (qc(j,i), j=1,nvar)
      enddo
      close(50)

      return
      end

C-----------------------------------------------------------------------------
C Time-step from cfl condition
C-----------------------------------------------------------------------------
      subroutine time_step(drmin, q, dt)
      implicit none
      include 'param.h'
      double precision drmin(ntmax), q(nvar,ntmax), dt(ntmax)

      integer          i
      double precision d, u, v, p, s, a, dtv, mutot, mul, mut,
     +                 sutherland
      external         sutherland

      dtglobal = 1.0d20
      if(iflow .ne. inviscid)then
         do i=1,nt
            d        = q(1,i)
            u        = q(2,i)
            v        = q(3,i)
            p        = q(4,i)
            s        = dsqrt(u*u + v*v)
            a        = dsqrt(gamma*p/d)
            mul      = sutherland(q(1,i))
            mutot    = mul
            dtv      = 2.0d0*GAMMA*mutot/(d*prandtl)
            dt(i)    = cfl*drmin(i)**2/(drmin(i)*(s + a) + dtv)
            dtglobal = dmin1(dtglobal, dt(i))
         enddo
      else
         do i=1,nt
            d        = q(1,i)
            u        = q(2,i)
            v        = q(3,i)
            p        = q(4,i)
            s        = dsqrt(u*u + v*v)
            a        = dsqrt(gamma*p/d)
            dt(i)    = cfl*drmin(i)/(s + a)
            dtglobal = dmin1(dtglobal, dt(i))
         enddo
      endif

      return
      end

C-----------------------------------------------------------------------------
C Time-step from cfl condition: refined version
C-----------------------------------------------------------------------------
      subroutine time_step2(edge, tedge, tarea, coord, q, dt)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), tedge(2,nemax)
      double precision tarea(ntmax), coord(2,npmax), q(nvar,ntmax), 
     +                 dt(ntmax)

      integer          i, j, t1, t2, e1, e2
      double precision sumspeed(nt), qa(nvar), d, u, v, p, a, nx, ny, 
     +                 nl, un, ll, lr

      do i=1,nt
         sumspeed(i) = 0.0d0
      enddo

C Interior edges
      do i=1,nin
         t1 = tedge(1,i)
         t2 = tedge(2,i)
         e1 =  edge(1,i)
         e2 =  edge(2,i)
         do j=1,nvar
            qa(j) = 0.5d0*( q(j,t1) + q(j,t2) )
         enddo
         d  = qa(1)
         u  = qa(2)
         v  = qa(3)
         p  = qa(4)

         nx =  ( coord(2,e2) - coord(2,e1) )
         ny = -( coord(1,e2) - coord(1,e1) )
         nl =  dsqrt(nx*nx + ny*ny)

         un = u*nx + v*ny
         a  = dsqrt(gamma*p/d)*nl
         ll = dabs(un) + a

         sumspeed(t1) = sumspeed(t1) + ll
         sumspeed(t2) = sumspeed(t2) + ll
      enddo

C Boundary edges
      do i=nin+1,ne
         t1 = tedge(1,i)
         e1 =  edge(1,i)
         e2 =  edge(2,i)

         d  = qa(1)
         u  = qa(2)
         v  = qa(3)
         p  = qa(4)

         nx =  ( coord(2,e2) - coord(2,e1) )
         ny = -( coord(1,e2) - coord(1,e1) )
         nl =  dsqrt(nx*nx + ny*ny)

         un = u*nx + v*ny
         a  = dsqrt(gamma*p/d)*nl
         ll = dabs(un) + a

         sumspeed(t1) = sumspeed(t1) + ll
      enddo

      dtglobal = 1.0d20
      do i=1,nt
         dt(i)    = cfl*tarea(i)/sumspeed(i)
         dtglobal = dmin1(dtglobal, dt(i))
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Save old solution
C-----------------------------------------------------------------------------
      subroutine save_old(q, qold)
      implicit none
      include 'param.h'
      double precision q(nvar,ntmax), qold(nvar,ntmax)

      integer          i, j

      do j=1,nvar
         do i=1,nt
            qold(j,i) = q(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------------
C convert primitive variable to conserved variable
C-----------------------------------------------------------------------------
      subroutine prim2con(prim, con)
      implicit none
      include 'common.h'
      double precision prim(nvar), con(nvar)

      double precision q2

      q2     = prim(2)**2 + prim(3)**2
      con(1) = prim(1)
      con(2) = prim(1)*prim(2)
      con(3) = prim(1)*prim(3)
      con(4) = prim(4)/GAMMA1 + 0.5d0*prim(1)*q2

      return
      end

C-----------------------------------------------------------------------------
C L2 and Linf norm of the finite volume residual
C-----------------------------------------------------------------------------
      subroutine residue(res)
      implicit none
      include 'param.h'
      double precision res(nvar,ntmax)

      integer          i, j
      double precision fr

      fres  = 0.0d0
      fresi = 0.0d0
      iresi = 0

      do i=1,nt
         fr = 0.0d0
         do j=1,nvar
            fr = fr + res(j,i)**2
         enddo
         fres = fres + fr
         fr   = dsqrt(fr)
         if(fr .gt. fresi)then
            fresi = fr
            iresi = i
         endif
      enddo

      fres = dsqrt(fres/nt)

      if(iter .eq. 1)then
         fres1 = fres
         print*,'Residue in first iteration =',fres1
      endif

      if(fres1 .ne. 0.0d0) fres = fres/fres1

      return
      end

C-----------------------------------------------------------------------------
C.....Error function, from Abromovitz and Stegun
C-----------------------------------------------------------------------------
      double precision function ERRF(X)
      double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

      ARG = X*X
      if(ARG .lt. 20.0d0)then
            E = dexp(-ARG)
      else
            E = 0.0d0
      endif
      VB = dabs(X)
      T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
      tmp1 = 1.061405429d0*T
      tmp2 = (tmp1 - 1.453152027d0)*T
      tmp3 = (tmp2 + 1.421413741d0)*T
      tmp4 = (tmp3 - 0.284496736d0)*T
      tmp5 = (tmp4 + 0.254829592d0)*T
      tmp6 = 1.0d0 - tmp5*E
      if(X .lt. 0.0d0)then
            ERRF = -tmp6
      else
            ERRF =  tmp6
      endif
      return
      end

C-----------------------------------------------------------------------------
C.....Prints data into a file for visualizing contours using gnuplot
C.....Taken from NSC2KE of Bijan Mohammadi
C-----------------------------------------------------------------------------
      SUBROUTINE isocont(ifile,F,COOR,NVAL,VAL)
      implicit none
      integer          ifile, nval
      double precision F(3),COOR(2,3),VAL(100)

      integer          IP1(3), ival, itr, k
      double precision epsi, ff1, ff2, ff3, ffma, d12, d23, val1, fk, 
     &                 fk1, fmi, fma, dif, eps, hh, x, y
C
      epsi   = 1.0d-5
      IP1(1) = 2
      IP1(2) = 3
      IP1(3) = 1
      FF1    = F(1)
      FF2    = F(2)
      FF3    = F(3)
      FFMA   = DMAX1(DABS(FF1),DABS(FF2))
      FFMA   = DMAX1(ffma,DABS(FF3))
      D12    = DABS(FF1-FF2)
      D23    = DABS(FF2-FF3)
      IF(D12+D23.LT.DMAX1(epsi,epsi*FFMA)) GOTO 1000
C  PAS DE RESTRICTION
C  ******************    
      DO 100 IVAL=1,NVAL
            VAL1 = VAL(IVAL)
            ITR  = 0
            DO 110 K=1,3
                  FK  = F(K)
                  FK1 = F(IP1(K))
                  FMI = DMIN1(FK,FK1)
                  FMA = DMAX1(FK,FK1)
                  DIF = FMA-FMI
                  IF(DIF.LT.epsi) GOTO 110
                  EPS = epsi*DIF
                  IF(VAL1.LT.FMI-EPS .OR. VAL1.GT.FMA+EPS) GOTO 110
                  HH  = DABS(FK-VAL1)/DIF
                  X   = COOR(1,K) + HH*(COOR(1,IP1(K))-COOR(1,K))
                  Y   = COOR(2,K) + HH*(COOR(2,IP1(K))-COOR(2,K))
                  IF(ITR.EQ.0) GOTO 115
                  write(ifile,*) x,y
                  write(ifile,*) 
                  GOTO 100
115               ITR = 1
                  write(ifile,*) x,y
110         CONTINUE
100   CONTINUE
1000  return      
      END
