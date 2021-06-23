!-----------------------------------------------------------------------------
!.....Read parameters from an input file and set freestream values
!-----------------------------------------------------------------------------
      subroutine read_input
      implicit none
      include 'param.h'
      integer :: inp, iargc, n, inpstatus
      character sdummy*32

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
      read(inp,*)sdummy, flow_type
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
      read(inp,*)sdummy, scrinterval
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

      if(flow_type /= 'inviscid' .and. flow_type /= 'laminar' .and. &
         flow_type /= 'rans')then
         print*,'Unknown flow type',flow_type
         print*,'Possible values: inviscid, laminar, rans'
         inpstatus = no
      endif

      if(timemode /= 'rk3' .and. timemode /= 'lusgs' .and. &
         timemode /= 'gmres')then
         print*,'Unknown time integration scheme ', timemode
         print*,'Possible values: rk3, lusgs, gmres'
         inpstatus = no
      endif

      if(flux_type /= 'roe' .and. flux_type == 'kfvs')then
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
      read(10,*) np, nt, nbe
      close(10)

      ntmax = nt
      nemax = 3*np

      end
