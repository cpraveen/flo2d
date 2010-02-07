C-----------------------------------------------------------------------------
C.....Read parameters from an input file and set freestream values
C-----------------------------------------------------------------------------
      subroutine read_input
      implicit none
      include 'param.h'
      integer :: inp, iargc, n, inpstatus
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
      close(10)
      write(*, '( " Number of points    :", i8)') np
      write(*, '( " Number of triangles :", i8)') nt

      ntmax = nt
      nemax = 3*np

      end
