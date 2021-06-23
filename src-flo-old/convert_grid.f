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
     +                 qv(nvar,npmax), tarea(ntmax),
     +                 drmin(ntmax), res(nvar,ntmax), c1(nvar),
     +                 c2(nvar), c3(nvar)
      double precision cl, cd
      double precision qx(3,npmax), qy(3,npmax), qcd(nvar,ntmax)
      real etime, elapsed(2), totaltime

      integer          i, j, irk

c     call math
c     call read_input
      gridfile = 'grid.in'
      timemode = 1
      call geometric(elem, edge, tedge, esue, vedge, spts, ptype,
     +               bdedge, esubp, coord, drmin, tarea, af)

c     Write out new grid
      open(20, file="grid.fm")
      write(20,*) np, nt, nbe
      do i=1,np
         write(20,'(i8,2e24.15)') i, coord(:,i)
      enddo
      do i=1,nt
         write(20,*) i, elem(:,i)
      enddo
      j = 0
      do i=nsw1,nsw2
         j = j+1
         write(20,*) j, edge(1,i), edge(2,i), solid
      enddo
      do i=nff1,nff2
         j = j+1
         write(20,*) j, edge(1,i), edge(2,i), farfield
      enddo
      close(20)
      print*,'Wrote converted grid into file grid.fm'

      stop
      end
