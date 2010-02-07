C------------------------------------------------------------------------------
C Save the solution in various formats
C------------------------------------------------------------------------------
      subroutine write_result(coord, elem, edge, qc, qv, cl, cd)
      implicit none
      include 'param.h'
      integer  :: elem(3,*), edge(2,*)
      real(dp) :: coord(2,*), qc(nvar,*), qv(nvar,*), cl, cd

      integer  :: ifile, ii, is, it, i
      real(dp) :: ft1(3), ft2(3), coort(2,3), val1(100), 
     &            val2(100), delro1, delro2,
     &            deltat, xx0, yy0, xx1, yy1, q2, mach, cp, ent, 
     &            r, u, v, p
 
      call vigie(coord, elem, qv)
      call mayavi(coord, elem, qv)
      call save_flow(qc)

C     Write pressure coefficient
      open(unit=10, file='WALL.DAT')
      do i=nsw1,nsw2
         is = edge(1,i)
         r  = qv(1,is)
         u  = qv(2,is)
         v  = qv(3,is)
         q2 = u**2 + v**2
         p  = qv(4,is)
         cp = -(p - p_inf)/(0.5d0*r_inf*q_inf**2)
         ent= p/r**gamma/ent_inf - 1.0d0
         write(10,'(3e16.6)') coord(1,is), cp, ent
      enddo
      close(10)

c     iso-pressure/mach/visc for gnuplot
      open(22,file='FLO.P',status='unknown')
      open(23,file='FLO.M',status='unknown')
      rewind(22)
      rewind(23)
         
      delro1 = (pmax-pmin)/niso
      delro2 = (mmax-mmin)/niso
      do ii=1,niso+1
         val1(ii) = pmin       + (ii-1)*delro1
         val2(ii) = mmin       + (ii-1)*delro2
      enddo

      do it=1,nt
         do i=1,3
            is         = elem(i,it)
            coort(1,i) = coord(1,is)
            coort(2,i) = coord(2,is)
            r          = qv(1,is)
            u          = qv(2,is)
            v          = qv(3,is)
            q2         = u**2 + v**2
            p          = qv(4,is)
            ft1(i)     = p
            ft2(i)     = dsqrt(q2*r/(GAMMA*p))
         enddo
         call isocont(22,ft1,coort,niso,val1)
         call isocont(23,ft2,coort,niso,val2)
      enddo
      close(22)
      close(23)

cc    velocity vector for gnuplot
      ifile  = 26
      open(ifile,file='FLO.VECT',status='unknown')
      rewind(ifile)
      deltat = 1.0d-2
      do is=1,np
         xx0 = coord(1,is)
         yy0 = coord(2,is)
         u   = qv(2,is)
         v   = qv(3,is)
         xx1 = u*deltat
         yy1 = v*deltat
         write(ifile,*) xx0, yy0, xx1, yy1
      enddo
      close(ifile)

c     Run gnuplot and Start xv if not already started
      if(xvstatus .eq. no .and. iterlast .eq. 0)then
         call system('gnuplot flo.gnu')
         if(display.eq.1) call system('xv -poll flo.png &')
         if(display.eq.2) call system('display -update 10 flo.png &')
         xvstatus = yes
      else
         call system('gnuplot flo.gnu &')
      endif

      end
C------------------------------------------------------------------------------
C Save solution for vigie
C------------------------------------------------------------------------------
      subroutine vigie(coord, elem, prim)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), prim(nvar,*)

      integer  :: i, vig
      real(dp) :: q2, mach

      vig = 55
      open(unit=vig, file='FLO.VIG')

      write(vig,*)'points',np
      do i=1,np
         write(vig,*) coord(1,i), coord(2,i)
      enddo

      write(vig,*)'triangles',nt
      do i=1,nt
         write(vig,*) elem(1,i), elem(2,i), elem(3,i)
      enddo

      write(vig,*)'scalars  Mach'
      do i=1,np
         q2   = prim(2,i)**2 + prim(3,i)**2
         mach = dsqrt(q2*prim(1,i)/(GAMMA*prim(4,i)))
         write(vig,*) mach
      enddo

      write(vig,*)'scalars  Pressure'
      do i=1,np
         write(vig,*) prim(4,i)
      enddo

      write(vig,*)'scalars  Density'
      do i=1,np
         write(vig,*) prim(1,i)
      enddo

      write(vig,*)'scalars  Entropy'
      do i=1,np
         write(vig,*) dlog10(prim(4,i)/prim(1,i)**GAMMA/ent_inf )
      enddo

      write(vig,*)'vectors  vel  u  v  1e-02'
      do i=1,np
         write(vig,*) prim(2,i), prim(3,i)
      enddo
      write(vig,*)'end_block'

      close(vig)

      end
C------------------------------------------------------------------------------
C.....Result in VTK format for MayaVi/ParaView
C------------------------------------------------------------------------------
      subroutine mayavi(coord, elem, prim)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), prim(nvar,*)

      integer  :: i, maya
      real(dp) :: q2, mach

      maya = 55
      open(unit=maya, file='FLO.vtk')

      write(maya,10)
10    format('# vtk DataFile Version 3.0')

      if(iflow .eq. inviscid)then
      write(maya,111) mach_inf, aoa_deg
111   format('Mach =', f6.3, 2x, ' AOA = ', f6.3)
      else
      write(maya,112) mach_inf, aoa_deg, Rey
112   format('Mach =', f6.3, 2x, ' AOA = ', f6.3, ' Reynolds = ', e10.4)
      endif

      write(maya,12)
12    format('ASCII')
      write(maya,13)
13    format('DATASET UNSTRUCTURED_GRID')

      write(maya,14) np
14    format('POINTS ', i8, 2x, 'float')
      do i=1,np
         write(maya,'(3e18.8)') coord(1,i), coord(2,i), 0.0
      enddo

      write(maya,15) nt, 4*nt
15    format('CELLS ', 2i10)
      do i=1,nt
         write(maya,'(i4,3i10)') 3, elem(1,i)-1, elem(2,i)-1, 
     &                           elem(3,i)-1
      enddo

      write(maya,16) nt
16    format('CELL_TYPES ', i10)
      do i=1,nt
         write(maya,'(i5)') 5
      enddo

      write(maya,17) np
17    format('POINT_DATA ', i10)
      write(maya,18) 'Mach'
18    format('SCALARS ', a10, '   float 1')
      write(maya,19)
19    format('LOOKUP_TABLE default')
      do i=1,np
         q2   = prim(2,i)**2 + prim(3,i)**2
         mach = dsqrt(q2*prim(1,i)/(GAMMA*prim(4,i)))
         write(maya,20) mach
      enddo
20    format(e18.8)

      write(maya,18) 'Pressure'
      write(maya,19)
      do i=1,np
         write(maya,20) prim(4,i)
      enddo

      write(maya,18) 'Density'
      write(maya,19)
      do i=1,np
         write(maya,20) prim(1,i)
      enddo

      write(maya,18) 'Entropy'
      write(maya,19)
      do i=1,np
         write(maya,20) dlog10(prim(4,i)/prim(1,i)**GAMMA/ent_inf)
      enddo

      write(maya,25) 'Velocity'
      write(maya,19)
      do i=1,np
         write(maya,'(3e18.8)') prim(2,i), prim(3,i), 0.0
      enddo

25    format('VECTORS ', a10, '   float')

      close(maya)

      end
c------------------------------------------------------------------------------
c     Write some solution values: Used for optimization
c------------------------------------------------------------------------------
      subroutine write_sol(iter, residue, cl, cd)
      implicit none
      include 'common.h'
      integer  :: iter
      real(dp) :: residue, cl, cd

      integer  :: fid

      fid = 20
      open(fid, file="FLO.OUT")
      write(fid,*) iter, residue, cl, cd
      close(fid)

      end
