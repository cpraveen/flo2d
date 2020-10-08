!------------------------------------------------------------------------------
!.....Result in VTK format for fidVi/ParaView
!------------------------------------------------------------------------------
      subroutine write_vtk(coord, elem, var)
      implicit none
      include 'param.h'
      integer  :: elem(3,*)
      real(dp) :: coord(2,*), var(nvar,*)

      integer  :: i, fid

      fid = 55
      open(unit=fid, file='out.vtk')

      write(fid,10)
10    format('# vtk DataFile Version 3.0')

      write(fid,111) mach_inf, aoa_deg
111   format('Mach =', f6.3, 2x, ' AOA = ', f6.3)

      write(fid,12)
12    format('ASCII')
      write(fid,13)
13    format('DATASET UNSTRUCTURED_GRID')

      write(fid,14) np
14    format('POINTS ', i8, 2x, 'float')
      do i=1,np
         write(fid,'(3e18.8)') coord(1,i), coord(2,i), 0.0
      enddo

      write(fid,15) nt, 4*nt
15    format('CELLS ', 2i10)
      do i=1,nt
         write(fid,'(i4,3i10)') 3, elem(1,i)-1, elem(2,i)-1, &
                                 elem(3,i)-1
      enddo

      write(fid,16) nt
16    format('CELL_TYPES ', i10)
      do i=1,nt
         write(fid,'(i5)') 5
      enddo

      write(fid,17) np
17    format('POINT_DATA ', i10)
      write(fid,18) 'Density'
18    format('SCALARS ', a10, '   float 1')
      write(fid,19)
19    format('LOOKUP_TABLE default')
      do i=1,np
         write(fid,20) var(1,i)
      enddo
20    format(f18.8)

      write(fid,18) 'xmomentum'
      write(fid,19)
      do i=1,np
         write(fid,20) var(2,i)
      enddo

      write(fid,18) 'ymomentum'
      write(fid,19)
      do i=1,np
         write(fid,20) var(3,i)
      enddo

      write(fid,18) 'Energy'
      write(fid,19)
      do i=1,np
         write(fid,20) var(4,i)
      enddo

      close(fid)

      end
