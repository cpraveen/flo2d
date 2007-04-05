C------------------------------------------------------------------------------
C Apply block diagonal preconditioner
C b = Dinv * a
C------------------------------------------------------------------------------
      subroutine precond(a, b)
      implicit none
      include 'param.h'
      include 'data.h'
      double precision a(*), b(*)

      integer          ns, i, j, k, jj, kk, c

      ns = nvar*nt
      c  = 1
      do i=1,ns,nvar
         do j=1,nvar
            jj = i + j - 1
            b(jj) = 0.0d0
            do k=1,nvar
               kk    = i + k - 1
               b(jj) = b(jj) + dmat(j,k,c)*a(kk)
            enddo
         enddo
         c = c + 1
      enddo

      return
      end
