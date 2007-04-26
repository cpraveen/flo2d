*************************************************************************
**                 DRIVER FOR THE GMRes CODE
*************************************************************************
      subroutine gmres(elem, edge, tedge, vedge, spts, bdedge,
     +                 coord, qc, qv, qx, qy, af, carea, dt, cl, cd,
     +                 res, qcd)
      implicit none
      include 'param.h'

      integer          elem(3,ntmax), edge(2,nemax), tedge(2,nemax),
     +                 vedge(2,nemax), spts(nspmax), bdedge(2,nbpmax)
      double precision coord(2,npmax), qc(nvar,ntmax), af(3,npmax),
     +                 qv(nvar,npmax), carea(ntmax), res(nvar,ntmax),
     +                 qcd(nvar,ntmax), dt(ntmax), qx(3,npmax),
     +                 qy(3,npmax), cl, cd

      integer lda, ldstrt, lwork
      parameter (lda = nvar*ntmax, ldstrt = 20)
      parameter (lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1)

      integer i, j, n, m, icount
      integer revcom, colx, coly, colz, nbscal
      integer irc(5), icntl(8), info(3)

      integer matvec, precondLeft, precondRight, dotProd
      parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)

      integer nout

      double precision  work(lwork)
      double precision  cntl(5), rinfo(2)

      double precision ZERO, ONE
      parameter (ZERO = 0.0d0, ONE = 1.0d0)

      double precision resd(nvar,ntmax)

c     Number of unknowns = 4 * number of triangles
      n = nvar*nt

***************************************************************
** Set the right-hand side in positions (n+1) to 2n of the array work.
** The right-hand side is chosen such that the exact solution
** is the vector of all ones.
***************************************************************

      icount = 0
      do i=1,nt
         do j=1,nvar
            icount         = icount + 1
            work(n+icount) = -res(j,i)
         enddo
      enddo

*********************************
** Choose the restart parameter
*********************************
      m = ldstrt

*******************************************************
** Initialize the control parameters to default value
*******************************************************

      call init_dgmres(icntl,cntl)

*************************
*c Tune some parameters
*************************
c     Residual tolerance
      cntl(1) = gerrtol

* Save the convergence history on standard output
      icntl(1) = 0
      icntl(2) = 0
      icntl(3) = 0
* Maximum number of iterations
      icntl(7) = gmaxiter
*
* preconditioner location
      icntl(4) = prectype
* orthogonalization scheme
      icntl(5)=0
* initial guess
c     icntl(6) = 0
* residual calculation strategy at restart
c     icntl(8) = 1
** Initialise the solution to zero
      do j = 1,n
        work(j) = ZERO
      enddo

      call diag_matrix(edge, tedge, coord, qc, carea, dt)

*****************************************
** Reverse communication implementation
*****************************************

10    call drive_dgmres(n,n,m,lwork,work,irc,icntl,cntl,info,rinfo)
      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)

      if (revcom.eq.matvec) then
c        perform the matrix vector product
c        work(colz) <-- A * work(colx)
         icount = 0
         do i=1,nt
            do j=1,nvar
               qcd(j,i) = work(colx + icount)
               icount   = icount + 1
            enddo
         enddo
         call fvresidual_q(elem, edge, tedge, vedge, spts, bdedge,
     +                     coord, qc, qv, qx, qy, af, carea, cl, cd, 
     +                     qcd, resd)
         icount = 0
         do i=1,nt
            do j=1,nvar
               work(colz+icount) = carea(i)*qcd(j,i)/dt(i) + resd(j,i)
               icount            = icount + 1
            enddo
         enddo
         goto 10

      else if (revcom.eq.precondLeft) then
c        perform the left preconditioning
c        work(colz) <-- M^{-1} * work(colx)
         call precond(work(colx),work(colz))
         goto 10

      else if (revcom.eq.precondRight) then
c        perform the right preconditioning
         call precond(work(colx),work(colz))
         goto 10

      else if (revcom.eq.dotProd) then
c        perform the scalar product
c        work(colz) <-- work(colx) work(coly)
         call dgemv('C',n,nbscal,ONE,work(colx),n,
     &               work(coly),1,ZERO,work(colz),1)
         goto 10
      endif

*******************************
* dump the solution on a file
*******************************

      icount = 0
      do i=1,nt
         do j=1,nvar
            icount   = icount + 1
            qcd(j,i) = work(icount)
         enddo
      enddo

      nout = 11
      open(nout,FILE='sol_dTest',STATUS='unknown')
      if (icntl(5).eq.0) then
        write(nout,*) 'Orthogonalisation : MGS'
      elseif (icntl(5).eq.1) then
        write(nout,*) 'Orthogonalisation : IMGS'
      elseif (icntl(5).eq.2) then
        write(nout,*) 'Orthogonalisation : CGS'
      elseif (icntl(5).eq.3) then
        write(nout,*) 'Orthogonalisation : ICGS'
      endif
      write(nout,*) 'Restart : ', m
      write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
      write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
      write(nout,*) 'Optimal workspace = ', info(3)
      write(nout,*) ' ************************************************ '
      write(nout,*)
      close(nout)
      return
      end
