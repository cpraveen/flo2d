c     Variables for gmres
c     integer          ja(4*nvar*nvar*ntmax), ia(nvar*ntmax+1)
c     integer          jlu(4*nvar*nvar*ntmax), ju(nvar*ntmax+1)
c     integer          levs(nvar*ntmax)
c     integer          jw(3*nvar*ntmax)
c     double precision jmat(4*nvar*nvar*ntmax), alu(4*nvar*nvar*ntmax)
c     double precision w(nvar*ntmax)

      integer, parameter :: lfil=2

      integer,pointer,dimension(:)  :: ja,ia,jlu,ju,levs,jw
      real(dp),pointer,dimension(:) :: jmat,alu,w

      common/variables/ja,ia,jlu,ju,levs,jw,jmat,alu,w
