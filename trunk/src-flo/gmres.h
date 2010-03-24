!     Variables for gmres
      !integer, parameter :: ntrimax=50000
      !integer          ja(4*nvar*nvar*ntrimax), ia(nvar*ntrimax+1)
      !integer          jlu(4*nvar*nvar*ntrimax), ju(nvar*ntrimax+1)
      !integer          levs(10*nvar*ntrimax)
      !integer          jw(3*nvar*ntrimax)
      !double precision jmat(4*nvar*nvar*ntrimax), alu(4*nvar*nvar*ntrimax)
      !double precision w(nvar*ntrimax)

      integer, parameter :: lfil=2

      integer,pointer,dimension(:)  :: ja, ia, jlu, ju, levs, jw
      real(dp),pointer,dimension(:) :: jmat, alu, w

      common/varint/ja,ia,jlu,ju,levs,jw
      common/varflt/jmat,alu,w
