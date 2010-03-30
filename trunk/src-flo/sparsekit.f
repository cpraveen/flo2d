c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(*),iao(*),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= dimension of A.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
c-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(*),iao(*),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
      end

c######################################################################
      subroutine ilu0(n, a, ja, ia, alu, jlu, ju, iw, ierr)
c######################################################################
      implicit real*8 (a-h,o-z)
      real*8 a(*), alu(*)
      integer ja(*), ia(*), ju(*), jlu(*), iw(*)
c------------------ right preconditioner ------------------------------*
c                    ***   ilu(0) preconditioner.   ***                *
c----------------------------------------------------------------------*
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of the L+U matrix is
c the same as that the A matrix, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c ILU0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c-----------------------------------------------------------------------
c
c on entry:
c---------
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c-----------
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju	  = pointer to the diagonal elements in alu, jlu.
c
c ierr	  = integer indicating error code on return
c	     ierr = 0 --> normal return
c	     ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c-------------
c iw	    = integer work array of length n.
c------------
c IMPORTANT
c-----------
c it is assumed that the the elements in the input matrix are stored
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L sorted by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling ilu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c
c-----------------------------------------------------------------------
        ju0 = n+2
        jlu(1) = ju0
c
c initialize work vector to zero's
c
      do 31 i=1, n
           iw(i) = 0
 31     continue
c
c main loop
c
      do 500 ii = 1, n
           js = ju0
c
c generating row number ii of L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              endif
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c
c     exit if diagonal element is reached.
c
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c
c     perform  linear combination
c
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
 140          continue
 150       continue
c
c     invert  and store diagonal element.
c
           if (alu(ii) .eq. 0.0d0) goto 600
           alu(ii) = 1.0d0/alu(ii)
c
c     reset pointer iw to zero
c
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c
c     zero pivot :
c
 600       ierr = ii
c
           return
           end

c######################################################################
      subroutine lusol(n, y, x, alu, jlu, ju)
c######################################################################
      real*8 x(*), y(*), alu(*)
      integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c This routine solves the system (LU) x = y, 
c given an LU decomposition of a matrix stored in (alu, jlu, ju) 
c modified sparse row format 
c
c-----------------------------------------------------------------------
c on entry:
c n   = dimension of system 
c y   = the right-hand-side vector
c alu, jlu, ju 
c     = the LU matrix as provided from the ILU routines. 
c
c on return
c x   = solution of LU x = y.     
c-----------------------------------------------------------------------
c 
c Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
c       will solve the system with rhs x and overwrite the result on x . 
c
c-----------------------------------------------------------------------
c local variables
c
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
      do 90 i = n, 1, -1
      do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91   continue
           x(i) = alu(i)*x(i)
 90   continue
c
      return
      end
c######################################################################
      subroutine lutsol(n, y, x, alu, jlu, ju) 
c######################################################################
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c This routine solves the system  Transp(LU) x = y,
c given an LU decomposition of a matrix stored in (alu, jlu, ju) 
c modified sparse row format. Transp(M) is the transpose of M. 
c----------------------------------------------------------------------- 
c on entry:
c n   = dimension of system 
c y   = the right-hand-side vector
c alu, jlu, ju 
c     = the LU matrix as provided from the ILU routines. 
c
c on return
c x   = solution of transp(LU) x = y.   
c-----------------------------------------------------------------------
c
c Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju) 
c       will solve the system with rhs x and overwrite the result on x . 
c 
c-----------------------------------------------------------------------
c local variables
c
        integer i,k
c
        do 10 i = 1, n
           x(i) = y(i)
 10     continue
c
c forward solve (with U^T)
c
        do 20 i = 1, n
           x(i) = x(i) * alu(i)
           do 30 k=ju(i),jlu(i+1)-1
              x(jlu(k)) = x(jlu(k)) - alu(k)* x(i)
 30        continue
 20     continue
c     
c     backward solve (with L^T)
c     
      do 40 i = n, 1, -1 
         do 50 k=jlu(i),ju(i)-1
              x(jlu(k)) = x(jlu(k)) - alu(k)*x(i)
 50      continue
 40   continue
c
      return
      end
c----------------end of lutsol -----------------------------------------
c-----------------------------------------------------------------------
      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
      implicit none 
      integer n
      real*8 a(*),alu(*),w(n)
      integer ja(*),ia(n+1),jlu(*),ju(n),levs(*),jw(3*n),lfil,iwk,ierr
c----------------------------------------------------------------------* 
c     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
c----------------------------------------------------------------------*
c
c on entry:
c========== 
c n       = integer. The row dimension of the matrix A. The matrix 
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.              
c
c lfil    = integer. The fill-in parameter. Each element whose
c           leve-of-fill exceeds lfil during the ILU process is dropped.
c           lfil must be .ge. 0 
c
c tol     = real*8. Sets the threshold for dropping small terms in the
c           factorization. See below for details on dropping strategy.
c  
c iwk     = integer. The minimum length of arrays alu, jlu, and levs.
c
c On return:
c===========
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c levs    = integer (work) array of size iwk -- which contains the 
c           levels of each element in alu, jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered in A or U.
c
c work arrays:
c=============
c jw      = integer work array of length 3*n.
c w       = real work array of length n 
c
c Notes/known bugs: This is not implemented efficiently storage-wise.
c       For example: Only the part of the array levs(*) associated with
c       the U-matrix is needed in the routine.. So some storage can 
c       be saved if needed. The levels of fills in the LU matrix are
c       output for information only -- they are not needed by LU-solve. 
c        
c----------------------------------------------------------------------
c w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
c jw(n+1:2n)  stores the nonzero indicator. 
c 
c Notes:
c ------
c All the diagonal elements of the input matrix must be  nonzero.
c
c----------------------------------------------------------------------* 
c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2,
     *     jlev, min 
      real*8 t, s, fact 
      if (lfil .lt. 0) goto 998
c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      n2 = n+n 
      ju0 = n+2
      jlu(1) = ju0
c
c     initialize nonzero indicator array + levs array -- 
c
      do 1 j=1,2*n 
         jw(j)  = 0
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
c     
c     unpack L-part and U-part of row of A in arrays w 
c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (t .eq. 0.0) goto 170 
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n2+lenl) = 0 
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
               jw(n2+ii) = 0 
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n2+jpos) = 0 
               jw(n+k) = jpos
            endif
 170     continue
c
         jj = 0
c
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jw(n+  (pointers/ nonzero indicator).
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in jw(n2+  (levels) 
            j = jw(n2+jj) 
            jw(n2+jj)  = jw(n2+k) 
            jw(n2+k) = j
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c
c     zero out element in row by resetting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c     
c     get the multiplier for row to be eliminated (jrow) + its level
c     
         fact = w(jj)*alu(jrow)
         jlev = jw(n2+jj) 
         if (jlev .gt. lfil) goto 150
c
c     combine current row and row jrow
c
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
                  jw(n2+i) = jlev+levs(k)+1 
               else
c
c     this is not a fill-in element 
c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            else
c     
c     dealing with lower part.
c     
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
                  jw(n2+lenl) = jlev+levs(k)+1 
               else
c
c     this is not a fill-in element 
c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            endif
 203     continue
         w(jj) = fact
         jw(jj)  = jrow
         goto 150 
 160     continue 
c     
c     reset double-pointer to zero (U-part) 
c     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c
c     update l-matrix
c         
         do 204 k=1, lenl 
            if (ju0 .gt. iwk) goto 996
            if (jw(n2+k) .le. lfil) then
               alu(ju0) =  w(k)
               jlu(ju0) =  jw(k)
               ju0 = ju0+1
            endif
 204     continue
c     
c     save pointer to beginning of row ii of U
c     
         ju(ii) = ju0
c
c     update u-matrix
c
         do 302 k=ii+1,ii+lenu-1 
            if (jw(n2+k) .le. lfil) then
               jlu(ju0) = jw(k)
               alu(ju0) = w(k)
               levs(ju0) = jw(n2+k) 
               ju0 = ju0+1
            endif
 302     continue

         if (w(ii) .eq. 0.0) goto 999 
c     
         alu(ii) = 1.0d0/ w(ii) 
c     
c     update pointer to beginning of next row of U.
c     
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
c
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      return
c     
c     zero row encountered in A or U. 
c     
 999  ierr = -5
      return
c----------------end-of-iluk--------------------------------------------
      end
