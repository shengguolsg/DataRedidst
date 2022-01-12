!
    Subroutine TwoTo2D_DCplxRedistBLACS( rowcol, ma, na, A, lda, DescA, B, ldb, DescB )
!
!      use iso_c_binding
!  
      implicit none
! ..
! .. Scalar Parameter ..
      character(len=1), intent(in) :: rowcol
      integer, intent(in)          :: ma, na, lda, ldb
! ..      
! .. Arrary Parametes ..
      integer, intent(in)  ::   DescA(*)
      integer, intent(in)  ::   DescB(*)
      complex*16, intent(inout)   :: A(lda,*)
      complex*16, intent(out)     :: B(ldb,*)
!      
! .. Purpose ..
! ==============
!      
! This routine redistribute a BCDD matrix A to matrix B in BCDD form. These two
! data distribution have different block size, nblk1 and nblk2, respectively. 
! 
! **** This routine assume that the process grid does not changes. ****
! **** Matrix A and B are used in the same processes only with different blocksize. ****
!
! This routine is used to simulate our algorithm by using BLACS routine since their routines
! are robust and efficient for 1D process grid. It uses LCM or GCD, and communication pipelining.
! 
! The central difference from BLACS is that we further use row and column communicators,      
! and BLACS only uses one communicator. Therefor, our approach can be much faster than
! BLACS when using many processes. We use BLACS_GRIDMAP to construct row and column 
! communicators.
!
!  .. Parameters ..      
!  =================
!  
!
!        
!      
! =========================
! Written by Shengguo Li, Oct. 2nd, 2020
!
! =============================================

!     .. Parameters ..
      INTEGER BLOCK_CYCLIC_2D,DLEN_,DTYPE_,CTXT_,M_,N_,MB_,NB_,RSRC_, &
              CSRC_,LLD_,ITHVAL
      PARAMETER (BLOCK_CYCLIC_2D=1,DLEN_=9,DTYPE_=1,CTXT_=2,M_=3,N_=4, &
              MB_=5,NB_=6,RSRC_=7,CSRC_=8,LLD_=9,ITHVAL=10)
      integer, parameter :: ione = 1
!      
      integer   :: nprow,npcol,rank,myrow,mycol,info,context, na_rows,&
                   na_cols, na_rows2, na_cols2, row_major, itmp
      integer   :: i,j,k,IP, lwork, nprocs, nprowt, npcolt
      integer   :: nblkr1, nblkr2, nblkc1, nblkc2
!      
      integer   :: desc_row1(9), desc_col1(9), desc_row2(9), desc_col2(9)
      
      ! Used for construct column and row communincators. 
      integer, allocatable  :: umapr(:,:), umapc(:,:)
      integer, allocatable  :: contxtr(:), contxtc(:)
      integer, allocatable  :: myrowr(:), mycolc(:)
      
      ! Allocate two temple matrices, one for row communication, and
      ! another for column communication. 
      complex*16, allocatable :: AR(:,:)

      logical, external     :: lsame
      integer, external     :: numroc, indxg2l
      
! ************** Get the process topology *********

      context = DESCB(2)
      call BLACS_GRIDINFO( context,nprow,npcol,myrow,mycol )
      call BLACS_PINFO( rank, nprocs )

      nblkr1 = DESCA( MB_ )
      nblkc1 = DESCA( NB_ )
      nblkr2 = DESCB( MB_ )
      nblkc2 = DESCB( NB_ )

      IF(lsame(rowcol, 'R')) THEN
         row_major= 1
      ELSE
         row_major= 0
      END IF
      
!*************************************************************
!** Stage 1: Column Communication stage                      *
!** Each process column perform 1D data redistribution       *
!*************************************************************
      na_rows  = numroc( ma, nblkr1, myrow, 0, nprow )
      na_cols  = numroc( na, nblkc1, mycol, 0, npcol )
      na_rows2 = numroc( ma, nblkr2, myrow, 0, nprow )

      ! Allocate workspace for row communications (along column process)
      lwork = max( nprow, npcol )
      Allocate( AR(na_rows2,na_cols) )
      Allocate( umapr(nprow,npcol),myrowr(lwork),mycolc(lwork),contxtr(npcol) )

      IF( row_major .eq. 1 ) THEN
         DO IP = 1, npcol  ! each column is a process group
            DO J=1, nprow
               umapr(J,IP) = (J-1)*npcol+IP-1
            END DO
            CALL BLACS_GET( context, 0, contxtr(IP) )
            CALL BLACS_GRIDMAP( contxtr(IP),umapr(:,IP),nprow,nprow,ione )
         END DO
      ELSE
         ! Column Major 
         DO IP = 1, npcol ! each column is a process group
            DO J=1, nprow
               umapr(J,IP) = (IP-1)*nprow+J-1
            END DO
            CALL BLACS_GET( context, 0, contxtr(IP) )
            CALL BLACS_GRIDMAP( contxtr(IP),umapr(:,IP),nprow, nprow, ione )
         END DO
      END IF
      
      deallocate( umapr )

      myrowr(1:lwork) = -1
      mycolc(1:lwork) = -1
      DO IP =1, npcol
         CALL BLACS_Gridinfo( contxtr(IP),nprowt,npcolt,myrowr(IP),mycolc(IP) )
      END DO

      ! Set up a Scalapack row communicator
      IP = mycol+1
      IF( myrowr(IP).GT.-1 .AND. mycolc(IP).GT.-1 ) THEN
         
         CALL DESCINIT( desc_row1,ma,na_cols,nblkr1,na_cols,0, 0,contxtr(IP),na_rows,info )
         CALL DESCINIT( desc_row2,ma,na_cols,nblkr2,na_cols,0, 0,contxtr(IP),na_rows2,info )

         ! Execute 1D row data redistribution
         CALL PZGEMR2D( ma,na_cols,A,1,1,desc_row1,AR,1,1,desc_row2,contxtr(IP)  )
      END IF


!*************************************************************      
!* Stage 2: Communications along the same process row        *
!*          Among the process columns                        *
!*************************************************************

      ! Allocate workspace for column communicator (along row process)
      Allocate( umapc(npcol,nprow), contxtc(nprow) )

      IF( row_major .eq. 1 ) THEN
         DO IP = 1, nprow  ! each column is a process group
            DO J=1, npcol
               umapc(J,IP) = (IP-1)*npcol+J-1
            END DO
            CALL BLACS_GET( context, 0, contxtc(IP) )
            CALL BLACS_GRIDMAP( contxtc(IP),umapc(:,IP),ione,ione,npcol )
         END DO
      ELSE
         ! Column Major 
         DO IP = 1, nprow ! each column is a process group
            DO J=1, npcol
               umapc(J,IP) = (J-1)*nprow+IP-1
            END DO
            CALL BLACS_GET( context, 0, contxtc(IP) )
            CALL BLACS_GRIDMAP( contxtc(IP),umapc(:,IP),ione,ione,npcol )
         END DO
      END IF

      deallocate( umapc )

      myrowr(1:lwork) = -1
      mycolc(1:lwork) = -1
      DO IP =1, nprow
         CALL BLACS_Gridinfo( contxtc(IP),nprowt,npcolt,myrowr(IP),mycolc(IP) )
      END DO

      ! Set up a Scalapack column communicator
      IP = myrow+1
      IF( myrowr(IP).GT.-1 .AND. mycolc(IP).GT.-1 ) THEN
         
         CALL DESCINIT( desc_col1,na_rows2,na,na_rows2,nblkc1,0, 0,contxtc(IP),na_rows2,info )
         CALL DESCINIT( desc_col2,na_rows2,na,na_rows2,nblkc2,0, 0,contxtc(IP),na_rows2,info )

         ! Execute 1D column data redistribution
         CALL PZGEMR2D( na_rows2,na,AR,1,1,desc_col1,B,1,1,desc_col2,contxtc(IP)  )
      END IF
      
      deallocate( contxtr, contxtc )
      deallocate( myrowr,mycolc )
      deallocate( AR )

    end Subroutine TwoTo2D_DCplxRedistBLACS
