!
    Subroutine OneTo2D_DoubleDataRedist5( rowcol, na, A, lda, nca, B, ldb, ncb, mblk, nblk, DescB, &
                                  Ms, nprocs, mpi_comm_rows, mpi_comm_cols, mpi_comm_global, Rtimes )
!
      use mpi
      use iso_c_binding
!  
      implicit none
! ..
! .. Scalar Parameter ..
      character(len=1), intent(in)     :: rowcol
      integer, intent(in)  :: na, lda, nca, ldb, ncb, mblk, nblk, nprocs, mpi_comm_rows, &
                              mpi_comm_cols, mpi_comm_global
! ..      
! .. Arrary Parametes ..
      integer, intent(in)  ::   Ms(*)
      integer, intent(in)  ::   Descb(*)
!#ifdef USE_ASSUMED_SIZE
      double precision, intent(inout)   :: A(lda,*)
      double precision, intent(out)     :: B(ldb,*)
      double precision, intent(out)     :: Rtimes(*)
!#else
!      complex(kind=(kind(1.0d0))), intent(inout)   :: A(lda,nca)
!      complex(kind=(kind(1.0d0))), intent(out)     :: B(ldb,ncb)
!#endif
!      
! .. Purpose ..
! ==============
!      
! This routine redistribute matrix A to matrix B in BCDD form. On input, A is a na-by-nca matrix, 
! which is stored in 1D block column form. Different processes may have different columns. 
! The numbers of processes are stored in array Ms(1:nprocs), where nprocs is the total number
! of processes. This routine consists of two phases. It first communicates among process rows in 
! the same column, and then it communicates among the process columns. 
!
! This version is designed for the case that some processes may not have matrix A, i.e., nca==0. 
! Some entries of vector Ms are also zero. 
!
! This routine allows the row block size mblk and column block size nblk of matrix B are different, i.e.,
! mblk != nblk.    
!      
!  .. Parameters ..
!  =================
!  rowcol   (input)     It equals 'R' or 'C', and shows the ordering of process grid in 
!                       row-major or column-major.
!     
!  mblk     (input)     The row block size of B for BCDD.   
!
!  nblk     (input)     The column block size of B for BCDD.
!      
! =========================
! Written by Shengguo Li, Sept. 17, 2020
!
! 1) This version works for process with no columns.
! ===================================================

      integer   :: nprow,npcol,rank,myrow,mycol,lrankr,lrankc,mpierr,context,   &
                   lwork, nc_local, itmp, istart,iend, lnseg, ptype,lenGs,lbnd, &
                   lbnd2, na_rows, it,lstart,partner,msend,mreceive,    &
                   icomrow,mseg_recv,ipcolumn,msendT,mreceiveT,nc_tmp,lnwork,   &
                   comm_cycles, nc_local1
      integer   :: i,j,k,j1,row_major,msi,n_seg,i_s,irow_type,irow_start,isend_start
	  double precision      :: times(4), time0, time1 
							! times(1): Index computation; 
							! times(2): message packing and unpacking; 
							! times(3): communication schedule;
							! times(4): data transfer
	  double precision, parameter :: zero = 0.0D+0

      integer, allocatable  :: Na_locals(:), PT_locals(:), sendcounts(:),sdispls(:),        &
                               recvcounts(:),rdispls(:),Gindex(:),Gs_index(:),Gp_index(:),  &
                               Gl_index(:),Gs_index2(:),V_send(:),V_receive(:),Gp_index2(:),&
                               Gp_recv_index(:)
      double precision, allocatable :: Asend(:), Arev(:,:), Asend2(:,:), Arev2(:,:) 

      logical, external     :: lsame
      integer, external     :: numroc, indxl2g, indxg2l
      
! ************** Get the process topology *********
      call MPI_Comm_rank( mpi_comm_global, rank, mpierr )

      context = DESCB(2)
      call BLACS_GRIDINFO( context,nprow,npcol,myrow,mycol )

!#ifdef TIMING
	time0 = MPI_Wtime()
	DO i = 1, 4
		times(i) = zero
	END DO
!#endif
!*************************************************************
!** Each process reorder the local vectors to row BCDD form. *
!** Data for other processes are stored continuously.        *
!*************************************************************

      lwork = max( na*nca,1 )  ! To be safe, lwork=1 when nca=0.
      allocate( Asend(lwork),Na_locals(nprow),PT_locals(nprow)  )
      DO i = 1, nprow
         Na_locals( i ) = numroc( na, mblk, i-1, 0, nprow )
      END DO

      IF( lsame(rowcol, 'R') ) THEN
         row_major= 1
      ELSE
         row_major= 0
      END IF

      ! Compute the sum of numbers of columns of the processes in the same column.
      nc_local=0
      IF (row_major .eq. 1 ) THEN  ! Row-major
         DO i = 1, nprow
            itmp = (i-1)*npcol+mycol+1
            nc_local = nc_local+Ms(itmp)
         END DO
      ELSE ! Column-major
         DO i = 1, nprow
            itmp = mycol*nprow+i
            nc_local = nc_local+Ms(itmp)
         END DO
      END IF

      na_rows = Na_locals( myrow+1 )
      nc_local1 = max(nc_local, 1)
      allocate( Arev(na_rows,nc_local1) )

      ! Initialize the starting position of each process row
      PT_locals(1) = 1
      DO i=2, nprow
         PT_locals(i) = PT_locals(i-1) + Na_locals(i-1)*nca
      END DO
	  
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(1) = times(1) + time1 
	time0 = MPI_Wtime()
!#endif
      ! Scan each column of local matrix
      n_seg = na/mblk
      IF (n_seg*mblk .lt. na) n_seg = n_seg + 1
      Do j = 1, nca
         Do i_s = 0, n_seg-2  ! The first n_seg-1 segments
            irow_type = mod(i_s,nprow)+1
            isend_start = PT_locals(irow_type) 
            irow_start = i_s*mblk+1
            
            ! Copy the next mblk entries of A to Asend
            DO i=irow_start, irow_start+mblk-1
               Asend(isend_start+i-irow_start) = A(i,j)
            END DO
            PT_locals(irow_type) = PT_locals(irow_type)+mblk
         END DO ! i_s
         
         ! The last segment may have less than mblk rows
         i_s = n_seg-1
         irow_type = mod(i_s,nprow)+1
         isend_start = PT_locals(irow_type)
         irow_start = i_s*mblk+1

         ! Copy the next segment of A to Asend
         DO i=irow_start, na
            Asend(isend_start+i-irow_start) = A(i,j)
         END DO
         PT_locals(irow_type) = PT_locals(irow_type)+(na-irow_start+1)
         
      END DO ! col j
	  
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(2) = times(2) + time1 
	time0 = MPI_Wtime()
!#endif
	  
!**********************************************************
!* Each process column performs MPI_AlltoAllv and gather  *
!* all the required data from other processes.            *
!**********************************************************

      allocate( sendcounts(nprow),sdispls(nprow),recvcounts(nprow),rdispls(nprow) )

      DO i=1, nprow
         sendcounts(i) = Na_locals(i)*nca
      END DO
      sdispls(1) = 0
      DO i=2, nprow
         sdispls(i) = sdispls(i-1) + Na_locals(i-1)*nca
      END DO

      IF( row_major .eq. 1) THEN
         DO i =1, nprow
            ! get the number of columns in the i-th process row
            itmp = (i-1)*npcol+mycol+1
            nc_tmp = Ms(itmp)
            recvcounts(i) = na_rows*nc_tmp
         END DO

         rdispls(1) = 0
         DO i=1, nprow-1
            itmp = (i-1)*npcol+mycol+1
            nc_tmp = Ms(itmp)
            rdispls(i+1) = rdispls(i) + na_rows*nc_tmp
         END DO
      ELSE                ! Column-major
         DO i =1, nprow
            ! get the number of columns in the i-th process row
            itmp = mycol*nprow+i
            nc_tmp = Ms(itmp) 
            recvcounts(i) = na_rows*nc_tmp
         END DO
         
         rdispls(1) = 0
         DO i=1, nprow-1
            itmp = mycol*nprow+i
            nc_tmp = Ms(itmp)
            rdispls(i+1) = rdispls(i) + na_rows*nc_tmp
         END DO
      END IF

!#if TIMING
	time1 = MPI_Wtime() -time0
	times(1) = times(1) + time1; 
	time0 = MPI_Wtime()
!#endif
	  
      CALL MPI_Alltoallv( Asend,sendcounts,sdispls,MPI_DOUBLE, &
           Arev,recvcounts,rdispls,MPI_DOUBLE,mpi_comm_rows, mpierr )
		   
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1
	time0 = MPI_Wtime()
!#endif

      deallocate( sendcounts, sdispls )
      deallocate( recvcounts, rdispls )
      deallocate( Na_locals, PT_locals )
      deallocate( Asend )

!*************************************************************      
!* Stage 2: Communications along the same process row        *
!*          Among the process columns                        *
!*************************************************************

      ! Allocate the workspace
      allocate( Gindex(2*nprow) )  ! There are nprow initial segments.

      n_seg = 0
      IF( row_major .eq. 1 ) THEN
         DO i = 1, nprow
            itmp = (i-1)*npcol+mycol+1  ! process's global ID +1
            
            istart=1
            DO j=1, itmp-1
               istart= istart + Ms(j) 
            END DO
            IF( Ms(itmp) > 0 ) THEN
               n_seg = n_seg + 1
               Gindex(2*n_seg-1) = istart
               Gindex(2*n_seg) = istart+Ms(itmp)-1
            END IF
         END DO
      ELSE 
         ! Column-major
         DO i=1, nprow
            itmp = mycol*nprow+i

            istart=1
            DO j=1, itmp-1
               istart = istart+ Ms(j)
            END DO
            IF( Ms(itmp) > 0 ) THEN
               n_seg = n_seg + 1     ! Total number segments of global continuous vectors
               Gindex(2*n_seg-1) = istart
               Gindex(2*n_seg) = istart+Ms(itmp)-1
            END IF
         END DO
      END IF

      DO i = 1, n_seg-1
         ! Check whether segments are neighbor
         if( Gindex(2*i)+1 == Gindex(2*i+1) ) then
            Gindex(2*i) = Gindex(2*(i+1))
            Do j = i+1, n_seg-1 ! move the right segments forward
               j1 = j+1
               Gindex(2*j-1) = Gindex(2*j1-1)
               Gindex(2*j) = Gindex(2*j1) 
            END DO
            
            n_seg = n_seg-1
         end if
      END DO

      ! Allocate the workspace. 
      lnwork = 2*n_seg + (nc_local/nblk) +1    ! n_seg and nc_local may both be zero ***
!      lnwork = 10*lnwork
      lnwork = na
      Allocate( Gs_index(2*lnwork),Gp_index(lnwork),Gl_index(lnwork+1) )

      Gp_index(1:lnwork) = -1
      ! Deal with each segment defined in [Gindex(2*i-1), Gindex(2*i)]
      lenGs = 0
      Gl_index(1) = 1
      DO i=1, n_seg
         istart = ceiling( real( Gindex(2*i-1)) / real(nblk) )
         iend   = ceiling( real(Gindex(2*i)) / real(nblk) )
         lnseg = iend-istart ! There are lnseg boundaries inside current segment.

         if( lnseg .eq. 0) then  ! current segment belongs to the same process
            ptype = mod( istart-1,npcol )
            lenGs = lenGs + 1
            Gs_index(2*lenGs-1) = Gindex(2*i-1)
            Gs_index(2*lenGs) = Gindex(2*i)

            Gp_index(lenGs)  = ptype
            Gl_index(lenGs+1) = Gl_index(lenGs)+( Gindex(2*i)-Gindex(2*i-1) )+1   ! The starting position of next segment. 

         else
            ! Deal with the first sub-segment
            ptype = mod( istart-1,npcol )
            lenGs = lenGs+1
            lbnd2 = istart*nblk  ! right boundary
            Gs_index(2*lenGs-1) = Gindex(2*i-1)
            Gs_index(2*lenGs) = lbnd2

            Gp_index(lenGs) = ptype
            Gl_index(lenGs+1) = Gl_index(lenGs) + (lbnd2-Gindex(2*i-1))+1
            
            DO j =1, lnseg-1  ! The subsegment between j-th and (j+1)-th boundaries.
               ptype = mod( ptype+1,npcol )
               lenGs = lenGs+1
               lbnd  = (istart-1+j)*nblk+1  ! left  boundary
               lbnd2 = (istart+j)*nblk      ! right boundary
               Gs_index(2*lenGs-1) = lbnd
               Gs_index(2*lenGs)   = lbnd2

               Gp_index(lenGs)   = ptype
               Gl_index(lenGs+1) = Gl_index(lenGs)+nblk
            END DO

            ! Deal with the last sub-segment
            j = lnseg
            ptype = mod( ptype+1,npcol )
            lenGs = lenGs+1
            lbnd  = (istart-1+j)*nblk+1
            lbnd2 = Gindex(2*i)

            Gs_index(2*lenGs-1) = lbnd
            Gs_index(2*lenGs)   = lbnd2

            Gp_index(lenGs) = ptype
            Gl_index(lenGs+1) = Gl_index(lenGs)+ (lbnd2-lbnd)+1 
         end if

      END DO  ! i=1, n_seg

!      call MPI_Barrier( mpi_comm_global, mpierr )
      
      deallocate( Gindex )
	  
!#if TIMING
	time1 = MPI_Wtime() - time0
	times(1) = times(1) + time1
	time0 = MPI_Wtime()
!#endif

      Allocate(Asend2(na_rows,nc_local1),V_send(npcol),V_receive(npcol),Gp_index2(npcol) )
      Allocate( Gp_recv_index(npcol),Gs_index2(2*lnwork) )

      ! Initialize these vectors to be zero
      V_send(1:npcol) = 0
      Gs_index2(1:2*lnwork) = 0
      Gp_index2(1:npcol) = 0
      
      it = 1
      istart = 1
      DO ptype =0, npcol-1
         V_send(ptype+1) = 0     ! Index starts from ONE.
         Gp_index2(ptype+1) = 0  ! Record the number of segments of each ptype
         DO i=1,lenGs
            if(Gp_index(i) .eq. ptype )  then
               itmp = Gl_index(i+1)-Gl_index(i)
               call dlacpy('A',na_rows,itmp,Arev(1,Gl_index(i)),na_rows,Asend2(1,istart),na_rows)
               
               V_send(ptype+1)    = V_send(ptype+1) + itmp  ! Record the total number of columns to be sent of certain type
               Gs_index2(2*it-1)  = Gs_index(2*i-1)         ! Record the global index of certain segment
               Gs_index2(2*it)    = Gs_index(2*i)
               Gp_index2(ptype+1) = Gp_index2(ptype+1) + 1  ! Record the number of segments of certain process
               it = it+1
               istart = istart+itmp  ! the column pointer
            end if
         END DO
      END DO

      ! Construct V_receive via communication. Each process holds one V_send vector.
      call mpi_alltoall(V_send,1,MPI_INT,V_receive,1,MPI_INT,mpi_comm_cols, mpierr)
      ! Construct Gp_recv_index, which records the number of segments to be received.
      call mpi_alltoall( Gp_index2,1,MPI_INT,Gp_recv_index,1,MPI_INT,mpi_comm_cols,mpierr )

      DO i=1, lenGs
         if( Gp_index(i) .eq. mycol ) then
            istart = Gs_index(2*i-1)
            itmp = Gl_index(i+1) - Gl_index(i)  ! Column size

            lstart = indxg2l( istart, nblk, mycol, 0, npcol )
            call dlacpy( 'A',na_rows,itmp,Arev(1,Gl_index(i)), na_rows, B(1,lstart),ldb )
         end if
      END DO

      deallocate( Arev )

!#if TIMING
	time1 = MPI_Wtime() -time0
	times(2) = times(2) + time1 
!#endif

!************************************************************************      
!*  Perform the Caterpillar approach for all pairwise communication     *
!*
!************************************************************************

      comm_cycles = npcol

      ! Allocate workspace Arev2 for receiving matrix entries
      lwork = maxval( V_receive(1:npcol) )  ! This is not optimal since self may be larger. 
      Allocate( Arev2(na_rows,lwork) )

      ! The communication is performed npcol times based on the Caterpillar algorithm.
      DO it =1, comm_cycles
         
         ! Each process find its corresponding communication partner
         partner = mod(npcol-it-mycol,npcol)
         partner = mod(partner+npcol, npcol)

10       CONTINUE
         ! if partner is itself, exit
         if ( partner .eq. mycol ) then
            Goto 20
         end if
         msend = V_send( partner+1 )         ! number of columns to send to partner
         mreceive  = V_receive( partner+1 )  ! number of columns to receive from partner
         
         if( msend .eq. 0 ) then
            if(mreceive .eq. 0 ) then
               GOTO 20
            else
!#if TIMING			   
	time0 = MPI_Wtime()	
!#endif
               mreceiveT = mreceive*na_rows
               call MPI_Recv( Arev2, mreceiveT, MPI_DOUBLE, partner, 1,&
                    mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1 
	time0 = MPI_Wtime()
!#endif
               ! *** if 2*mseg_recv > lwork, it will make trouble. ***
               mseg_recv = Gp_recv_index(partner+1) 
               call MPI_Recv(Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                    mpi_comm_cols,MPI_STATUS_IGNORE, mpierr )

               ! Compute the local index of each segment.
               ipcolumn = 1  ! point to Arev2
               DO i =1, mseg_recv
                  istart = Gs_index(2*i-1)
                  itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                  lstart = indxg2l( istart,nblk,mycol,0,npcol )
                  call dlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                  ipcolumn = ipcolumn+itmp
               END DO
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(2) = times(2) + time1
!#endif
            end if
         else ! (msend .ne. 0)
            if (mreceive .eq. 0) then
!#if TIMING			   
	time0 = MPI_Wtime()	
!#endif
               ! Send Asend2 first.
               istart = 1
               DO i=0, partner-1
                  istart = istart+V_send(i+1)
               END DO
               msendT = msend*na_rows
               call mpi_send( Asend2(1,istart),msendT,MPI_DOUBLE,partner,1,&
                    mpi_comm_cols,mpierr )
               
               ! Send Gs_index2 second.
               istart = 1
               DO i=0, partner-1
                  istart = istart + Gp_index2( i+1 )
               END DO
               call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                    mpi_comm_cols,mpierr )
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1 
!#endif
       
            else ! (mreceive .ne. 0)

               if( mycol .lt. partner ) then
!#if TIMING			   
	time0 = MPI_Wtime()	
!#endif				  
                  ! Send Asend2 first.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart+V_send(i+1)
                  END DO
                  msendT = msend*na_rows
                  call mpi_send( Asend2(1,istart),msendT,MPI_DOUBLE,partner,1,&
                       mpi_comm_cols,mpierr )

                  ! Send Gs_index2 secondly.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart + Gp_index2( i+1 )
                  END DO
                  call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                       mpi_comm_cols,mpierr )

                  ! Receive the columns of matrix.
                  ! *** If nc_locals < mreceive, it will break. ***
                  mreceiveT = mreceive*na_rows
                  call MPI_Recv( Arev2, mreceiveT, MPI_DOUBLE, partner, 1,&
                       mpi_comm_cols, MPI_STATUS_IGNORE, mpierr )

                  ! Copy the vectors into local matrix B.
                  ! *** if 2*mseg_recv > lwork, it will make trouble. ***
                  mseg_recv = Gp_recv_index(partner+1) 
                  call MPI_Recv(Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                       mpi_comm_cols,MPI_STATUS_IGNORE, mpierr )
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1
	time0 = MPI_Wtime()	
!#endif
                  ! Compute the local index of each segment.
                  ipcolumn = 1  ! point to Arev2
                  DO i =1, mseg_recv
                     istart = Gs_index(2*i-1)
                     itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                     lstart = indxg2l( istart,nblk,mycol,0,npcol )
                     call dlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                     ipcolumn = ipcolumn+itmp
                  END DO
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(2) = times(2) + time1
!#endif
               else
                  ! The larger process first receive and then send
!#if TIMING			   
	time0 = MPI_Wtime()	
!#endif				                    
                  ! Receive the columns of matrix.
                  ! *** If nc_locals < mreceive, it will break. ***
                  mreceiveT = mreceive*na_rows
                  call MPI_Recv( Arev2, mreceiveT, MPI_DOUBLE, partner, 1,&
                       mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

                  ! Copy the vectors into local matrix B.
                  ! *** if mseg_recv > lnwork, it will make trouble. ***
                  mseg_recv = Gp_recv_index(partner+1) 
                  call MPI_Recv( Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                       mpi_comm_cols,MPI_STATUS_IGNORE, mpierr )
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1
	time0 = MPI_Wtime()	
!#endif
                  ! Compute the local index of each segment.
                  ipcolumn = 1  ! point to Arev2
                  DO i =1, mseg_recv
                     istart = Gs_index(2*i-1)
                     itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                     lstart = indxg2l( istart,nblk,mycol,0,npcol )
                     call dlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                     ipcolumn = ipcolumn+itmp
                  END DO
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(2) = times(2) + time1
	time0 = MPI_Wtime()	
!#endif
                  ! Send Asend2 first.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart+V_send(i+1)
                  END DO
                  msendT = msend*na_rows
                  call mpi_send( Asend2(1,istart),msendT,MPI_DOUBLE,partner,1,&
                       mpi_comm_cols,mpierr )
               
                  ! Send Gs_index2 second.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart + Gp_index2( i+1 )
                  END DO
                  call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                       mpi_comm_cols,mpierr )
!#if TIMING
	time1 = MPI_Wtime() -time0
	times(4) = times(4) + time1
!#endif                  
               end if ! (mycol .lt. partner)
               
            end if ! (mreceive .eq. 0)
            
         end if  ! (msend .eq. 0)

         !********* Current communication is finished ****************
         !**************************************************************************

20       CONTINUE

         ! Call mpi_barrier
         call mpi_barrier( mpi_comm_cols, mpierr  )
         
      END DO ! it=1, npcol-1

!#if TIMING	  
	  call mpi_reduce( times, Rtimes, 4, mpi_double_precision, mpi_max, 0, mpi_comm_global,mpierr )
!#endif 

      deallocate( V_send, V_receive  )
      deallocate( Gp_recv_index, Gs_index2, Gs_index  )
      deallocate( Gp_index2, Gp_index, Gl_index )
      deallocate( Arev2, Asend2 )

    end Subroutine OneTo2D_DoubleDataRedist5
