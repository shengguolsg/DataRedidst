!
    Subroutine OneTo1D_ComplexDataRedistBip( na, A, lda, nca, B, ldb, ncb, nblk, &
                                          Ms, nprocs, mpi_comm_global )
!
      use mpi
      use mod_cardmpp
!      use iso_c_binding
!  
      implicit none
! ..
! .. Scalar Parameter ..
      integer, intent(in)  :: na, lda, nca, ldb, ncb, nblk, nprocs, mpi_comm_global
! ..      
! .. Arrary Parametes ..
      integer, intent(in)  ::   Ms(*)
      complex*16, intent(inout)   :: A(lda,*)
      complex*16, intent(out)     :: B(ldb,*)
!      
! .. Purpose ..
! ==============
!      
! This routine redistribute matrix A to matrix B in BCDD form. On input, A is a na-by-nca matrix, 
! which is stored in 1D block column form. Different processes may have different columns. 
! The numbers of processes are stored in array Ms(1:nprocs), where nprocs is the total number
! of processes. This routine assumes all the processes are in one dimension.
! We assume that      
! 1) Firstly, matrix A is stored irregularly, but the indexes of processes are 
!    continuous from P0 to nprocs-1.
! 2) Matrix A is stored columnwise, and its rows do not change. 
! 3) Note that A only contain ONE segment of continuous eigenvectors. 
!
! Matrix B is in 1D BCDD form with block size nblk. 
! There are two cases that may happen
! 1) nblk is small enough such that every process has some local columns of B
! 2) nblk is large, and some processes may have none of matrix B. 
!      
!  .. Parameters ..
!  =================
!  Ms   (input)     Integer Arrays of dimension nprocs. 
!        
!      
! =========================
! Written by Shengguo Li, 12th. Oct. 2020
!
! This version uses Bipartite matching for scheduling. 
! 
! =============================================

      integer   :: nprow,npcol,rank,mycol,mpierr,context,nc_local,    &
                   lwork, itmp, istart,iend, lnseg, ptype,lenGs,lbnd, &
                   lbnd2, na_rows, it,lstart,partner,msend,mreceive,  &
                   mseg_recv,ipcolumn,msendT,mreceiveT,nc_tmp,lnwork, &
                   comm_cycles,ldm, nmatch
      integer   :: i,j,k,j1,row_major,msi,n_seg,i_s,irow_type,irow_start,isend_start
      integer   :: maxdeg, maxms, nedge, nedgef, nnode, jstart, jend

      integer   :: Gindex(2)

      integer, allocatable  :: inode(:,:), match(:)
!      
      integer, allocatable  :: Gs_index(:),Gp_index(:),  &
                               Gl_index(:),Gs_index2(:),V_send(:),V_receive(:),Gp_index2(:),&
                               Gp_recv_index(:)
      complex*16, allocatable :: Asend2(:,:), Arev2(:,:) 

      logical, external     :: lsame
      integer, external     :: numroc, indxl2g, indxg2l
      
! ************** Get the process topology *********
      call MPI_Comm_rank( mpi_comm_global,rank,mpierr )

      npcol = nprocs  ! The process grid is 1D
      mycol = rank   

!** Compute the global indexes of local eigenvectors of each process

      istart = 1
      nc_local = Ms(rank+1)
      DO j = 0, rank-1
         istart = istart + Ms(j+1) 
      END DO
      Gindex(1) = istart
      Gindex(2) = istart+nc_local-1

      ! Allocate the workspace. 
      lnwork = (nca/nblk) +2 +1
      lnwork = 2*lnwork
      Allocate( Gs_index(2*lnwork),Gp_index(lnwork),Gl_index(lnwork+1) )

      ! Deal with each segment defined in [Gindex(2*i-1), Gindex(2*i)]
      n_seg = 1  ! We assume there is at most one segment. 
      lenGs = 0
      Gl_index(1) = 1
      i = 1
      IF (Gindex(2) .ge. Gindex(1) ) THEN  ! current segment is not empty
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
            ! There are lnseg boundaries inside the current segment.
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
         
      END IF  ! not empty

! **** Asend2 is used for communications ****

      Allocate( V_send(npcol),V_receive(npcol),Gp_index2(npcol) )
      Allocate( Gp_recv_index(npcol),Gs_index2(2*lnwork) )
      IF( nc_local .gt. 0 ) Allocate( Asend2(na,nc_local)  )

      ! na_local = nca, na_rows = na
      na_rows = na
      IF( Gindex(2) .ge. Gindex(1) ) THEN
         it = 1
         istart = 1
         DO ptype =0, npcol-1
            V_send(ptype+1) = 0     ! Index starts from ONE.
            Gp_index2(ptype+1) = 0  ! Record the number of segments of each ptype
            DO i=1,lenGs
               if(Gp_index(i) .eq. ptype )  then
                  itmp = Gl_index(i+1)-Gl_index(i)
                  call zlacpy('A',na_rows,itmp,A(1,Gl_index(i)),na_rows,Asend2(1,istart),na_rows)
                  
                  V_send(ptype+1)    = V_send(ptype+1) + itmp  ! Record the total number of columns to be sent of certain type
                  Gs_index2(2*it-1)  = Gs_index(2*i-1)         ! Record the global index of certain segment
                  Gs_index2(2*it)    = Gs_index(2*i)
                  Gp_index2(ptype+1) = Gp_index2(ptype+1) + 1  ! Record the number of segments of certain process
                  it = it+1
                  istart = istart+itmp  ! the column pointer
               end if
            END DO
         END DO
      ELSE
         V_send(1:npcol) = 0
         Gs_index2 = 0
         Gp_index2(1:npcol) = 0
      END IF

      ! Construct V_receive via communication. Each process holds one V_send vector.
      call mpi_alltoall( V_send,1,MPI_INT,V_receive,1,MPI_INT,mpi_comm_global, mpierr )
      ! Construct Gp_recv_index, which records the number of segments to be received.
      call mpi_alltoall( Gp_index2,1,MPI_INT,Gp_recv_index,1,MPI_INT,mpi_comm_global,mpierr )

      ! Copy the local vectors to the final appropriate position
      ! The following approach uses A instead of Asend2
      DO i=1, lenGs
         if( Gp_index(i) .eq. mycol ) then
            istart = Gs_index(2*i-1)
            itmp = Gl_index(i+1) - Gl_index(i)  ! Column size

            lstart = indxg2l( istart, nblk, mycol, 0, npcol )
            call zlacpy( 'A',na_rows,itmp,A(1,Gl_index(i)), na_rows, B(1,lstart),ldb )
         end if
      END DO

!************************************************************************      
!*  Perform the Bipartite matching for all pairwise communication       *
!*  Require smaller than npcol sychronization.                          *
!************************************************************************

      maxms = maxval( Ms(1:npcol) )
      maxdeg= ceiling( real(maxms)/real(nblk) ) + 2

      nnode = npcol
      nedge = maxdeg*nnode

      ldm = nedge
      allocate( inode(nedge,2),match(nnode) )

      jstart = 1    ! the global index of the first column of current segment
      nedgef = 1
      DO i =1, npcol
         jend = jstart + Ms(i) -1 ! the global index of the last column of current segment

         istart = ceiling( real(jstart)/real(nblk) )
         iend   = ceiling( real(jend)/real(nblk) )
         lnseg  = iend-istart

         if(lnseg .eq. 0 ) then  ! current segment belongs to the same process
            ptype = mod( istart-1, npcol)
            if( ptype .ne. (i-1) ) then
               inode(nedgef,1) = i              ! the index of node starts from one. 
               inode(nedgef,2) = ptype+1        ! using bipartite should add npcol
               nedgef = nedgef + 1
            end if
         else
            ! There are lnseg boundaries inside the current segment.

            ! Deal with the first sub-segment
            ptype = mod(istart-1, npcol)
            if( ptype .ne. (i-1) ) then    !! To avoid current sub-segment is local. 
               inode(nedgef,1) = i
               inode(nedgef,2) = ptype+1 ! add npcol+1
               nedgef = nedgef + 1
            end if

            DO j=1, lnseg-1  ! The subsegment between j-th and (j+1)-th boundaries
               ptype = mod(ptype+1, npcol)
               if( ptype .ne. (i-1) ) then
                  inode(nedgef,1) = i
                  inode(nedgef,2) = ptype+1 ! add npcol+1
                  nedgef = nedgef + 1
               end if
            END DO

            ! Deal with the last sub-segment
            ptype = mod(ptype+1, npcol)
            if( ptype .ne. (i-1) ) then
               inode(nedgef,1) = i
               inode(nedgef,2) = ptype+1 ! add npcol+1
               nedgef = nedgef + 1
            end if
            
         end if ! (lnseg .ne. 0)

         jstart = jend+1
      END DO  ! (i=1, npcol)

      nedgef = nedgef -1  ! Total number of edges
      nedge = nedgef

      ! Allocate workspace Arev2 for receiving matrix entries
      lwork = maxval( V_receive(1:npcol) )  ! This is not optimal since self may be larger. 
      Allocate( Arev2(na_rows,lwork) )

      ! The communication is performed comm_cycles times based on the caterpillar algorithm.
      it = 0
      DO while( nedgef.GT. 0 .AND. it < 20 )

         it = it+ 1
         call max_match( nnode,nedge,inode,ldm,match,nmatch )

         partner = match(rank+1)-1

10       CONTINUE
         ! if partner is invalid, exit
         if (partner .lt. 0 .OR. partner.eq.mycol) then
            Goto 20
         end if
         msend = V_send(partner+1)         ! number of columns to send to partner
         mreceive  = V_receive(partner+1)  ! number of columns to receive from partner
         
         if( msend .eq. 0 ) then
            if(mreceive .eq. 0 ) then
               ! These two processes do not need communication
               ! exit
               GOTO 20
            else
               mreceiveT = mreceive*na_rows
               call MPI_Recv( Arev2, mreceiveT, MPI_C_DOUBLE_COMPLEX, partner, 1,&
                    mpi_comm_global, MPI_STATUS_IGNORE, mpierr)

               ! *** if 2*mseg_recv > lwork, it will make trouble. ***
               mseg_recv = Gp_recv_index(partner+1) 
               call MPI_Recv(Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                    mpi_comm_global,MPI_STATUS_IGNORE, mpierr )

               ! Compute the local index of each segment.
               ipcolumn = 1  ! point to Arev2
               DO i =1, mseg_recv
                  istart = Gs_index(2*i-1)
                  itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                  lstart = indxg2l( istart,nblk,mycol,0,npcol )
                  call zlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                  ipcolumn = ipcolumn+itmp
               END DO
            end if
         else ! (msend .ne. 0)
            if (mreceive .eq. 0) then

               istart = 1
               DO i=0, partner-1
                  istart = istart+V_send(i+1)
               END DO
               msendT = msend*na_rows
               call mpi_send( Asend2(1,istart),msendT,MPI_C_DOUBLE_COMPLEX,partner,1,&
                    mpi_comm_global,mpierr )
               
               ! Send Gs_index2 second.
               istart = 1
               DO i=0, partner-1
                  istart = istart + Gp_index2( i+1 )
               END DO
               call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                    mpi_comm_global,mpierr )
       
            else ! (mreceive .ne. 0)

               if( mycol .lt. partner ) then
                  ! first send and then receive
                  ! Data is stored in Asend2, and the process just as above.

                  ! Send Asend2 first.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart+V_send(i+1)
                  END DO
                  msendT = msend*na_rows
                  call mpi_send( Asend2(1,istart),msendT,MPI_C_DOUBLE_COMPLEX,partner,1,&
                       mpi_comm_global,mpierr )

                  ! Send Gs_index2 secondly.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart + Gp_index2( i+1 )
                  END DO
                  call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                       mpi_comm_global,mpierr )

                  ! Receive the columns of matrix.
                  ! *** If nc_locals < mreceive, it will break. ***
                  mreceiveT = mreceive*na_rows
                  call MPI_Recv( Arev2, mreceiveT, MPI_C_DOUBLE_COMPLEX, partner, 1,&
                       mpi_comm_global, MPI_STATUS_IGNORE, mpierr)

                  ! Copy the vectors into local matrix B.
                  ! *** if 2*mseg_recv > lwork, it will make trouble. ***
                  mseg_recv = Gp_recv_index(partner+1) 
                  call MPI_Recv(Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                       mpi_comm_global,MPI_STATUS_IGNORE, mpierr )

                  ! Compute the local index of each segment.
                  ipcolumn = 1  ! point to Arev2
                  DO i =1, mseg_recv
                     istart = Gs_index(2*i-1)
                     itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                     lstart = indxg2l( istart,nblk,mycol,0,npcol )
                     call zlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                     ipcolumn = ipcolumn+itmp
                  END DO

               else
                  ! first receive and then send
                  
                  ! Receive the columns of matrix.
                  ! *** If nc_locals < mreceive, it will break. ***
                  mreceiveT = mreceive*na_rows
                  call MPI_Recv( Arev2, mreceiveT, MPI_C_DOUBLE_COMPLEX, partner, 1,&
                       mpi_comm_global, MPI_STATUS_IGNORE, mpierr)

                  ! Copy the vectors into local matrix B.
                  ! *** if mseg_recv > lnwork, it will make trouble. ***
                  mseg_recv = Gp_recv_index(partner+1) 
                  call MPI_Recv(Gs_index,2*mseg_recv,MPI_INT,partner,1,&
                       mpi_comm_global,MPI_STATUS_IGNORE, mpierr )

                  ! Compute the local index of each segment.
                  ipcolumn = 1  ! point to Arev2
                  DO i =1, mseg_recv
                     istart = Gs_index(2*i-1)
                     itmp = Gs_index(2*i)-istart+1 ! number of columns of current segment

                     lstart = indxg2l( istart,nblk,mycol,0,npcol )
                     call zlacpy( 'A',na_rows,itmp,Arev2(1,ipcolumn),na_rows,B(1,lstart),ldb )
                     ipcolumn = ipcolumn+itmp
                  END DO

                  ! Send Asend2 first.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart+V_send(i+1)
                  END DO
                  msendT = msend*na_rows
                  call mpi_send( Asend2(1,istart),msendT,MPI_C_DOUBLE_COMPLEX,partner,1,&
                       mpi_comm_global,mpierr )
               
                  ! Send Gs_index2 second.
                  istart = 1
                  DO i=0, partner-1
                     istart = istart + Gp_index2( i+1 )
                  END DO
                  call mpi_send( Gs_index2(2*istart-1),2*Gp_index2(partner+1),MPI_INT,partner,1,&
                       mpi_comm_global,mpierr )
                  
               end if ! ( mycol .lt. partner)
               
            end if ! (mreceive .eq. 0)
            
         end if  ! (msend .eq. 0)

         !********* Current communication is finished ****************
         !**************************************************************************

20       CONTINUE

         ! Delete the matched edges
         call graph_del(nedge, inode, ldm, match, nedgef)
         nedge = nedgef 

         ! Call mpi_barrier
         call mpi_barrier( mpi_comm_global, mpierr  )
         
      END DO ! it=1, npcol-1

      if( nc_local .gt. 0 ) deallocate( Asend2 )
      
      deallocate( V_send, V_receive  )
      deallocate( Gp_recv_index, Gs_index2, Gs_index  )
      deallocate( Gp_index2, Gp_index, Gl_index )
      deallocate( Arev2 )
      deallocate( inode, match)


    end Subroutine OneTo1D_ComplexDataRedistBip
