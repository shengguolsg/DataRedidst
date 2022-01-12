!
      PROGRAM TEST_DataRedist2D
!
        USE MPI
        implicit none
!
!     This routine test PDGEMR2D, and tries to redistribute data from 1D 
!     block column distribution form to 2D block cyclic form. At first,
!     each process has local matrix A of dimension N-by-K, and finally
!     it is transformed to matrix B, distributed with blocksize MB and
!     NB for row and column distribution, respectively. 
!
!     The aim of this routine is to test whether PDGEMR2D can finish
!     this job, and further test its performance when using many
!     processes.         
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, ONE
      PARAMETER          ( ZERO = 0.0D+0, TWO=2.0D+0, ONE=1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CONTEXT, I, IAM, INFO, M, MYCOL, MYROW, N,NB, &
                         NPCOL, NPROCS, NPROW, NZ, myid, NP, J,  &
                         mpierr, na_rows, na_cols, nb_rows, nb_cols
      INTEGER            LWORK, LIWORK, LDB, MB, MA, NA
      DOUBLE PRECISION   time1, time2, time
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. MPI Parameters
      INTEGER            mpi_comm_rows, mpi_comm_cols
!     ..
!     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCB( 50 )
      Complex*16,  ALLOCATABLE :: A(:,:),B(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Xr(:,:)
      character*16        :: arg1, arg2, arg3
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO, &
                         BLACS_SETUP, DESCINIT, PDLAPRNT
!     ..
!     .. Executable Statements ..
!
!     Set up the problem

      call GET_COMMAND_ARGUMENT(1, arg1)
      call GET_COMMAND_ARGUMENT(2, arg2)
      call GET_COMMAND_ARGUMENT(3, arg3)
      read(arg1, *)  N
      read(arg2, *)  NA
      read(arg3, *)  NB

!      N = 30000
!      NA = 19
!      NB = 38 

      M = N
      MA = NA
      MB = NB

!      M = 9273 
!      N = 9273
!
!      MA = 64 
!      NA = 64 
!      NB = 89
!      MB = 89
!
      call MPI_Init( mpierr )
      call MPI_Comm_rank( MPI_COMM_WORLD,myid,mpierr )
      call MPI_Comm_size( MPI_COMM_WORLD,NPROCS,mpierr )

      do nprow = NINT(SQRT(REAL(nprocs))),2,-1
          if(mod(nprocs,nprow) == 0 ) exit
      enddo
      ! at the end of the above loop, nprocs is always divisible by np_cols
      npcol = nprocs/nprow
!
      IF( myid .eq. 0 ) write(*,*) "N,NPROCS,NA,NB=", N,nprocs,NA,NB
!
!     Initialize a single BLACS context
!
      CALL BLACS_GET( 0, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
!
      ! Determine the necessary size of the distributed matrices,
      ! we use the Scalapack tools routine NUMROC for that.

       na_rows = numroc( M, MA, myrow, 0, nprow )
       na_cols = numroc( N, NA, mycol, 0, npcol )

       nb_rows = numroc( M, MB, myrow, 0, nprow )
       nb_cols = numroc( N, NB, mycol, 0, npcol )

       allocate( A(na_rows,na_cols), B(nb_rows,nb_cols))
       allocate( Xr(na_rows,na_cols) )

       call random_number( Xr )
       
       A=Xr+(0.D0,1.0D0)*Xr
!
!     These are basic array descriptors
!
      CALL DESCINIT( DESCB, M, N, MB, NB, 0, 0, CONTEXT, nb_rows, INFO )
      CALL DESCINIT( DESCA, M, N, MA, NA, 0, 0, CONTEXT, na_rows, INFO )
!
!     Redistribute matrix A to matrix B

!     Create mpi_comm_rows and mpi_comm_cols. In ELPA,
!     mpi_comm_rows is used for communicating within rows, i.e., all processes
!     having the same COLUMN coordinate share one mpi_comm_rows. 
      call mpi_comm_split(mpi_comm_world,mycol,myrow,mpi_comm_rows,mpierr)
      call mpi_comm_split(mpi_comm_world,myrow,mycol,mpi_comm_cols,mpierr)
      
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time = MPI_Wtime()
      CALL PZGEMR2D( M,N,A,1,1,DESCA,B,1,1,DESCB,CONTEXT )
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time1 = MPI_Wtime() - time

      IF( myrow.eq.0 .and. mycol.eq.0 ) &
             write(*,*) "BLACS Redistribution costs", time1

      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time = MPI_Wtime()
!      CALL TwoTo2D_DcplxRedistSimple( M,N,A,na_rows,DESCA,B,nb_rows,DESCB, &
!           nprocs, mpi_comm_rows, mpi_comm_cols, mpi_comm_world )
      CALL TwoTo2D_DCplxRedistBLACS( 'R',M,N,A,na_rows,DESCA,B,nb_rows,DESCB  )
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time2 = MPI_Wtime() - time

      IF( myrow.eq.0 .and. mycol.eq.0 ) THEN 
             write(*,*) "Ours redistribution costs", time2
             write(*,*) "Speedup over blacs is ", time1 / time2
      END IF


      deallocate(A)
      deallocate(B)

!      CALL BLACS_GRIDEXIT( CONTEXT )
!
   20 CONTINUE
!
      CALL BLACS_EXIT( 0 )
!
!
      STOP
      
    END PROGRAM TEST_DataRedist2D
