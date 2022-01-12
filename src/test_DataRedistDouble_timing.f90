!
      PROGRAM TEST_DataRedistDouble_Timing
!
        USE MPI
        implicit none
!
!     This routine tests our new data redistribution algorithm, and write the basic 
!     functionalities. Later, we will transform it to a callable subroutine.
!        
!     This routine is to be incorporate with ELPA. It uses row and column communicators.  
!

        
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, ONE
      PARAMETER          ( ZERO = 0.0D+0, TWO=2.0D+0, ONE=1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CONTEXT, I, IAM, INFO, M, MYCOL, MYROW, N,NB, &
                         NPCOL, NPROCS, NPROW, NZ, myid, NP, J,  &
                         mpierr, na_rows, na_cols, nb_rows, nb_cols 
      INTEGER            LWORK, LIWORK, LDB, MB, Msi, nev
      DOUBLE PRECISION   time1, time
      CHARACTER(len=256) file_name
      CHARACTER(len=1)   rowcol
      character*16        :: arg1
      character*16        :: arg2
!     ..
!     .. MPI Parameters
      INTEGER            mpi_comm_rows, mpi_comm_cols
      
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. Local Arrays ..
      INTEGER            DESCA( 10 ), DESCB( 10 )
	  DOUBLE PRECISION    :: Rtimes(4)
      INTEGER,ALLOCATABLE :: Ms(:)                       ! Number of columns of each process
      Double precision,  ALLOCATABLE :: A(:,:),B(:,:),Z(:,:)
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO, &
                         BLACS_SETUP, DESCINIT, PDLAPRNT
!     ..
!     .. Executable Statements ..
!
!     Set up the problem
!
      call GET_COMMAND_ARGUMENT(1, arg1)
!      call GET_COMMAND_ARGUMENT(2, arg2)
!      read(arg1, *)  file_name
      read(arg1, *)  N

      !N = 16384 
!
      MB = 64 
      NB = 64
      rowcol = 'C'
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
!     Initialize a single BLACS context
!
      CALL BLACS_GET( 0, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, rowcol, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )

!      write(*,*) "rank,myrow,mycol=",myid,myrow,mycol
!
!     Create mpi_comm_rows and mpi_comm_cols. In ELPA,
!     mpi_comm_rows is used for communicating within rows, i.e., all processes
!     having the same COLUMN coordinate share one mpi_comm_rows. 
      call mpi_comm_split(mpi_comm_world,mycol,myrow,mpi_comm_rows,mpierr)
      call mpi_comm_split(mpi_comm_world,myrow,mycol,mpi_comm_cols,mpierr)
!
      allocate( Ms(NPROCS) ) ! No. of columns of each process

      ! Process 0 read Ms from a file, and then broadcast to others.
      IF( myid .eq. 0 ) THEN
!         file_name = 'output10k_ex1.txt'

!         open(unit=1, file=file_name )
!         read(1,*) ( Ms(i), i=1,nprocs )
!         close(1)
          Ms = N/nprocs
      END IF

      CALL MPI_BCAST( Ms,nprocs,MPI_INT,0,MPI_COMM_WORLD,mpierr )

      nev = 0
      DO i=1, nprocs
         nev = nev + Ms(i)
      END DO

      ! Determine the necessary size of the distributed matrices,
      ! we use the Scalapack tools routine NUMROC for that.
      nb_rows = numroc( N, MB, myrow, 0, nprow )
      nb_cols = numroc( nev, NB, mycol, 0, npcol )
!      
      CALL DESCINIT( DESCB, N, N, MB, NB, 0, 0, CONTEXT, nb_rows, INFO )
!      
      Msi = Ms(myid+1) 
      allocate( A(N,Msi), B(nb_rows,nb_cols), Z(nb_rows,nb_cols)  )

      if(myid.eq.0) write(*,*) "N, Nprocs, msi=", N, nprocs, Msi
      
      call random_number(A)

!     Call the routine 1DTo2D_DataRedist to obtain the BCDD form

      call mpi_barrier( mpi_comm_world, mpierr )
      time = MPI_Wtime()
      
      Call OneTo2D_DoubleDataRedist5( rowcol, N, A, N, msi,B, nb_rows, nb_cols, MB, NB, DESCB, Ms, nprocs, &
           mpi_comm_rows, mpi_comm_cols, mpi_comm_world, Rtimes )
      
      call mpi_barrier( mpi_comm_world, mpierr )
      time1 = MPI_Wtime() - time

      if( myid .eq. 0 ) then
		write(*,*) "Data redistribution costs", time1
		write(*,*) "Index computation costs", Rtimes(1)
		write(*,*) "Message packing and unpacking costs", Rtimes(2)
		write(*,*) "Data transfer costs", Rtimes(4)
	  end if

!      call PZLAPRNT( n,nev,B,1,1,descB,1,1,'A0',6, Z )
      
      deallocate( A )
      deallocate( B )

!      CALL BLACS_GRIDEXIT( CONTEXT )
!
   20 CONTINUE
!
      CALL BLACS_EXIT( 0 )
!
!
      STOP
      
    END PROGRAM TEST_DataRedistDouble_Timing
