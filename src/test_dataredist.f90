!
      PROGRAM TEST_DataRedist
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
                         mpierr, na_rows, na_cols, NPROWA, NPCOLA, &
                         CONTEXTA, nb_rows, nb_cols, MYROWA, MYCOLA
      INTEGER            LWORK, LIWORK, LDB, MB, K, MA, NA
      DOUBLE PRECISION   time1, time2
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCB( 50 )
      DOUBLE PRECISION,  ALLOCATABLE :: A(:,:),B(:,:)
      character*16       :: arg1
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
      read(arg1, *)  N

!      N = 32768 
      K = N 
!
      NB = 64
      MB = 64
	  
      call MPI_Init( mpierr )
      call MPI_Comm_rank( MPI_COMM_WORLD,myid,mpierr )
      call MPI_Comm_size( MPI_COMM_WORLD,NPROCS,mpierr )
	  

      do nprow = NINT(SQRT(REAL(nprocs))),2,-1
          if(mod(nprocs,nprow) == 0 ) exit
      enddo
      ! at the end of the above loop, nprocs is always divisible by np_cols
      npcol = nprocs/nprow
!
!      IF( myid .eq. 0 ) write(*,*) "test_pdstedc", MPI_COMM_WORLD
!
      MA = N
      NA = ceiling( real(k)/real(nprocs) ) 
!
!     Initialize a single BLACS context
!
      CALL BLACS_GET( 0, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
!
      NPROWA = 1
      NPCOLA = NPROCS
      CALL BLACS_GET( CONTEXT, 0, CONTEXTA )
      CALL BLACS_GRIDINIT( CONTEXTA, 'R', NPROWA, NPCOLA )
      CALL BLACS_GRIDINFO( CONTEXTA, NPROWA, NPCOLA, MYROWA, MYCOLA )
!      
      ! Determine the necessary size of the distributed matrices,
      ! we use the Scalapack tools routine NUMROC for that.

       na_rows = numroc( N, MA, myrowa, 0, nprowA )
       na_cols = numroc( N, NA, mycola, 0, npcolA )

       nb_rows = numroc( N, MB, myrow, 0, nprow )
       nb_cols = numroc( N, NB, mycol, 0, npcol )

      allocate( A(na_rows,na_cols), B(nb_rows,nb_cols))      

      call random_number(A)
!
!     These are basic array descriptors
!
      CALL DESCINIT( DESCB, N, N, MB, NB, 0, 0, CONTEXT, nb_rows, INFO )
      CALL DESCINIT( DESCA, N, N, MA, NA, 0, 0, CONTEXTA,na_rows, INFO )
!      DESCA(1) = 502       ! DTYPE
!      DESCA(2) = CONTEXTA   ! CTXT
!      DESCA(3) = N         ! N_A
!      DESCA(4) = NA        ! NB
!      DESCA(5) = 0         ! SRC_A
!      DESCA(6) = MA        ! LLD_A

!     Redistribute matrix A to matrix B
   
      if(myid.eq.0) write(*,*) "N, nprocs, msi=", N, nprocs, NB

      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time1 = MPI_Wtime()
      CALL PDGEMR2D( N,N,A,1,1,DESCA,B,1,1,DESCB,CONTEXT  ) 
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      time2 = MPI_Wtime() - time1

      IF( myrow.eq.0 .and. mycol.eq.0 ) &
             write(*,*) "Redistribution costs", time2


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
      END
