!
      PROGRAM TEST_DataRedist1D
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
!     ..
!     .. MPI Parameters
      INTEGER            mpi_comm_rows, mpi_comm_cols
      
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. Local Arrays ..
      INTEGER            DESCA( 10 ), DESCC( 10 ), descb(9) 
      INTEGER,ALLOCATABLE :: Ms(:)                       ! Number of columns of each process
      Complex*16,  ALLOCATABLE :: A(:,:),B(:,:),C(:,:)
      Double precision, allocatable :: Xr(:,:)
      character*16  :: arg1
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
!      call GET_COMMAND_ARGUMENT(1, arg1)
!      read(arg1, *)  file_name

      N = 256
!
      MB = N
      NB = 64 
!
      call MPI_Init( mpierr )
      call MPI_Comm_rank( MPI_COMM_WORLD,myid,mpierr )
      call MPI_Comm_size( MPI_COMM_WORLD,NPROCS,mpierr )

!      do nprow = NINT(SQRT(REAL(nprocs))),2,-1
!          if(mod(nprocs,nprow) == 0 ) exit
!      enddo
      ! at the end of the above loop, nprocs is always divisible by np_cols
!      npcol = nprocs/nprow
!
!     Initialize a single BLACS context
!
      nprow = 1
      npcol = nprocs
      CALL BLACS_GET( 0, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )

!      write(*,*) "rank,myrow,mycol=",myid,myrow,mycol
!
      allocate( Ms(NPROCS) ) ! No. of columns of each process

      ! Process 0 read Ms from a file, and then broadcast to others.
      IF( myid .eq. 0 ) THEN
!         file_name = 'bnd_mult.txt'

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

      if(myid .eq. 0) write(*,*) "nev, n=", nev, n

      ! Determine the necessary size of the distributed matrices,
      ! we use the Scalapack tools routine NUMROC for that.
      nb_rows = numroc( N, MB, myrow, 0, nprow )
      nb_cols = numroc( nev, NB, mycol, 0, npcol )
!      
      CALL DESCINIT( DESCC, N, N, MB, NB, 0, 0, CONTEXT, nb_rows, INFO )
      CALL DESCINIT( DESCB, N, nev, MB, NB, 0, 0, CONTEXT, N, INFO )
      
!      write(*,*) "rank,nb_rows,nb_cols=", myid, nb_rows, nb_cols, nprocs
!      
      Msi = Ms(myid+1) 
      allocate( A(N,Msi),B(n,nb_cols),C(n,nb_cols)  )
      allocate( Xr(N,Msi) )
      
      call random_number(Xr)

      A=Xr+(0.D0,1.0D0)*Xr

!      if (myid .eq. 3 ) then
!         Do j=1, msi
!            DO i=1,n
!               write(*,*) "A(",i,j,")=", A(i,j)
!            end do
!         END DO
!      end if

!     Call the routine 1DTo2D_DataRedist to obtain the BCDD form

      call mpi_barrier( mpi_comm_world, mpierr )
      time = MPI_Wtime()
      
      Call OneTo1D_ComplexDataRedistBip( N, A, N, msi,B, N, nb_cols, NB, Ms, nprocs, &
                                      mpi_comm_world )
      
      call mpi_barrier( mpi_comm_world, mpierr )
      time1 = MPI_Wtime() - time

      if( myid .eq. 0 ) write(*,*) "Data redistribution costs", time1

!      call PZLAPRNT( n,nev,B,1,1,descB,1,1,'A0',6, C )

!      if (myid .eq. 3 ) then
!         Do j=1, nb_cols
!            DO i=1,n
!               write(*,*) "B(",i,j,")=", B(i,j)
!           end do
!         END DO
!      end if
      
      deallocate(A)
      deallocate(B)
      deallocate(C)

!      CALL BLACS_GRIDEXIT( CONTEXT )
!
   20 CONTINUE
!
      CALL BLACS_EXIT( 0 )
!
!
      STOP
      
    END PROGRAM TEST_DataRedist1D
