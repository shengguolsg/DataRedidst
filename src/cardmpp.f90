!***********************************************************************

module mod_cardmpp

  interface cardmpp
     module procedure max_match
     module procedure search 
     module procedure blsaug
     module procedure addmem
     module procedure fndpth
     module procedure ancest
     module procedure STATUS
     module procedure addbrg
     module procedure addstk
     module procedure push
     module procedure bastar
     module procedure ssort
     module procedure graph_del
  end interface cardmpp

contains

  subroutine max_match( n,m,edges,ldm, match,nmatch  )

    implicit none

    integer :: n, m, nmatch, ldm
    integer :: edges(ldm,2), match(n)

! There are n vertices and m edges. The vector match(n) gives
! the output of matched edges. 

    ! Local variable
    integer  :: i, j, p, u, v, card, np1, mt2
    
    integer, allocatable :: degree(:), nindex(:), nbrll(:,:), ptr(:), nghbr(:)

    np1 = n+1
    mt2 = m*2
    allocate( degree(n),nbrll(mt2,2),ptr(n),nindex(np1) )
    allocate( nghbr(mt2) ) 
    
!   INPUT EDGES, STORE ADJACENCY LISTS IN NBRLL (NEIGHBOR LINKED LISTS)

    DO  i = 1, n
       degree(i) = 0
       ptr(i) = 0
    END DO
    p = 0
    DO  i = 1, m
       u = edges(i, 1)  ! left vertex
       v = edges(i, 2)  ! right vertex
       
       p = p+1
       nbrll(p,1) = v
       nbrll(p,2) = ptr(u)
       ptr(u) = p
       degree(u) = degree(u) + 1
       p = p + 1
       nbrll(p,1) = u
       nbrll(p,2) = ptr(v)
       ptr(v) = p
       degree(v) = degree(v) + 1
    END DO

!     NOW STORE ADJACENCY LISTS SEQUENTIALLY

    nindex(1) = 1
    DO  i = 1, n
       nindex(i+1) = nindex(i) + degree(i)
    END DO
    DO  i = 1, n
       p = ptr(i)
       DO  j = nindex(i), nindex(i+1)-1
          nghbr(j) = nbrll(p,1)
          p = nbrll(p,2)
       END DO
    END DO

    card = 0
!     CALL SEARCH TO FIND OPTIMAL MATCHING AND CARDINALITY
    CALL search(n,m,degree,nindex,nghbr,match,card)

    nmatch = card
    
!    WRITE(6,60)  CARD
!   60 FORMAT('s  ', I10)

!    DO 68 J = 1 , N
!       WRITE(6,66)  J, MATCH(J)
!66     FORMAT('m  ', 2I10)
!68     CONTINUE

    deallocate( degree,nbrll,ptr,nindex,nghbr )
    
    return
  end subroutine max_match

  subroutine graph_del( m,edges,ldm,match,nedgef )

!*****************************************************************************
!
!! GRAPH_DEL_ARC_MATCH deletes the matched edges from inode and jnode. 
!
!
!  Author:
!    Shengguo Li
!******************************************************************************
    
  implicit none

  integer  ::  m       ! number of edges
  integer  ::  ldm     ! leading dimension of edges
  integer  ::  nedgef  ! The number of remaining edges

  integer  ::  match(*)
  integer  ::  edges(ldm,*)
  
  integer  ::  i, pedge, nstart, nend
  
!!! Search inode and jnode, and delete the matched edges
  nedgef = 0
  pedge = 1
  Do i =1, m
     nstart = edges(i,1)   ! left vertex of current edge for search
     nend   = edges(i,2)   ! right vertex of current edge 

     !! Search the whole matched edges
     if( match(nstart) .ne. nend) then  ! This edge is not matched, and keep it.
        if(pedge .LT. i ) then  ! move i-th entry above, else pedge is in the right position.
           edges(pedge,1) = edges(i,1)
           edges(pedge,2) = edges(i,2)
        end if
        pedge = pedge +1
     !else
        ! current edge is matched and just ignored. 
     end if

  END DO
  nedgef = pedge-1   ! The total number of remaining edges.

  return
end subroutine graph_del

SUBROUTINE search(n,m,degree,nindex,nghbr,match,card )

!     PERFORMS BREADTH-FIRST SEARCHES FROM FREE (UNMATCHED) VERTICES.
!     IDENTIFIES BRIDGES (EDGES JOINING ALTERNATING PATHS FROM
!     FREE VERTICES)
!-----------------------------------------------------------------------
  implicit none

  INTEGER, INTENT(IN)     :: n    ! number of vertices
  INTEGER, INTENT(IN)     :: m    ! number of edges
  INTEGER, INTENT(INOUT)  :: card !
!  
  INTEGER, INTENT(INOUT)  :: degree(*)   ! degree of each vertice
  INTEGER, INTENT(INOUT)  :: nindex(*)   !
  INTEGER, INTENT(INOUT)  :: nghbr(*)    !
  INTEGER, INTENT(INOUT)  :: match(*)    ! 
  
!----------------------------------------------------------
!----------------------------------------------------------  
!     
  INTEGER  :: np1, mt2, infty, undef, nd2
  

!     ARGUMENTS

!     LOCAL VARIABLES

  INTEGER :: b, u, v, e, i, uu, vv, temp, next, nlevi, code
  INTEGER :: vp1, afree, bfree, pfree, dfree, calls
  LOGICAL :: didaug

  integer :: ni, mi

! Local arrays  
  integer, allocatable :: leveli(:), mem(:), evlev(:), odlev(:), bloss(:), &
       anom(:,:), aindex(:), base(:), bstar(:), peakl(:), &
       peakr(:), f(:), path(:), bridge(:,:), pred(:,:), pindex(:), &
       derp(:,:), dindex(:), prdctr(:), edge(:), lr(:)
  integer,dimension(:), allocatable :: bindex
  logical, allocatable          :: visitd(:)
  character(len=1), allocatable :: mark(:)
  
!     COMMON STORAGE

!INTEGER :: ni, mi, match(n), card, nindex(np1), nghbr(mt2)
!COMMON ni, mi, match, card, nindex, nghbr

!CHARACTER (LEN=1) :: mark(n)
!COMMON /blk0/ mark

!INTEGER :: evlev(n), odlev(n), bloss(n)
!COMMON /blk1/ evlev, odlev, bloss

!INTEGER :: anom(2,n), aindex(n), afree
!COMMON /blk2/ anom, aindex, afree

!LOGICAL :: visitd(n)
!COMMON /blk3/ visitd

!INTEGER :: base(nd2), bstar(nd2), peakl(nd2), peakr(nd2), f(n), path(n)
!COMMON /blk4/ base, bstar, peakl, peakr, f, path

!INTEGER :: bridge(3,m), bindex(nd2), bfree
!COMMON /blk5/ bridge, bindex, bfree

!INTEGER :: pred(2,m), pindex(n), pfree, derp(2,m), dindex(n), dfree, prdctr(n)
!COMMON /blk6/ pred, pindex, pfree, derp, dindex, dfree, prdctr

!INTEGER :: edge(mt2)
!COMMON /blk7/ edge

!INTEGER :: lr(n), calls
!COMMON /blk9/ lr, calls
!-----------------------------------------------------------------------

  np1=n+1
  mt2=m*2
  infty=np1
  undef=0
  nd2=n/2

  ni = n
  mi = m

  ! Allocate workspace and variables
  allocate( leveli(n), mem(n) )
  allocate( mark(n) )                         ! blk0
  allocate( evlev(n), odlev(n), bloss(n) )    ! blk1
  allocate( anom(2,n), aindex(n) )            ! blk2
  allocate( visitd(n) )                       ! blk3
  allocate( base(nd2), bstar(nd2), peakl(nd2), peakr(nd2), f(n), path(n) ) ! blk4
  allocate( bridge(3,m), bindex(0:nd2) )        ! blk5
  allocate( pred(2,m), pindex(n), derp(2,m), dindex(n), prdctr(n) )  ! blk6
  allocate( edge(mt2) )                       ! blk7
  allocate( lr(n) )                           ! blk9

!  *** STEP -1 ***  Compute initial matching

!     INITIALIZE MATCH ARRAY

DO  v = 1, ni
  match(v) = 0
  mem(v) = v
END DO
card = 0

CALL ssort(degree, mem, ni)
loop55:  DO  vv = 1, ni
  v = mem(vv)
  IF (match(v) == 0) THEN
    DO  uu = nindex(v), nindex(v+1)-1
      u = nghbr(uu)
      IF (match(u) == 0) THEN
        match(u) = v
        match(v) = u
        card = card + 1
        CYCLE loop55
      END IF
    END DO
  END IF
END DO loop55
IF (card == ni/2) RETURN

! *** STEP 0 *** Initialization

5 DO  v = 1, ni
  evlev(v) = infty
  odlev(v) = infty
  bloss(v) = undef
  pindex(v) = 0
  dindex(v) = 0
  aindex(v) = 0
  visitd(v)   = .false.
  lr(v) = 0
  mark(v) = ' '
  f(v) = 0
  prdctr(v) = 0
END DO
pfree = 1
dfree = 1
afree = 1

DO  e = 1, mi*2
  edge(e) = 0
END DO

DO  v = 0, ni/2
  bindex(v) = 0
END DO
bfree = 1

DO  v = 1, mi-1
  vp1 = v + 1
  pred(1,v) = vp1
  derp(1,v) = vp1
END DO
pred(1, mi) = 0
derp(1, mi) = 0
DO  v = 1, ni - 1
  vp1 = v + 1
  anom(1,v) = vp1
END DO
anom(1, ni) = 0

DO  v = 1, mi-1
  bridge(1,v) = v+1
END DO
bridge(1,mi) = 0

i = -1
didaug = .false.
calls = 0
b = 0

! *** STEP 1 ***  Check for free vertices

DO  v = 1, ni
  IF (match(v) == 0) evlev(v) = 0
END DO

! *** STEP 2 ***  Proceed to next level (I) of search

200 i = i + 1
IF (i > n/2) THEN
  PRINT *, 'SEARCH:  SEARCH LEVEL LIMIT EXCEEDED'
  STOP
END IF
nlevi = 0

!     Loop 210 searches for vertices with level=I

DO  v = 1, ni
  IF ((evlev(v) == i) .OR. (odlev(v) == i)) THEN
!        IF (MIN(EVLEV(V), ODLEV(V)) .EQ. I) THEN
    nlevi = nlevi + 1
    leveli(nlevi) = v
  END IF
END DO

!     If no more vertices have level I then HALT (return to main)

IF (nlevi == 0)  RETURN

!     Is search level I even or odd?

IF (MOD(i,2) == 0) THEN
  
! *** STEP 3 ***  For even search levels
  
  DO  vv = 1, nlevi
    v = leveli(vv)
    IF ((i == 0) .AND. (match(v) /= 0)) CYCLE
    
!           Loop 320 processes all neighbors of current V
    
    DO  uu = nindex(v), nindex(v+1)-1
      u = nghbr(uu)
      
!              EDGE UV IS UNUSED IF CODE IS EVEN
      
      CALL STATUS(u, v, edge, code, 'R', nindex, nghbr)
      IF ((match(u) /= v) .AND. (MOD(code,2) == 0)) THEN
        IF (evlev(u) /= infty) THEN
          temp = (evlev(u) + i) / 2
          IF ((evlev(u) /= i) .OR. (u < v)) THEN
            CALL addbrg(u, v, temp, bridge, bindex, bfree)
          END IF
        ELSE
          IF (odlev(u) == infty) odlev(u) = i + 1
          IF (odlev(u) == i + 1) THEN
            CALL addstk(pred, pindex, pfree, u, v)
            CALL addstk(derp, dindex, dfree, v, u)
            prdctr(u) = prdctr(u) + 1
          ELSE IF (odlev(u) < i) THEN
            CALL addstk(anom, aindex, afree, u, v)
          END IF
        END IF
      END IF
    END DO
  END DO
  
ELSE
  
! *** STEP 4 ***  For odd search levels
  
  DO  vv = 1, nlevi
    v = leveli(vv)
    IF (bloss(v) == undef) THEN
      u = match(v)
      IF (odlev(u) == i) THEN
        IF (u < v) THEN
          CALL addbrg(u, v, i, bridge, bindex, bfree)
        END IF
      ELSE IF (odlev(u) == infty) THEN
        
!                   MAKE V THE ONLY PREDECESSOR OF U
        
        evlev(u) = i + 1
        next = pindex(u)
        
!                   CHECK FOR ERROR CONDITION
        
        IF (next /= 0) PRINT *, 'warning FROM SEARCH: ',  &
            u, ' SHOULD NOT HAVE ANY PREDECESSORS'
        408               IF (next /= 0) THEN
          pindex(u) = pred(1, next)
          pred(1, next) = pfree
          pfree = next
          next = pindex(u)
          GO TO 408
        END IF
        
        CALL addstk(pred, pindex, pfree, u, v)
        CALL addstk(derp, dindex, dfree, v, u)
        prdctr(u) = prdctr(u) + 1
      END IF
    END IF
  END DO
END IF

! *** STEP 5 ***

next = bindex(i)
510 IF (next /= 0) THEN
  u = bridge(2, next)
  v = bridge(3, next)
  next = bridge(1, next)
  IF ((mark(u) == 'e') .OR. (mark(v) == 'e')) GO TO 510
  IF ((bloss(u) /= undef) .AND. (bloss(v) == bloss(u))) GO TO 510
  CALL blsaug(u,v,didaug,b,i,n,m,match,nindex,nghbr,card, &
       mark,evlev,odlev,bloss,anom,aindex, &
       base,bstar,peakl,peakr,f, path,bridge,bindex,bfree, &
       pred,pindex,derp,dindex,prdctr, edge,lr,calls )
  IF (card == ni/2) RETURN
  GO TO 510
END IF
IF (didaug) THEN
  GO TO 5
ELSE
  GO TO 200
END IF

  deallocate( leveli, mem )
  deallocate( mark )                   ! blk0
  deallocate( evlev, odlev, bloss )    ! blk1
  deallocate( anom, aindex )           ! blk2
  deallocate( visitd )                 ! blk3
  deallocate( base,bstar,peakl,peakr,f,path )  ! blk4
  deallocate( bridge, bindex )         ! blk5
  deallocate( pred, pindex, derp, dindex, prdctr )  ! blk6
  deallocate( edge )                   ! blk7
  deallocate( lr )                     ! blk9


END SUBROUTINE search
!***********************************************************************

SUBROUTINE blsaug(w1, w2, didaug, b, i, n,m,match,nindex,nghbr,card, &
       mark,evlev,odlev,bloss,anom,aindex, &
       base,bstar,peakl,peakr,f, path,bridge,bindex,bfree, &
       pred,pindex,derp,dindex,prdctr, edge,lr,calls )

!     FINDS AN AUGMENTING PATH CONTAINING BRIDGE W1-W2 OR CREATES A NEW
!     BLOSSOM.  THE COMMON VARIABLE CALLS IS INCREMENTED EACH TIME THIS
!     ROUTINE IS CALLED.  I IS THE SEARCH LEVEL.  B COUNTS BLOSSOMS.

!     DEPTH-FIRST SEARCHES ARE CONDUCTED SIMULTANEOUSLY FROM W1 AND W2.
!     THE LEFT SEARCH STARTS AT W1 AND THE RIGHT SEARCH STARTS AT W2.
!     LR(N) = +CALLS IF VERTEX N IS VISITED BY THE LEFT SEARCH.
!     LR(N) = -CALLS IF VERTEX N IS VISITED BY THE RIGHT SEARCH.
!-----------------------------------------------------------------------

  implicit none
  
  INTEGER, INTENT(INOUT)                      :: w1
  INTEGER, INTENT(INOUT)                      :: w2
  LOGICAL, INTENT(OUT)                     :: didaug
  INTEGER, INTENT(OUT)                     :: b
  INTEGER, INTENT(IN)                      :: i
!  IMPLICIT CHARACTER (a-z)

  integer            :: n,m,card,bfree,calls
  integer            :: match(*),nindex(*),nghbr(*)
  integer            :: evlev(*), odlev(*), bloss(*)
  integer            :: anom(2,*), aindex(*)
  integer            :: base(*), bstar(*), peakl(*), peakr(*),f(*),path(*)
  integer            :: bridge(3,*),bindex(*)
  integer            :: pred(2,*), pindex(*), derp(2,*), dindex(*), prdctr(*)
  integer            :: edge(*), lr(*)
  character(len=1)   :: mark(*)
  
!     PARAMETERS

  integer :: np1, mt2, undef, nd2

!     ARGUMENTS

!     LOCAL VARIABLES

  INTEGER :: vl, vr, dcv, barier, p1, p2, u, v, temp, next,  &
       INDEX, nmem, uu, lindex, rindex, etop, ljump, rjump
  integer :: ni, mi, itmp, ijob, icode

  integer, allocatable :: member(:), estk(:)
  
!     COMMON STORAGE

!INTEGER :: ni, mi, match(n), card, nindex(np1), nghbr(mt2)
!COMMON ni, mi, match, card, nindex, nghbr

!CHARACTER (LEN=1) :: mark(n)
!COMMON /blk0/ mark

!INTEGER :: evlev(n), odlev(n), bloss(n)
!COMMON /blk1/ evlev, odlev, bloss

!INTEGER :: anom(2,n), aindex(n), afree
!COMMON /blk2/ anom, aindex, afree

!INTEGER :: base(nd2), bstar(nd2), peakl(nd2), peakr(nd2), f(n), path(n)
!COMMON /blk4/ base, bstar, peakl, peakr, f, path

!INTEGER :: bridge(3,m), bindex(nd2), bfree
!COMMON /blk5/ bridge, bindex, bfree

!INTEGER :: pred(2,m), pindex(n), pfree, derp(2,m), dindex(n), dfree, prdctr(n)
!COMMON /blk6/ pred, pindex, pfree, derp, dindex, dfree, prdctr

!INTEGER :: edge(mt2)
!COMMON /blk7/ edge

!INTEGER :: lr(n), calls
!COMMON /blk9/ lr, calls
!-----------------------------------------------------------------------

  calls = calls + 1
  ni = n
  mi = m
  np1=n+1
  mt2=m*2
  undef=0
  nd2=n/2

  allocate( member(n),estk(n) )
  
! *** STEP 0 ***

vl = w1
ljump=0
itmp = 0
IF ( bloss(w1) /= undef) THEN
  CALL bastar(itmp, vl, bstar, f, bloss)
END IF

vr = w2
rjump=0
IF (bloss(w2) /= undef) THEN
  CALL bastar(itmp, vr, bstar, f, bloss)
END IF

!     QUIT IF AN EMPTY BLOSSOM IS DETECTED

IF (vl == vr) GO TO 500

lindex = pindex(vl)
rindex = pindex(vr)
lr(vl) = calls
lr(vr) = -calls
member(1) = vl
member(2) = vr
nmem = 2
f(w1) = undef
dcv = undef
!     IN ORIGINAL ALGORITHM W2 IS REPLACED BY VR BELOW
barier = w2

! *** STEP 1.1 ***  If VL and VR are both free vertices

100 IF ((match(vl) == 0) .AND. (match(vr) == 0) ) THEN
   ijob = -1
  CALL fndpth(w1, vl, undef, ijob, n,m, nindex,nghbr, &
       mark,evlev,odlev,bloss,base,peakl,peakr,f,path,&
       pred,pindex,edge,lr )
  !n,m,match,nindex,nghbr,card, &
  !     mark,evlev,odlev,bloss,anom,aindex,afree, &
  !     base,bstar,peakl,peakr,f, path,bridge,bindex,bfree, &
  !     pred,pindex,pfree,derp,dindex,dfree,prdctr, edge,lr,calls )

  ijob = 1
  CALL fndpth(w2, vr, undef, ijob, n,m,nindex,nghbr, &
       mark,evlev,odlev,bloss,base,peakl,peakr,f,path,&
       pred,pindex,edge,lr )
  
!        Concatenate paths
  
  path(w1) = w2
  
!        AUGMENT MATCHING
  
  p1 = vl
  111    p2 = path(p1)
  match(p1) = p2
  match(p2) = p1
  etop = etop + 1
  estk(etop) = p1
  etop = etop + 1
  estk(etop) = p2
  p1 = path(p2)
  IF (p2 /= vr) GO TO 111
  card = card + 1
  didaug = .true.
  
!        PERFORM TOPOLOGICAL ERASE AND RETURN
  
  122    IF (etop > 0) THEN
    p1 =  estk(etop)
    etop = etop - 1
    IF (mark(p1) /= 'e') THEN
      mark(p1) = 'E'
      next = dindex(p1)
      127          IF (next /= 0) THEN
        p2 = derp(2, next)
        next = derp(1, next)
        prdctr(p2) = prdctr(p2) - 1
        IF (prdctr(p2) == 0) THEN
          etop = etop + 1
          estk(etop) = p2
        END IF
        GO TO 127
      END IF
    END IF
    GO TO 122
  END IF
  GO TO 500
  
! *** STEP 1.2 ***  If VL and VR are not both free
  
ELSE
  IF (MIN(evlev(vl),odlev(vl)) >= MIN(evlev(vr),odlev(vr))) THEN
    GO TO 200
  ELSE
    GO TO 300
  END IF
END IF

! Does VL have unused ancestor edges?

200 IF (lindex == 0) THEN
  u = 0
ELSE
  CALL ancest(vl, 1, u, lindex, pred, edge, mark, nindex, nghbr)
END IF
IF (u == 0) THEN
  
! *** STEP 2.1 ***  VL has no more unused ancestors
  
  IF (f(vl) == undef) THEN
    IF (dcv /= undef) THEN
      GO TO 400
    ELSE
      PRINT *, 'WARNING - BLSAUG(',w1,',',w2,') QUITTING ',  &
          'BECAUSE DCV UNDEFINED AT STEP 2.1'
      GO TO 500
    END IF
  ELSE
    vl = f(vl)
    lindex = pindex(vl)
    GO TO 100
  END IF
ELSE
  
! ***  STEP 2.2 *** VL has some unused ancestors
  
!        MARK VL-U 'USED'
   icode = 1
  CALL STATUS(vl, u, edge, icode, 'W', nindex, nghbr)
  ljump = vl
  IF (bloss(u) /= undef) THEN
    CALL bastar(vl, u, bstar, f, bloss)
    lindex = pindex(vl)
  END IF
  
!    *** Is U marked? ***
  
  IF (ABS(lr(u)) /= calls) THEN
    
!        *** STEP 2.2 (a) ***   U is unmarked
    
    lr(u) = calls
    CALL addmem(member, nmem, u)
    f(u) = vl
    vl = u
    lindex = pindex(vl)
    GO TO 100
  ELSE
    
!        *** STEP 2.2 (b) ***   U is marked
    
    IF ((u == barier) .OR. (u /= vr)) THEN
      GO TO 100
    ELSE
      lr(u) = calls
      CALL addmem(member, nmem, u)
      vr = f(vr)
      IF (rjump /= 0) vr = rjump
      rindex = pindex(vr)
      f(u) = vl
      vl = u
      lindex = pindex(vl)
      dcv = u
      GO TO 100
    END IF
  END IF
END IF

! Does VR have unused ancestor edges?

300 IF (rindex == 0) THEN
  u = 0
ELSE
  CALL ancest(vr, 1, u, rindex, pred, edge, mark, nindex, nghbr)
END IF
IF (u == 0) THEN
  
!    *** STEP 3.1 ***  VR has no more unused ancestors
  
  IF (vr == barier) THEN
    IF (dcv == 0) THEN
      PRINT *, 'WARNING - BLSAUG(',w1,',',w2,') QUITTING ',  &
          'BECAUSE DCV UNDEFINED AT STEP 3.1'
      GO TO 500
    END IF
    vr = dcv
    rindex = pindex(vr)
    barier = dcv
    lr(vr) = -calls
    CALL addmem(member, nmem, vr)
    IF (f(vl) == undef) THEN
      GO TO 400
    ELSE
      vl = f(vl)
      IF (ljump /= 0) vl = ljump
      lindex = pindex(vl)
    END IF
  ELSE
    vr = f(vr)
    rindex = pindex(vr)
  END IF
  GO TO 100
ELSE
  
!     *** STEP 3.2 ***  VR has some unused ancestors
  
!        MARK VR-U 'USED'
   icode = 1
  CALL STATUS(vr, u, edge, icode, 'W', nindex, nghbr)
  rjump = vr
  IF (bloss(u) /= undef) THEN
    CALL bastar(vr, u, bstar, f, bloss)
    rindex = pindex(vr)
  END IF
  
!    *** Is U unmarked? ***
  
  IF (ABS(lr(u)) /= calls) THEN
    
!        *** STEP 3.2 (a) ***
    
    lr(u) = -calls
    CALL addmem(member, nmem, u)
    f(u) = vr
    vr = u
    rindex = pindex(vr)
    GO TO 100
  ELSE
    
!        *** STEP 3.2 (b) ***
    
    IF (u == vl) THEN
      dcv = u
    END IF
    GO TO 100
  END IF
END IF

! *** STEP 4 ***

400 lr(dcv) = 0

!     create new blossom using all vertices marked L or R during call

b = b + 1
IF (b > nd2) THEN
  PRINT *,'BLSAUG: OUT OF STORAGE FOR BASE, BSTAR, PEAKL, PEAKR'
  STOP
END IF
DO  uu = 1, nmem
  u = member(uu)
  IF (u == dcv) CYCLE
  IF (bloss(u) == undef) THEN
    bloss(u) = b
    IF (evlev(u) < odlev(u)) THEN
      odlev(u) = 2 * i + 1 - evlev(u)
    ELSE
      evlev(u) = 2 * i + 1 - odlev(u)
      INDEX = aindex(u)
      421          IF (INDEX /= 0) THEN
        v = anom(2, INDEX)
        INDEX = anom(1, INDEX)
        temp = (evlev(u) + evlev(v))/2
        CALL addbrg(u, v, temp, bridge, bindex, bfree)
        
!                 MARK U-V 'USED'
        icode = 1
        CALL STATUS(u, v, edge, icode, 'W', nindex, nghbr)
        GO TO 421
      END IF
    END IF
  END IF
END DO

peakl(b) = w1
peakr(b) = w2
base(b)  = dcv
bstar(b) = dcv

! *** STEP 5 ***  Return to SEARCH

deallocate( member, estk) 

500 RETURN
END SUBROUTINE blsaug
!***********************************************************************

SUBROUTINE addmem(member, nmem, u)

!     KEEPS TRACK OF VERTICES MARKED L OR R DURING A CALL TO BLSAUG
!     IN CASE A NEW BLOSSOM IS FORMED

!-----------------------------------------------------------------------
  implicit none

  INTEGER, INTENT(OUT)                     :: member(1)
  INTEGER, INTENT(OUT)                     :: nmem
  INTEGER, INTENT(IN)                      :: u

!-----------------------------------------------------------------------
  nmem = nmem + 1
  member(nmem) = u
  RETURN
END SUBROUTINE addmem
!***********************************************************************

SUBROUTINE fndpth(high, low, b, job, n,m, nindex,nghbr, &
       mark,evlev,odlev,bloss,base,peakl,peakr,f,path,&
       pred,pindex,edge,lr )

!     JOB .NE. -1:  RETURNS A LINKED LIST FROM HIGH TO LOW
!     JOB .EQ. -1:  RETURNS A LINKED LIST FROM LOW TO HIGH
!     JOB .EQ.  2:  ALLOWS TRANSFERS BETWEEN LEFT AND RIGHT SEARCHTREES
!                   (OCCURS ONLY WHEN OPENING A BLOSSOM)

!-----------------------------------------------------------------------
  implicit none
  
  INTEGER, INTENT(IN OUT)                  :: high
  INTEGER, INTENT(IN OUT)                  :: low
  INTEGER, INTENT(OUT)                     :: b
  INTEGER, INTENT(OUT)                     :: job
!  IMPLICIT CHARACTER (a-z)

  integer            :: n,m
  integer            :: nindex(*),nghbr(*)
  integer            :: evlev(*), odlev(*), bloss(*)
  integer            :: base(*), peakl(*), peakr(*),f(*),path(*)
  integer            :: pred(2,*), pindex(*)
  integer            :: edge(*), lr(*)
  character(len=1)   :: mark(*)
  logical            :: visitd(n)

!     PARAMETERS

  integer :: np1, nd8, mt2, nd2

!INTEGER, PARAMETER :: n=8192
!INTEGER, PARAMETER :: m=32768
!INTEGER, PARAMETER :: np1=n+1
!INTEGER, PARAMETER :: nd8=n/8
!INTEGER, PARAMETER :: mt2=m*2
!INTEGER, PARAMETER :: nd2=n/2

!     ARGUMENTS


!     LOCAL VARIABLES

  INTEGER :: entrnc, bass, lb, retadd
  INTEGER :: u, v, pree, pntr, succ, vindex

  integer :: stklim, stktop, ni, mi, icode
  integer, allocatable :: stack(:,:) 

!     COMMON STORAGE

!INTEGER :: ni, mi, match(n), card, nindex(np1), nghbr(mt2)
!COMMON ni, mi, match, card, nindex, nghbr

!CHARACTER (LEN=1) :: mark(n)
!COMMON /blk0/ mark

!INTEGER :: evlev(n), odlev(n), bloss(n)
!COMMON /blk1/ evlev, odlev, bloss

!LOGICAL :: visitd(n)
!COMMON /blk3/ visitd

!INTEGER :: base(nd2), bstar(nd2), peakl(nd2), peakr(nd2), f(n), path(n)
!COMMON /blk4/ base, bstar, peakl, peakr, f, path

!INTEGER :: pred(2,m), pindex(n), pfree, derp(2,m), dindex(n), dfree, prdctr(n)
!COMMON /blk6/ pred, pindex, pfree, derp, dindex, dfree, prdctr

!INTEGER :: edge(mt2)
!COMMON /blk7/ edge

!INTEGER :: stklim, stktop
!COMMON /blk8/ stklim, stktop

!INTEGER :: lr(n), calls
!COMMON /blk9/ lr, calls
!-----------------------------------------------------------------------

  ni = n
  mi = m
  np1 = n+1
  nd8 = n/8
  mt2 = m*2
  nd2 = n/2

  stklim = nd8
  stktop = 0

  allocate( stack(8,nd8) ) 

! *** STEP 0 ***

10 IF (high == low) THEN
  
!        WHOLE PATH IS JUST "high"
  
  GO TO 820
END IF

v = high
vindex = pindex(v)

! *** STEP 1 ***

100 IF (vindex == 0) THEN
  u = 0
ELSE
  CALL ancest(v, 2, u, vindex, pred, edge, mark, nindex, nghbr)
END IF
IF (u == 0) THEN
  IF (f(v) == 0) THEN
    PRINT *, 'ERROR --  FNDPTH CANNOT FIND PATH'
    STOP
  ELSE
    v = f(v)
  END IF
  vindex = pindex(v)
  GO TO 100
END IF

! *** STEP 2 ***

IF (bloss(v) ==  b) THEN
  
!        MARK U-V 'VISITED'

   icode = 2
  CALL STATUS(u, v, edge, icode, 'W', nindex, nghbr)
ELSE
  u = base(bloss(v))
END IF

! *** STEP 3 ***

IF (u == low) GO TO 600

! *** STEP 4 ***

IF  (visitd(u)) GO TO 100
IF ( MIN(evlev(u),odlev(u)) <= MIN(evlev(low),odlev(low)) ) GO TO 100
IF (job == 2) GO TO 501
IF ((bloss(u) == b) .AND. (lr(u) == -lr(high))) GO TO 100
! *** STEP 5 ***

501 visitd(u) = .true.
f(u) = v
v = u
vindex = pindex(v)
GO TO 100

! *** STEP 6 ***

600 path(v) = low
620 IF (v /= high) THEN
  u = v
  v = f(v)
  path(v) = u
  GO TO 620
END IF

! *** STEP 7 ***

entrnc = high
710 IF (entrnc /= low) THEN
  bass = path(entrnc)
  IF (bloss(entrnc) /= b) THEN
    
!            SIMULATE CALL TO OPEN(ENTRNC, BASS)
    
    GO TO 900
    
  END IF
777 entrnc = bass
  GO TO 710
END IF

! *** STEP 8 ***

IF (job == -1) THEN
  
!     INVERT PATH IF NECESSARY
  
  pree = 0
  pntr = high
  succ = path(high)
  810 IF (pntr /= low) THEN
    path(pntr) = pree
    pree = pntr
    pntr = succ
    succ = path(pntr)
    GO TO 810
  END IF
  path(pntr) = pree
END IF

!     CHECK STACK BEFORE RETURNING

820 IF (stktop /= 0) THEN
  high   = stack(1, stktop)
  low    = stack(2, stktop)
  b      = stack(3, stktop)
  job    = stack(4, stktop)
  entrnc = stack(5, stktop)
  bass   = stack(6, stktop)
  lb     = stack(7, stktop)
  retadd = stack(8, stktop)
  stktop = stktop - 1
ELSE
  RETURN
END IF

IF (retadd == 777) THEN
  GO TO 777
ELSE IF (retadd == 902) THEN
  GO TO 902
ELSE
  GO TO 904
END IF

!***********************************************************************

!     FUNCTION OPEN STARTS AT LINE 900
!     EMBEDDED INSIDE FNDPTH TO REMOVE RECURSION

900 lb = bloss(entrnc)
IF (evlev(entrnc) <= odlev(entrnc) ) THEN
  CALL push(high, low, b, job, entrnc, bass, lb, 777, stack,stklim, stktop)
  high = entrnc
  low = bass
  job = 2
  b = lb
  GO TO 10
ELSE
  IF (lr(entrnc) > 0) THEN
    CALL push(high, low, b, job, entrnc, bass, lb, 902, stack,stklim,stktop)
    high = peakl(lb)
    low = entrnc
    job = -1
    b = lb
    GO TO 10
902 path(peakl(lb)) = peakr(lb)
    CALL push(high, low, b, job, entrnc, bass, lb, 777, stack,stklim, stktop)
    high= peakr(lb)
    low = bass
    job = 1
    b = lb
    GO TO 10
  ELSE
    CALL push(high, low, b, job, entrnc, bass, lb, 904, stack,stklim, stktop)
    high = peakr(lb)
    low = entrnc
    job = -1
    b = lb
    GO TO 10
904 path(peakr(lb)) = peakl(lb)
    CALL push(high, low, b, job, entrnc, bass,lb,777,stack, stklim,stktop)
    high = peakl(lb)
    low = bass
    job = 1
    b = lb
    GO TO 10
  END IF
END IF

deallocate( stack )

END SUBROUTINE fndpth
!***********************************************************************

SUBROUTINE ancest(v, job, u, INDEX, pred, edge, mark,nindex,nghbr)

!     SEARCHES PREDECESSOR LIST OF V
!         JOB = 1:  CHECK FOR UNUSED EDGES
!         JOB = 2:  CHECK FOR UNVISITED EDGES

!         U .EQ. 0:  NO SUITABLE EDGE WAS FOUND
!         U .NE. 0:  U-V IS A SUITABLE EDGE
!-----------------------------------------------------------------------

  implicit none
  
  INTEGER, INTENT(IN OUT)                  :: v
  INTEGER, INTENT(IN)                  :: job
  INTEGER, INTENT(OUT)                     :: u
  INTEGER, INTENT(OUT)                     :: INDEX
  INTEGER, INTENT(IN)                      :: pred(2,*)
  INTEGER, INTENT(IN OUT)                  :: edge(*)
  CHARACTER (LEN=1), INTENT(IN OUT)        :: mark(*)
  INTEGER, INTENT(IN OUT)                  :: nindex(*)
  INTEGER, INTENT(IN OUT)                  :: nghbr(*)
  !IMPLICIT CHARACTER (a-z)


  INTEGER :: w, code
!-----------------------------------------------------------------------
  u = 0
201 IF ((INDEX /= 0) .AND. (u == 0)) THEN
     w = pred(2, INDEX)
     INDEX = pred(1, INDEX)
     IF (mark(w) /= 'e') THEN
        CALL STATUS(w, v, edge, code, 'R', nindex, nghbr)
        IF (job == 1) THEN
           IF (MOD(code, 2) == 0) u = w
        ELSE
           IF (code < 2) u = w
        END IF
     END IF
     GO TO 201
  END IF
  RETURN
END SUBROUTINE ancest
!***********************************************************************

SUBROUTINE STATUS(u, v, edge, code, job, nindex, nghbr)

!     JOB .EQ. 'W' (WRITE):, ADD CODE TO STATUS(U,V)
!     JOB .EQ. 'R' (READ) :  SETS CODE =  STATUS(U,V)
!     STATUS(U,V) IS SYMMETRIC SO UPPER TRIANGLE IS STORED IN EDGE
!-----------------------------------------------------------------------
  implicit none
  
  INTEGER, INTENT(IN)                      :: u
  INTEGER, INTENT(IN)                      :: v
  INTEGER, INTENT(OUT)                     :: edge(*)
  INTEGER, INTENT(IN OUT)                  :: code
  CHARACTER (LEN=1), INTENT(IN)        :: job
  INTEGER, INTENT(IN)                      :: nindex(*)
  INTEGER, INTENT(IN OUT)                  :: nghbr(*)

  ! IMPLICIT CHARACTER (a-z)

  INTEGER :: i, j, k, INDEX
!-----------------------------------------------------------------------
  IF (u < v) THEN
     i = u
     j = v
  ELSE
     i = v
     j = u
  END IF
  INDEX = 0
  DO  k=nindex(i), nindex(i+1)-1
     IF (nghbr(k) == j) INDEX = k
  END DO
  IF (INDEX == 0) THEN
     PRINT *, 'STATUS:  CANNOT FIND EDGE ', i, j
     STOP
  END IF
  IF (job == 'w') THEN
     edge(INDEX) = edge(INDEX) + code
  ELSE
     code = edge(INDEX)
  END IF
  RETURN
END SUBROUTINE STATUS
!***********************************************************************

SUBROUTINE addbrg(u, v, br, bridge, bindex, bfree)

!     ADDS EDGE (U,V) TO BRIDGE(BR)

!     U, V   -  INTEGER VERTICIES
!     BR     -  INTEGER BRIDGE LEVEL NUMBER
!     BRIDGE -  ARRAY OF SIZE 3 x M  HOLDING ALL BRIDGES
!     BINDEX -  ARRAY OF SIZE N/2 HOLDING POINTERS TO EACH LEVEL
!     BFREE  -  INTEGER POINTING TO FIRST FREE LOCATION IN BRIDGE
!     NEXT   -  INTEGER POINTER

!-----------------------------------------------------------------------
  implicit none
  
  INTEGER, INTENT(IN)                      :: u
  INTEGER, INTENT(IN)                      :: v
  INTEGER, INTENT(IN OUT)                  :: br
  INTEGER, INTENT(IN OUT)                  :: bridge(3,*)
  INTEGER, INTENT(IN OUT)                  :: bindex(*)
  INTEGER, INTENT(IN OUT)                  :: bfree

  !  IMPLICIT CHARACTER (a-z)


  INTEGER :: next
!-----------------------------------------------------------------------
  IF (bfree == 0) THEN
     WRITE(6,20) br, u, v
  20    FORMAT(1X, 'ADDBRG:  STACK OVERFLOW ERROR', 5X, 'BR, U, V:',3I10)
     STOP
  ELSE
     next = bfree
     bfree = bridge(1, next)
     bridge(2, next) = u
     bridge(3, next) = v
     bridge(1, next) = bindex(br)
     bindex(br)       = next
  END IF
  RETURN
END SUBROUTINE addbrg
!***********************************************************************

SUBROUTINE addstk(stack, INDEX, free, u, v)

!     ADDS ELEMENT V TO SUBSTACK U

!     STACK  -  ARRAY OF SIZE 2 x M  HOLDING ALL SUBSTACKS
!     INDEX  -  ARRAY OF SIZE N HOLDING POINTERS TO EACH SUBSTACK
!     FREE   -  INTEGER POINTING TO FIRST FREE LOCATION IN STACK
!     U      -  INTEGER VERTEX
!     V      -  INTEGER VERTEX
!     NEXT   -  INTEGER POINTER
!-----------------------------------------------------------------------

  implicit none
  
  INTEGER, INTENT(IN OUT)                  :: stack(2,*)
  INTEGER, INTENT(IN OUT)                  :: INDEX(*)
  INTEGER, INTENT(IN OUT)                  :: free
  INTEGER, INTENT(IN OUT)                  :: u
  INTEGER, INTENT(IN)                      :: v
!  IMPLICIT CHARACTER (a-z)

  INTEGER :: next
!-----------------------------------------------------------------------
  IF (free == 0) THEN
     WRITE(6,20) u, v
  20    FORMAT(1X, 'ADDSTK:  STACK OVERFLOW ERROR', 5X, 'U, V:',2I10)
     STOP
  ELSE
     next = free
     free = stack(1, next)
     stack(2, next) = v
     stack(1, next) = INDEX(u)
     INDEX(u)       = next
  END IF
  RETURN
END SUBROUTINE addstk
!***********************************************************************

SUBROUTINE push(hi, lo, b, job, ent, bas, lb, ra, stack,stklim,stktop)

!     ADDS ELEMENTS TO STACK USED TO SIMULATE RECURSIVE CALLS
!     BETWEEN FNDPTH AND OPEN
!-----------------------------------------------------------------------

  implicit none
  
  INTEGER, INTENT(IN)                      :: hi
  INTEGER, INTENT(IN)                      :: lo
  INTEGER, INTENT(IN)                      :: b
  INTEGER, INTENT(IN)                      :: job
  INTEGER, INTENT(IN)                      :: ent
  INTEGER, INTENT(IN)                      :: bas
  INTEGER, INTENT(IN)                      :: lb
  INTEGER, INTENT(IN)                      :: ra
  INTEGER, INTENT(OUT)                     :: stack(8,*)
!  IMPLICIT CHARACTER (a-z)

  INTEGER :: stklim, stktop
  !COMMON /blk8/ stklim, stktop
!-----------------------------------------------------------------------
IF (stktop < stklim) THEN
  stktop = stktop + 1
  stack(1, stktop) = hi
  stack(2, stktop) = lo
  stack(3, stktop) = b
  stack(4, stktop) = job
  stack(5, stktop) = ent
  stack(6, stktop) = bas
  stack(7, stktop) = lb
  stack(8, stktop) = ra
ELSE
  WRITE(6, 9) hi, lo, b, job, ra
  9    FORMAT(1X, 'PUSH:  STACK OVERFLOW ERROR',  &
      5X, 'HI, LO, B, JOB, RA =', 5I10)
  STOP
END IF
RETURN
END SUBROUTINE push
!***********************************************************************

SUBROUTINE bastar(v, u, bstar, f, bloss)

!     SETS U = BASE*(U) AND PRESERVES PATH BACK TO OLDV
!-----------------------------------------------------------------------

  implicit none
  
  INTEGER, INTENT(IN OUT)                  :: v
  INTEGER, INTENT(IN OUT)                  :: u
  INTEGER, INTENT(IN OUT)                  :: bstar(*)
  INTEGER, INTENT(OUT)                     :: f(*)
  INTEGER, INTENT(IN OUT)                  :: bloss(*)
!  IMPLICIT CHARACTER (a-z)

  INTEGER :: oldv, w
!-----------------------------------------------------------------------
  oldv = v
10 f(u) = v
  v = u
  u = bstar(bloss(u))
  IF (bloss(u) /= 0) GO TO 10

  w = f(v)
  IF (oldv == 0) THEN
     f(u) = v
     v = 0
  END IF

!     PATH COMPRESSION

20 IF (w /= oldv) THEN
     bstar(bloss(w)) = u
     w = f(w)
     GO TO 20
  END IF
  RETURN
END SUBROUTINE bastar
!***********************************************************************

SUBROUTINE ssort(a,b,l)

!     *      SORTING OF THE VECTOR A(I) IN INCREASING ORDER BY THE     *
!     *      SHELL-ALGORITHM.                                          *
!     *                                                                *
!     *      PARAMETERS:                                               *
!     *      INPUT:                                                    *
!     *         L        DIMENSION OF THE VECTOR                       *
!     *         A(I)     VECTOR TO SORT (INTEGER)                      *
!     *         B(I)     = I    (INTEGER)  I=1,...,N                   *
!     *      OUTPUT:                                                   *
!     *         A(I)     THE SORTED VECTOR                             *
!     *         B(I)     PERMUTATION VECTOR OF THE SORTED VECTOR       *
!     *                                                                *
! *** ******************************************************************
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none
  
  INTEGER, INTENT(IN OUT)                  :: a(*)
  INTEGER, INTENT(IN OUT)                  :: b(*)
  INTEGER, INTENT(IN)                      :: l
  
  !IMPLICIT INTEGER (a-z)

  integer  :: f, n2, s, t, ls, is, ah, bh, js, i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
f = 1
IF (l <= f) RETURN
n2 = (l - f+1)/2
s = 1023
DO  t = 1, 10
  IF (s > n2) GO TO 90
  ls = l - s
  DO  i = f, ls
    is = i + s
    ah = a(is)
    bh = b(is)
    j = i
    js = is
    5 IF (ah >= a(j)) GO TO 10
    a(js) = a(j)
    b(js) = b(j)
    js = j
    j = j - s
    IF (j >= f) GO TO 5
    10 a(js) = ah
    b(js) = bh
  END DO
  90 s = s/2
END DO
RETURN
END SUBROUTINE ssort

end module mod_cardmpp
