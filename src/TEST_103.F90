subroutine rpn_mpi_test_103
!
! simple, standalone test for RPN_MPI_ez_halo_parms
! does not use the full RPN_MPI setup
!
  use ISO_C_BINDING
  implicit none
#include <RPN_MPI.hf>
  include 'RPN_MPI_mpif.inc'     ! interface prototypes for MPI routines

  type(RPN_MPI_mpi_definitions)     :: d
  type(RPN_MPI_mpi_definitions_raw) :: dr
  type(RPN_MPI_Comm) :: col_comm, row_comm
  type(RPN_MPI_Loc) :: p
  integer, parameter :: NXCH = 100
  integer :: NI, NJ, NK, halox, haloy
  integer, dimension(:,:,:), allocatable :: z
  integer :: rankx, sizex, ranky, sizey, ier, petot, ranktot, errors
  integer :: i, j, k, offx, offy, l, m, larg1, larg2, larg3, stat1, stat2, stat3
!   integer :: rowcomm, colcomm
  integer :: ilo, ihi, jlo, jhi
  character(len=128) :: argv1, argv2, argv3, mode
  logical :: printit, redblack, yfirst, async, barrier
  real(kind=8) :: t1, t2
  real(kind=8), dimension(NXCH) :: txch

  call MPI_Init(ier)
  call RPN_MPI_get_mpi_definitions(d, ier)       ! get wrapped definitions (for calls to MPI routines)
  call RPN_MPI_get_mpi_definitions_raw(dr, ier)  ! get "raw" definitions
  printit = .false.
  call get_command_argument(1,argv1,larg1,stat1)
  call get_command_argument(2,argv2,larg2,stat2)
  call get_command_argument(3,argv3,larg3,stat3)
  if(stat1 .ne. 0 .or. stat2 .ne. 0) goto 777
  read(argv1,*,err=777)sizex,sizey
  read(argv2,*,err=777)NI, NJ, NK, halox, haloy
  mode = 'DEFAULT'
  yfirst = .false.
  if(stat3 .eq. 0) then
    printit  = argv3(1:1) .eq. 't'  ! print arrays
!     redblack = argv3(2:2) .eq. 'r'  ! use red/black method
!     async    = argv3(2:2) .eq. 'a'  ! use async method
    barrier  = argv3(2:2) .eq. 'b'  ! activate barrier between E-W and N-S
    yfirst   = argv3(3:3) .eq. 'y'  ! y first, then x
  else
    argv3 = 't'  ! this will activate the rnak printout
  endif
!   if(redblack) mode = 'REDBLACK'
!   if(async)    mode = 'ASYNC'
  if(barrier)  mode = 'BARRIER'
  printit = printit .and. ni*nj < 30

  allocate (z(1-halox:NI+halox,1-haloy:NJ+haloy,NK))
!
! split d%MPI_COMM_WORLD into row and column
!
  call MPI_comm_size(d%MPI_COMM_WORLD, petot, ier)
  call MPI_comm_rank(d%MPI_COMM_WORLD, ranktot, ier)
  if(sizex*sizey .ne. petot) goto 777
!   if(ranktot == 0) write(6,*)'redblack =',redblack
  if(ranktot == 0) write(6,2)'sizex,sizey,NI,NJ,NK,halox,haloy=',sizex,sizey,NI,NJ,NK,halox,haloy
! storage_size(d) not equal to storage_size(dr) is an ERROR
  if(ranktot == 0) write(6,2)'storage size of RPN_MPI_mpi_definitions and RPN_MPI_mpi_definitions =', &
       storage_size(d), storage_size(dr)
  if(storage_size(d) .ne. storage_size(dr)) then
    write(6,*) 'ERROR: storage size of wrapped types is not equal to storage size of "raw" types'
    goto 777
  endif
  if(yfirst) then
    rankx = ranktot / sizey
    ranky = mod(ranktot,sizey)
  else
    ranky = ranktot / sizex
    rankx = mod(ranktot,sizex)
  endif
  do l=0,petot-1
    if(ranktot .eq. l .and. argv3(1:1) .eq. 't') write(6,2)'ranktot, petot, rankx, ranky =',ranktot+1, petot, rankx, ranky
    call MPI_Barrier(d%MPI_COMM_WORLD,ier)
  enddo
  call MPI_Comm_split(d%MPI_COMM_WORLD, ranky, ranktot, row_comm,ier)
  call MPI_Comm_split(d%MPI_COMM_WORLD, rankx, ranktot, col_comm,ier)
! set halo exchange timing parameters
  call RPN_MPI_ez_halo_parms(row_comm, col_comm, mode)

  offy = ranky * NJ
  offx = rankx * NI
  z = 0
  do k=1,NK
  do j=1,NJ
  do i=1,NI
    z(i,j,k) = 1000000*mod(i+offx,256) + 1000*mod(j+offy,256) + k
  enddo
  enddo
  enddo
  call MPI_Barrier(d%MPI_COMM_WORLD,ier)
!
! if we use
!      include 'RPN_MPI_mpi_definitions.inc'
!      call RPN_MPI_init(.......) 
!      type(mpi_layout) :: l
!      call RPN_MPI_get_mpi_layout(l, ierr)
! then we can use
!      rowcomm = p%comm%grid%row
!      colcomm = p%comm%grid%column
!
!
  p%a = loc(z)                             ! address of array subject to halo exchange
!   row_comm = RPN_MPI_Comm(rowcomm)
!   col_comm = RPN_MPI_Comm(colcomm)
  call RPN_MPI_ez_halo(LoC(z),1-halox,NI+halox,1-haloy,NJ+haloy,NI,NJ,NK,halox,haloy)
  call RPN_MPI_halo(LoC(z),1-halox,NI+halox,1-haloy,NJ+haloy,NI,NJ,NK,halox,haloy,row_comm,col_comm)
  call RPN_MPI_reset_halo_timings          ! ignore timings for first call
  call MPI_Barrier(d%MPI_COMM_WORLD,ier)     ! full sync

  do i = 1, NXCH
    call MPI_Barrier(d%MPI_COMM_WORLD,ier)
    t1 = MPI_Wtime()
    p%a = loc(z)
    call RPN_MPI_halo(LoC(z),1-halox,NI+halox,1-haloy,NJ+haloy,NI,NJ,NK,halox,haloy,row_comm,col_comm)
    txch(i) = MPI_Wtime() - t1
  enddo
  txch = txch * 1000000   ! convert into microseconds
  call MPI_Barrier(d%MPI_COMM_WORLD,ier)

  call RPN_MPI_print_halo_timings

  call MPI_Barrier(d%MPI_COMM_WORLD,ier)
  if(printit) then                       ! if we are printing arrays after exchange
    do m = sizey-1, 0, -1
    do l = sizex-1, 0, -1
      if(ranktot .eq. l+m*sizex)then
	write(6,*)"PE(",rankx,',',ranky,')'
	do j=NJ+haloy,1-haloy,-1
	  print 1,z(:,j,1)
	enddo
      endif
      call MPI_Barrier(d%MPI_COMM_WORLD,ier)
    enddo
    enddo
  endif
!
! error analysis
!
  errors = 0
  ilo = 1 - halox
  jlo = 1 - haloy
  ihi = NI + halox
  jhi = NJ + haloy
  if(ranky .eq. 0) jlo = 1         ! south PE
  if(ranky .eq. sizey-1) jhi = NJ  ! north PE
  if(rankx .eq. 0) ilo = 1         ! west PE
  if(rankx .eq. sizex-1) ihi = NI  ! east PE
  do k = 1, nk
  do j = jlo, jhi
  do i = ilo, ihi
    if(z(i,j,k) .ne. 1000000*mod(i+offx,256) + 1000*mod(j+offy,256) + k) errors = errors + 1
  enddo
  enddo
  enddo
  call flush(6)
  call MPI_Barrier(d%MPI_COMM_WORLD,ier)
  if(ranktot .eq. 0) write(6,*)'avg/min/max per exchange =',sum(txch)/NXCH,minval(txch),maxval(txch),' microseconds'
  if(errors .ne. 0)  write(6,*)'rank, ERRORS =',ranktot, errors
!   write(6,*)'rank, ERRORS =',ranktot, errors
777 continue
  call MPI_Finalize(ier)
1 format(20I10.9)
2 format(A,20i6)
end subroutine rpn_mpi_test_103
