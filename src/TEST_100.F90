subroutine rpn_mpi_test_100
! check that the RPN_MPI includes behave correctly
! especially when adding the MPI prototypes using RPN_MPI wrappers
  use ISO_C_BINDING
  implicit none
#include <RPN_MPI.hf>
#include <RPN_MPI_mpif.hf>
  print *,'INFO: compilation/load is O.K.'
  print *,'INFO: this is a non MPI test'
  call sub1          ! set values
  call sub2          ! check values from module
  call sub3          ! check values in "user" mode
  print *,'INFO: if no error message was issued, the test is O.K.'
end subroutine rpn_mpi_test_100
!============================================================================
! this test routine uses an internal module
! it initializes some variables in that module 
! in order to make the user mode test possible
!============================================================================
subroutine sub1
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
  type(RPN_MPI_Comm) :: dcomm = RPN_MPI_Comm(12345)

  call RPN_MPI_init_mpi_layout         ! get wrapped structure pointers mw and dw initialized

  mw%comm%wrld%all = dw%MPI_COMM_WORLD
  mw%rank%wrld%all = 1
  mw%size%wrld%all = 111
  mw%size%sgrd%all = 55
  mw%size%grid%all = 11
  mw%comm%grid%all = dcomm
  mw%comm%sgrd%all = RPN_MPI_Comm(23456)
  mw%comm%blck%all = RPN_MPI_Comm(34567)

  lw%wrld%all%comm = dw%MPI_COMM_WORLD
  lw%wrld%all%rank = 1
  lw%wrld%all%size = 111
  lw%sgrd%all%size = 55
  lw%grid%all%size = 11
  lw%grid%all%comm = dcomm
  lw%sgrd%all%comm = RPN_MPI_Comm(23456)
  lw%blck%all%comm = RPN_MPI_Comm(34567)

end subroutine sub1
!============================================================================
!  example of "library" code, uses internal module RPN_MPI_mpi_layout
!  to modify some of the internal variables
!============================================================================
subroutine sub2
! check that wrapped and raw structures contain the dame information
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none

  print *,'==================== FROM MODULE  ===================='
  print *,'MPI_COMM_WORLD    =', MPI_COMM_WORLD
  print *,'dr%MPI_COMM_WORLD =', dr%MPI_COMM_WORLD
  if(dr%MPI_COMM_WORLD .ne. MPI_COMM_WORLD) print *,'ERROR'
  print *,'dw%MPI_COMM_WORLD =', dw%MPI_COMM_WORLD
  if(transfer(dw%MPI_COMM_WORLD,1) .ne. MPI_COMM_WORLD) print *,'ERROR'

  print *,'MPI_COMM_NULL     =', MPI_COMM_NULL
  print *,'dr%MPI_COMM_NULL  =', dr%MPI_COMM_NULL
  if(dr%MPI_COMM_NULL .ne. MPI_COMM_NULL) print *,'ERROR'
  print *,'dw%MPI_COMM_NULL  =', dw%MPI_COMM_NULL
  if(dw%MPI_COMM_NULL%wrapped_value .ne. MPI_COMM_NULL) print *,'ERROR'

! NOTE: 
!   not all compilers behave correctly with print *,mw%comm%wrld%all
!   hence the use of transfer(mw%comm%wrld%all,1) to force a "honest" integer
! mw, ml : OLD style 
! lw, lr : NEW style 
  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, transfer(mw%comm%wrld%all,1)
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, transfer(mw%comm%grid%all,1)
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, transfer(mw%comm%sgrd%all,1)
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, transfer(mw%comm%appl%all,1)
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, transfer(mw%comm%blck%all,1)
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

  print *,'lw%wrld%all%comm  =', lr%wrld%all%comm, transfer(lw%wrld%all%comm,1)
  print *,'lw%grid%all%comm  =', lr%grid%all%comm, transfer(lw%grid%all%comm,1)
  print *,'lw%sgrd%all%comm  =', lr%sgrd%all%comm, transfer(lw%sgrd%all%comm,1)
  if(transfer(lw%sgrd%all%comm,1) .ne. 23456) print *,'ERROR'
  print *,'lw%appl%all%comm  =', lr%appl%all%comm, transfer(lw%appl%all%comm,1)
  print *,'lw%blck%all%comm  =', lr%blck%all%comm, transfer(lw%blck%all%comm,1)
  if(transfer(lw%blck%all%comm,1) .ne. 34567) print *,'ERROR'

  print *,'mw%rank%wrld%all  =', ml%rank%wrld%all, transfer(mw%rank%wrld%all,1)
  if(transfer(mw%rank%wrld%all,1) .ne. 1) print *,'ERROR'
  print *,'mw%size%wrld%all  =', ml%size%wrld%all, transfer(mw%size%wrld%all,1)
  if(transfer(mw%size%wrld%all,1) .ne. 111) print *,'ERROR'
  print *,'mw%rank%sgrd%all  =', ml%rank%sgrd%all, transfer(mw%rank%sgrd%all,1)
  print *,'mw%size%sgrd%all  =', ml%size%sgrd%all, transfer(mw%size%sgrd%all,1)
  if(transfer(mw%size%sgrd%all,1) .ne. 55) print *,'ERROR'
  print *,'mw%rank%grid%all  =', ml%rank%grid%all, transfer(mw%rank%grid%all,1)
  print *,'mw%size%grid%all  =', ml%size%grid%all, transfer(mw%size%grid%all,1)
  if(transfer(mw%size%grid%all,1) .ne. 11) print *,'ERROR'

  print *,'lw%wrld%all%rank  =', lr%wrld%all%rank, transfer(lw%wrld%all%rank,1)
  if(transfer(lw%wrld%all%rank,1) .ne. 1) print *,'ERROR'
  print *,'lw%wrld%all%size  =', lr%wrld%all%size, transfer(lw%wrld%all%size,1)
  if(transfer(lw%wrld%all%size,1) .ne. 111) print *,'ERROR'
  print *,'lw%sgrd%all%rank  =', lr%sgrd%all%rank, transfer(lw%sgrd%all%rank,1)
  print *,'lw%sgrd%all%size  =', lr%sgrd%all%size, transfer(lw%sgrd%all%size,1)
  if(transfer(lw%sgrd%all%size,1) .ne. 55) print *,'ERROR'
  print *,'lw%grid%all%rank  =', lr%grid%all%rank, transfer(lw%grid%all%rank,1)
  print *,'lw%grid%all%size  =', lr%grid%all%size, transfer(lw%grid%all%size,1)
  if(transfer(lw%grid%all%size,1) .ne. 11) print *,'ERROR'

end subroutine sub2

!============================================================================
!  example of "user" code
!  uses the access functions,
!  the prototypes of which are obtained from include files
!  RPN_MPI.hf        for the RPN_MPI library prototypes
!  RPN_MPI_mpif.hf  for the MPI function prototypes using RPN_MPI wrappers
!============================================================================
subroutine sub3
! perform the same check as sub2 but with copies of structures
! obtained using RPN_MPI_get_mpi_.... subroutines
  use ISO_C_BINDING
  implicit none
#include <RPN_MPI.hf>
#include <RPN_MPI_mpif.hf>
  type(RPN_MPI_mpi_definitions_raw) :: dr  ! "raw" MPI symbols
  type(RPN_MPI_mpi_definitions)     :: dw  ! "wrapped" MPI symbols
  type(mpi_layout_internal)         :: ml  ! "raw" RPN_MPI layout information
  type(mpi_layout)                  :: mw  ! "wrapped" RPN_MPI layout information
  type(mpi_layout_r)                :: lr  ! "raw" RPN_MPI layout information (replaces ml)
  type(mpi_layout_f)                :: lw  ! "wrapped" RPN_MPI layout information (replaces mw)
  integer :: ier, count, dest, tag
  logical :: false = .false.
  integer(C_INT), dimension(1,2,3), target :: array
  type(RPN_MPI_Loc)  :: buf        ! the address of any array
  type(RPN_MPI_Comm) :: comm       ! a communicator
  integer, dimension(:), allocatable :: mpistatus

! get the RPN_MPI "mirror" definitions for subsequent pure MPI calls
! get a copy of the "raw", everything an integer, structure
  call RPN_MPI_get_mpi_definitions_raw(dr, ier)  
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_definitions_raw)'
! get a copy of the "typed" definitions 
  call RPN_MPI_get_mpi_definitions(dw, ier)      
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_definitions)'

! get the RPN_MPI layout (application/supergrid/grid ...) information
!
! get a copy of the "raw", everything an integer, OLD structure
  call RPN_MPI_get_mpi_layout_raw(ml, ier)  ! 
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout_raw[ml])'
! get a copy of the "raw", everything an integer, NEW structure
  call RPN_MPI_get_mpi_layout_raw(lr, ier)  
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout_raw[lr])'

! get a copy of the "typed" OLO definitions 
  call RPN_MPI_get_mpi_layout(mw, ier)      
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout[mw])'
! get a copy of the "typed" NEW definitions 
  call RPN_MPI_get_mpi_layout(lw, ier)      
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout[lw])'

  print *,'==================== IN USER MODE ===================='
  print *,'MPI_COMM_WORLD    =', transfer(dw%MPI_COMM_WORLD,1), dr%MPI_COMM_WORLD
  if(transfer(dw%MPI_COMM_WORLD,1) .ne. dr%MPI_COMM_WORLD) print *,'ERROR'

  print *,'MPI_COMM_NULL     =', transfer(dw%MPI_COMM_NULL,1), dr%MPI_COMM_NULL
  if(transfer(dw%MPI_COMM_NULL,1) .ne. dr%MPI_COMM_NULL) print *,'ERROR'

! NOTE: 
!   not all compilers behave correctly with print *,mw%comm%wrld%all
!   hence the use of transfer(mw%comm%wrld%all,1) to force a "honest" integer
! mw, ml : OLD style 
! lw, lr : NEW style 
  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, transfer(mw%comm%wrld%all,1)
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, transfer(mw%comm%grid%all,1)
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, transfer(mw%comm%sgrd%all,1)
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, transfer(mw%comm%appl%all,1)
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, transfer(mw%comm%blck%all,1)
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

  print *,'lw%wrld%all%comm  =', lr%wrld%all%comm, transfer(lw%wrld%all%comm,1)
  print *,'lw%grid%all%comm  =', lr%grid%all%comm, transfer(lw%grid%all%comm,1)
  print *,'lw%sgrd%all%comm  =', lr%sgrd%all%comm, transfer(lw%sgrd%all%comm,1)
  if(transfer(lw%sgrd%all%comm,1) .ne. 23456) print *,'ERROR'
  print *,'lw%appl%all%comm  =', lr%appl%all%comm, transfer(lw%appl%all%comm,1)
  print *,'lw%blck%all%comm  =', lr%blck%all%comm, transfer(lw%blck%all%comm,1)
  if(transfer(lw%blck%all%comm,1) .ne. 34567) print *,'ERROR'

  print *,'mw%rank%wrld%all  =', ml%rank%wrld%all, transfer(mw%rank%wrld%all,1)
  if(transfer(mw%rank%wrld%all,1) .ne. 1) print *,'ERROR'
  print *,'mw%size%wrld%all  =', ml%size%wrld%all, transfer(mw%size%wrld%all,1)
  if(transfer(mw%size%wrld%all,1) .ne. 111) print *,'ERROR'
  print *,'mw%rank%sgrd%all  =', ml%rank%sgrd%all, transfer(mw%rank%sgrd%all,1)
  print *,'mw%size%sgrd%all  =', ml%size%sgrd%all, transfer(mw%size%sgrd%all,1)
  if(transfer(mw%size%sgrd%all,1) .ne. 55) print *,'ERROR'
  print *,'mw%rank%grid%all  =', ml%rank%grid%all, transfer(mw%rank%grid%all,1)
  print *,'mw%size%grid%all  =', ml%size%grid%all, transfer(mw%size%grid%all,1)
  if(transfer(mw%size%grid%all,1) .ne. 11) print *,'ERROR'

  print *,'lw%wrld%all%rank  =', lr%wrld%all%rank, transfer(lw%wrld%all%rank,1)
  if(transfer(lw%wrld%all%rank,1) .ne. 1) print *,'ERROR'
  print *,'lw%wrld%all%size  =', lr%wrld%all%size, transfer(lw%wrld%all%size,1)
  if(transfer(lw%wrld%all%size,1) .ne. 111) print *,'ERROR'
  print *,'lw%sgrd%all%rank  =', lr%sgrd%all%rank, transfer(lw%sgrd%all%rank,1)
  print *,'lw%sgrd%all%size  =', lr%sgrd%all%size, transfer(lw%sgrd%all%size,1)
  if(transfer(lw%sgrd%all%size,1) .ne. 55) print *,'ERROR'
  print *,'lw%grid%all%rank  =', lr%grid%all%rank, transfer(lw%grid%all%rank,1)
  print *,'lw%grid%all%size  =', lr%grid%all%size, transfer(lw%grid%all%size,1)
  if(transfer(lw%grid%all%size,1) .ne. 11) print *,'ERROR'

  allocate(mpistatus(dw%MPI_STATUS_SIZE))
!  put address of array into "wrapped pointer"
  buf   = RPN_MPI_Loc(C_LOC(array))    ! way #1, use derived type constructor
  buf   = LoC(array)                   ! way #2, use special macro Loc (CASE SENSITIVE)
  buf%p = C_LOC(array)                 ! way #3, not recommended, will break if name of internal element changes
  if(false) then  ! never executed, syntax test for wrapped types usage when calling the MPI library
    call MPI_Bcast(buf, count, dw%MPI_INTEGER, 0, comm, ier)
    call MPI_Barrier(dw%MPI_COMM_WORLD, ier)
    call MPI_comm_size(comm, count, ier)
    call MPI_send(buf, count, dw%MPI_REAL, dest, tag, comm, ier)
  endif

end subroutine sub3
