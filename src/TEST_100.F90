subroutine rpn_mpi_test_100
! check that the RPN_MPI includes behave correctly
! especially when adding the MPI prototypes using RPN_MPI wrappers
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
  include 'RPN_MPI_mpif.inc'
  print *,'INFO: compilation/load is O.K.'
  print *,'INFO: this is a non MPI test'
  call sub1          ! set values
  call sub2          ! check values from module
  call sub3          ! check values in "user" mode
  print *,'INFO: if no error message was issued, the test is O.K.'
end subroutine rpn_mpi_test_100

subroutine sub1
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
  type(RPN_MPI_Comm) :: dcomm = RPN_MPI_Comm(12345)
  type(RPN_MPI_mpi_definitions) :: dummyd
  call RPN_MPI_init_mpi_layout         ! get wrapped structure pointers mw and dw initialized
  mw%comm%wrld%all = dw%MPI_COMM_WORLD
  mw%rank%wrld%all = 1
  mw%size%wrld%all = 111
  mw%size%sgrd%all = 55
  mw%size%grid%all = 11
  mw%comm%grid%all = dcomm
  mw%comm%sgrd%all = RPN_MPI_Comm(23456)
  mw%comm%blck%all = RPN_MPI_Comm(34567)
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
  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, transfer(mw%comm%wrld%all,1)
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, transfer(mw%comm%grid%all,1)
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, transfer(mw%comm%sgrd%all,1)
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, transfer(mw%comm%appl%all,1)
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, transfer(mw%comm%blck%all,1)
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

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

end subroutine sub2

!============================================================================
!  example of "user" code
!  uses the access functions,
!  the prototypes of which are obtained from include files
!  RPN_MPI.inc       for the RPN_MPI library prototypes
!  RPN_MPI_mpif.inc  for the MPI function prototypes using RPN_MPI wrappers
!============================================================================
subroutine sub3
! perform the same check as sub2 but with copies of structures
! obtained using RPN_MPI_get_mpi_.... subroutines
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
  include 'RPN_MPI_mpif.inc'
  type(RPN_MPI_mpi_definitions_raw) :: dr  ! "raw" MPI symbols
  type(RPN_MPI_mpi_definitions)     :: dw  ! "wrapped" MPI symbols
  type(mpi_layout_internal)         :: ml  ! "raw" RPN_MPI layout information
  type(mpi_layout)                  :: mw  ! "wrapped" RPN_MPI layout information
  integer :: ier, count, dest, tag
  logical :: false = .false.
  integer, dimension(1,2,3) :: array
  type(RPN_MPI_Loc)  :: buf        ! the address of any array
  type(RPN_MPI_Comm) :: comm       ! a communicator
  integer, dimension(:), allocatable :: mpistatus

! get the RPN_MPI "mirror" definitions for subsequent pure MPI calls
  call RPN_MPI_get_mpi_definitions_raw(dr, ier)  ! get a copy of the "raw", everything an integer, structure
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_definitions_raw)'
  call RPN_MPI_get_mpi_definitions(dw, ier)      ! get a copy if the "typed" definitions 
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_definitions)'

! get the RPN_MPI layout (application/supergrid/grid ...) information
  call RPN_MPI_get_mpi_layout_raw(ml, ier)  ! get a copy of the "raw", everything an integer, structure
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout_raw)'
  call RPN_MPI_get_mpi_layout(mw, ier)      ! get a copy if the "typed" definitions 
  if(ier .ne. dr%MPI_SUCCESS) print *,'ERROR (RPN_MPI_get_mpi_layout)'

  print *,'==================== IN USER MODE ===================='
  print *,'MPI_COMM_WORLD    =', transfer(dw%MPI_COMM_WORLD,1), dr%MPI_COMM_WORLD
  if(transfer(dw%MPI_COMM_WORLD,1) .ne. dr%MPI_COMM_WORLD) print *,'ERROR'

  print *,'MPI_COMM_NULL     =', transfer(dw%MPI_COMM_NULL,1), dr%MPI_COMM_NULL
  if(transfer(dw%MPI_COMM_NULL,1) .ne. dr%MPI_COMM_NULL) print *,'ERROR'

  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, transfer(mw%comm%wrld%all,1)
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, transfer(mw%comm%grid%all,1)
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, transfer(mw%comm%sgrd%all,1)
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, transfer(mw%comm%appl%all,1)
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, transfer(mw%comm%blck%all,1)
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

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

  allocate(mpistatus(dw%MPI_STATUS_SIZE))
  buf   = RPN_MPI_Loc(loc(array))    ! recommended way to put address of array into "wrapped pointer"
  buf%a = loc(array)                 ! alternative way, could break if name of internal element changes
  if(false) then  ! never executed, only there to test syntax
    call MPI_Bcast(buf, count, dw%MPI_INTEGER, 0, comm, ier)
    call MPI_Barrier(dw%MPI_COMM_WORLD, ier)
    call MPI_comm_size(comm, count, ier)
    call MPI_send(buf, count, dw%MPI_REAL, dest, tag, comm, ier)
  endif

end subroutine sub3
