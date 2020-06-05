subroutine rpn_mpi_test_100
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
  print *,'INFO: compilation/load test O.K.'
  call sub1          ! set values
  call sub2          ! check values from module
  call sub3          ! check values in "user" mode
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
subroutine sub2
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

  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, mw%comm%wrld%all
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, mw%comm%wrld%all
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, mw%comm%wrld%all
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, mw%comm%wrld%all
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, mw%comm%wrld%all
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

  print *,'mw%rank%wrld%all  =', ml%rank%wrld%all, mw%rank%wrld%all
  if(transfer(mw%rank%wrld%all,1) .ne. 1) print *,'ERROR'
  print *,'mw%size%wrld%all  =', ml%size%wrld%all, mw%size%wrld%all
  if(transfer(mw%size%wrld%all,1) .ne. 111) print *,'ERROR'
  print *,'mw%rank%sgrd%all  =', ml%rank%sgrd%all, mw%rank%sgrd%all
  print *,'mw%size%sgrd%all  =', ml%size%sgrd%all, mw%size%sgrd%all
  if(transfer(mw%size%sgrd%all,1) .ne. 55) print *,'ERROR'
  print *,'mw%rank%grid%all  =', ml%rank%grid%all, mw%rank%grid%all
  print *,'mw%size%grid%all  =', ml%size%grid%all, mw%size%grid%all
  if(transfer(mw%size%grid%all,1) .ne. 11) print *,'ERROR'

end subroutine sub2

subroutine sub3
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
  type(RPN_MPI_mpi_definitions_raw) :: dr
  type(RPN_MPI_mpi_definitions) :: dw
  type(mpi_layout_internal) :: ml
  type(mpi_layout) :: mw
  integer :: ierr

  call RPN_MPI_get_mpi_definitions_raw(dr, ierr)
  call RPN_MPI_get_mpi_definitions(dw, ierr)
  call RPN_MPI_get_mpi_layout_raw(ml, ierr)
  call RPN_MPI_get_mpi_layout(mw, ierr)

  print *,'==================== IN USER MODE ===================='
  print *,'MPI_COMM_WORLD    =', dw%MPI_COMM_WORLD, dr%MPI_COMM_WORLD
  if(transfer(dw%MPI_COMM_WORLD,1) .ne. dr%MPI_COMM_WORLD) print *,'ERROR'

  print *,'MPI_COMM_NULL     =', dw%MPI_COMM_NULL, dr%MPI_COMM_NULL
  if(transfer(dw%MPI_COMM_NULL,1) .ne. dr%MPI_COMM_NULL) print *,'ERROR'

  print *,'mw%comm%wrld%all  =', ml%comm%wrld%all, mw%comm%wrld%all
  print *,'mw%comm%grid%all  =', ml%comm%grid%all, mw%comm%wrld%all
  print *,'mw%comm%sgrd%all  =', ml%comm%sgrd%all, mw%comm%wrld%all
  if(transfer(mw%comm%sgrd%all,1) .ne. 23456) print *,'ERROR'
  print *,'mw%comm%appl%all  =', ml%comm%appl%all, mw%comm%wrld%all
  print *,'mw%comm%blck%all  =', ml%comm%blck%all, mw%comm%wrld%all
  if(transfer(mw%comm%blck%all,1) .ne. 34567) print *,'ERROR'

  print *,'mw%rank%wrld%all  =', ml%rank%wrld%all, mw%rank%wrld%all
  if(transfer(mw%rank%wrld%all,1) .ne. 1) print *,'ERROR'
  print *,'mw%size%wrld%all  =', ml%size%wrld%all, mw%size%wrld%all
  if(transfer(mw%size%wrld%all,1) .ne. 111) print *,'ERROR'
  print *,'mw%rank%sgrd%all  =', ml%rank%sgrd%all, mw%rank%sgrd%all
  print *,'mw%size%sgrd%all  =', ml%size%sgrd%all, mw%size%sgrd%all
  if(transfer(mw%size%sgrd%all,1) .ne. 55) print *,'ERROR'
  print *,'mw%rank%grid%all  =', ml%rank%grid%all, mw%rank%grid%all
  print *,'mw%size%grid%all  =', ml%size%grid%all, mw%size%grid%all
  if(transfer(mw%size%grid%all,1) .ne. 11) print *,'ERROR'

end subroutine sub3
