#if defined(TEST)
program my_test
call mpi_init(ierr)
call rpn_mpi_test_101
call mpi_finalize(ierr)
stop
end
#endif
subroutine rpn_mpi_test_101
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_MPI.inc'
  integer :: ierr, narg, arglen, rank, pop
  external :: Userinit
  integer :: Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,Io
  character(len=128) :: AppID, string

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, pop, ierr)
  do narg = 1, 6
    if(narg == 1) then
      call get_command_argument(1 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Pex
!       print *,'pex =',pex
    endif
    if(narg == 2) then
      call get_command_argument(2 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Pey
!       print *,'pey =',pey
    endif
    if(narg == 3) then
      call get_command_argument(3 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Multigrids
!       print *,'Multigrids =',Multigrids
    endif
    if(narg == 4) then
      call get_command_argument(4 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Grids
!       print *,'Grids =',Grids
    endif
    if(narg == 5) then
      call get_command_argument(5 , AppID, arglen, ierr)
      if(ierr .ne. 0) goto 999
!       print *,'AppID = "',trim(AppID)//'"'
    endif
    if(narg == 6) then
      call get_command_argument(6 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Io
!       print *,'Io =',Io
    endif
  enddo

  if(Io >= 0) then
    if(rank >= pop/2) AppID = "<N02>"
    ierr = RPN_MPI_init(Userinit,Pelocal,Petotal,Pex,Pey,Multigrids,Grids,AppID,Io)
  else
    ierr = RPN_COMM_init_multi_level(Userinit,Pelocal,Petotal,Pex,Pey,Multigrids,Grids)
  endif

777 continue
  call RPN_MPI_finalize(ierr)
  stop
999 continue
  print *,"ERROR in arguments",ierr
  go to 777
  return
end subroutine rpn_mpi_test_101

subroutine Userinit
  return
end subroutine Userinit
