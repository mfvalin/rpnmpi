subroutine rpn_mpi_test_104
!
! simple, standalone test for RPN_MPI_ez_halo_parms
! does not use the full RPN_MPI setup
!
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
!   include 'RPN_MPI_mpif.inc'

  integer :: ier, npe, rank
  type(RPN_MPI_mpi_definitions_raw) :: dr
  integer :: i, j, k, larg1, larg2, larg3, stat1, stat2, stat3, lni, lnj, lnk, gni, gnj, gnk
  character(len=128) :: argv1, argv2, argv3, mode
  integer :: lnimax, lnjmax, lnkmax
  integer, dimension(:,:,:), allocatable   :: z
  integer, dimension(:,:,:,:), allocatable :: zt

  call MPI_Init(ier)

  call RPN_MPI_get_mpi_definitions_raw(dr, ier)  ! get "raw" definitions
  call MPI_Comm_size(dr%MPI_COMM_WORLD, npe, ier)
  call MPI_Comm_rank(dr%MPI_COMM_WORLD, rank, ier)

  call get_command_argument(1,argv1,larg1,stat1)
  call get_command_argument(2,argv2,larg2,stat2)
  call get_command_argument(3,argv3,larg3,stat3)
  if(stat1 .ne. 0 .or. stat2 .ne. 0 .or. stat3 .ne. 0) goto 777
  write(0,*) 'I am PE',rank+1,' of',npe
  read(argv1,*,err=777) gni
  read(argv2,*,err=777) gnj
  read(argv3,*,err=777) gnk
  lnimax = (gni + npe -1) / npe
  lnjmax = (gnj + npe -1) / npe
  lnkmax = (gnk + npe -1) / npe
  lni = min(lnimax, gni - (lnimax * rank))
  lni = max(0,lni)
  lnj = min(lnjmax, gnj - (lnjmax * rank))
  lnj = max(0,lnj)
  lnk = min(lnkmax, gnk - (lnkmax * rank))
  lnk = max(0,lnk)
  lnjmax = (gnj + npe -1) / npe
  lnkmax = (gnk + npe -1) / npe
  write(0,*) gni, gnj, gnk, lni, lnj, lnk
  allocate(z(lni,lnj,gnk))
  allocate(zt(lni,lnj,max(1,lnk),npe))
  z = 0
  zt = 0
  write(0,*) sum(z), sum(zt)
1 call MPI_Finalize(ier)
  return

777 continue
  write(0,*) 'ERROR in arguments'
  goto 1
end subroutine rpn_mpi_test_104
