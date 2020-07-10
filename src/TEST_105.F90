#if defined(STAND_ALONE)
program test_105
call rpn_mpi_test_105
end
#endif

subroutine rpn_mpi_test_105
!
! simple, standalone test for RPN_MPI_grid_topo
! does not use the full RPN_MPI setup
!
  use ISO_C_BINDING
  implicit none
#include <RPN_MPI.hf>
#define LoC(what) rpn_mpi_loc(loc(what))

  integer :: ier, npe, rank, larg
  type(RPN_MPI_mpi_definitions_raw) :: dr
  integer :: i, j, npex, npey, blkx, blky, stat, x_first, ez
  character(len=128) :: argv(6)
  type(RPN_MPI_Ftopo) :: topo
  integer, dimension(:,:), allocatable :: pe_matrix
  integer, dimension(:), allocatable :: pe_x, pe_y

  call MPI_Init(ier)

  call RPN_MPI_get_mpi_definitions_raw(dr, ier)       ! get "raw" definitions
  call MPI_Comm_size(dr%MPI_COMM_WORLD, npe, ier)     ! get world size
  call MPI_Comm_rank(dr%MPI_COMM_WORLD, rank, ier)    ! get world rank

  topo%grd%comm = RPN_MPI_Comm(dr%MPI_COMM_WORLD)     ! grid communicator (MPI_COMM_WORLD)
  topo%grd%comm%wrapped_value = dr%MPI_COMM_WORLD     ! alternative syntax (not recommended)

  do i=1,6
    call get_command_argument(i,argv(i),larg,stat)    ! get 5 program arguments
    if(stat .ne. 0) goto 777
  enddo
  write(0,*) 'I am PE',rank+1,' of',npe
  read(argv(1),*,err=777) npex                        ! PE topology (npex, npey)
  read(argv(2),*,err=777) npey
  read(argv(3),*,err=777) blkx                        ! block topology (blkx, blky)
  read(argv(4),*,err=777) blky
  read(argv(5),*,err=777) x_first                     ! spread first along x if 1
  read(argv(6),*,err=777) ez                          ! ez test flag (.ne. 0, use ez calls + setup
  if(npex*npey .ne. npe) then
    write(0,*) 'got',npe,' PEs, expected',npex*npey
    goto 777
  endif
  allocate( pe_matrix(0:npex-1,0:npey-1), pe_x(npe), pe_y(npe) )
  pe_x      = 0
  pe_y      = 0
  pe_matrix = -1
  call RPN_MPI_set_grid_topo(npex, npey, blkx, blky, x_first == 1)
  call RPN_MPI_grid_topo(topo, npex, npey, blkx, blky, x_first == 1, ier)
!   call MPI_Gather()

  if(rank == 0) then
    do j = npey-1, 0, -1
      write(0,2) j, pe_matrix(:,j)
    enddo
    write(0,3) (i, i=0,npex-1)
  endif

1 call MPI_Finalize(ier)
  return

777 continue
  write(0,*) 'ERROR in arguments'
  goto 1
2 format(30I10)
3 format(10X,29I10)
end subroutine rpn_mpi_test_105
