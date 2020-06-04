!   include 'RPN_MPI_mpi_symbols.inc'
!   include 'RPN_MPI_mpi_layout.inc'
!   interface
!     subroutine RPN_MPI_get_mpi_definitions(what, ierr)
!       import :: RPN_MPI_mpi_definitions
!       type(RPN_MPI_mpi_definitions), intent(INOUT) :: what
!       integer, intent(OUT) :: ierr
!     end subroutine RPN_MPI_get_mpi_definitions
!     subroutine RPN_MPI_get_mpi_layout(what, ierr)
!       import :: mpi_layout, C_INT
!       type(mpi_layout), intent(INOUT) :: what
!       integer(C_INT), intent(OUT) :: ierr
!     end subroutine RPN_MPI_get_mpi_layout
!     subroutine RPN_MPI_get_mpi_layout_raw(what, ierr)
!       import :: mpi_layout, C_INT
!       type(mpi_layout_internal), intent(INOUT) :: what
!       integer(C_INT), intent(OUT) :: ierr
!     end subroutine RPN_MPI_get_mpi_layout
!   end interface
