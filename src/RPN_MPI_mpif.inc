!
! interfaces to some MPI routines, using RPN_MPI wrappers
!
interface
  subroutine MPI_Send(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
    import :: RPN_MPI_Comm, RPN_MPI_Datatype, RPN_MPI_Loc
    type(RPN_MPI_Loc), intent(IN), value :: BUF
    integer, intent(IN) :: COUNT, DEST, TAG
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Send
  subroutine MPI_Comm_size(COMM, SIZE, IERROR)
    import :: RPN_MPI_Comm
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: SIZE
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_size
  subroutine MPI_Comm_rank(COMM, RANK, IERROR)
    import :: RPN_MPI_Comm
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: RANK
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_rank
  subroutine MPI_Barrier(COMM, IERROR)
    import :: RPN_MPI_Comm
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Barrier
  subroutine MPI_Comm_split(COMM, COLOR, KEY, NEWCOMM, IERROR)
    import :: RPN_MPI_Comm
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(IN) :: COLOR, KEY
    type(RPN_MPI_Comm), intent(OUT) :: NEWCOMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_split
  function MPI_Wtime() result(t)
    real(kind=8) :: t
  end function MPI_Wtime
end interface
