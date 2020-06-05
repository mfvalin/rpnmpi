!
! interfaces to some MPI routines, using RPN_MPI wrappers
!
interface
  subroutine MPI_Send(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
    import :: RPN_MPI_Comm, RPN_MPI_Datatype, RPN_MPI_Loc
    implicit none
    type(RPN_MPI_Loc), intent(IN), value :: BUF
    integer, intent(IN) :: COUNT, DEST, TAG
    type(RPN_MPI_Datatype), intent(IN) :: DATATYPE
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Send
  subroutine MPI_Recv(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
    import :: RPN_MPI_Comm, RPN_MPI_Datatype, RPN_MPI_Loc
    implicit none
    type(RPN_MPI_Loc), intent(IN), value :: BUF
    integer, intent(IN) :: COUNT, SOURCE, TAG
    type(RPN_MPI_Datatype), intent(IN) :: DATATYPE
    integer, intent(OUT), dimension(*) :: STATUS
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Recv
  subroutine MPI_Sendecv(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
    import :: RPN_MPI_Comm, RPN_MPI_Datatype, RPN_MPI_Loc
    implicit none
    type(RPN_MPI_Loc), intent(IN), value :: SENDBUF, RECVBUF
    integer, intent(IN) :: SENDCOUNT, RECVCOUNT, SOURCE, DEST, SENDTAG, RECVTAG
    type(RPN_MPI_Datatype), intent(IN) :: SENDTYPE, RECVTYPE
    integer, intent(OUT), dimension(*) :: STATUS
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Sendecv
  subroutine MPI_Comm_size(COMM, SIZE, IERROR)
    import :: RPN_MPI_Comm
    implicit none
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: SIZE
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_size
  subroutine MPI_Comm_rank(COMM, RANK, IERROR)
    import :: RPN_MPI_Comm
    implicit none
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: RANK
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_rank
  subroutine MPI_Barrier(COMM, IERROR)
    import :: RPN_MPI_Comm
    implicit none
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Barrier
  subroutine MPI_Comm_split(COMM, COLOR, KEY, NEWCOMM, IERROR)
    import :: RPN_MPI_Comm
    implicit none
    type(RPN_MPI_Comm), intent(IN) :: COMM
    integer, intent(IN) :: COLOR, KEY
    type(RPN_MPI_Comm), intent(OUT) :: NEWCOMM
    integer, intent(OUT) :: IERROR
  end subroutine MPI_Comm_split
  function MPI_Wtime() result(t)
    real(kind=8) :: t
  end function MPI_Wtime
end interface
