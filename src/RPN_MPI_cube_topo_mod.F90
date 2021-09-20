! RPN_MPI - Library of useful routines for C and FORTRAN programming
! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
!                          Environnement Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

module rpn_mpi_types_mod
  use ISO_C_BINDING
  use mpi_f08
  implicit none

  type :: comm_size_rank
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer(C_INT) :: size = -1
    integer(C_INT) :: rank = -1
  end type
  type(comm_size_rank), parameter :: MPI_COMM_SIZE_RANK_NULL = comm_size_rank(MPI_COMM_NULL, -1, -1)

  type :: neighbor
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer(C_INT) :: rank = -1
    logical :: flip = .false.
  end type
  type(neighbor), parameter :: NEIGHBOR_NULL = neighbor(MPI_COMM_NULL, -1, .false.)

  type :: cube_topology
    type(comm_size_rank) :: cube_comm = MPI_COMM_SIZE_RANK_NULL
    integer :: side
    type(comm_size_rank) :: side_comm = MPI_COMM_SIZE_RANK_NULL
    integer :: row
    type(comm_size_rank) :: row_comm = MPI_COMM_SIZE_RANK_NULL
    integer :: column
    type(comm_size_rank) :: col_comm = MPI_COMM_SIZE_RANK_NULL
    type(neighbor) :: north = NEIGHBOR_NULL
    type(neighbor) :: south = NEIGHBOR_NULL
    type(neighbor) :: east = NEIGHBOR_NULL
    type(neighbor) :: west = NEIGHBOR_NULL
  end type

  integer, parameter :: NORTH_EDGE = 1
  integer, parameter :: SOUTH_EDGE = 2
  integer, parameter :: EAST_EDGE  = 3
  integer, parameter :: WEST_EDGE  = 4

  contains

subroutine rpn_mpi_cube_exchange(cube, edges_s, edges_r, nv, lnij, nk)
  implicit none
  type(cube_topology), intent(IN) :: cube    ! cube topology
  integer, intent(IN) :: nv, lnij, nk
  integer, dimension(nv, lnij, nk, 4), intent(IN)  :: edges_s
  integer, dimension(nv, lnij, nk, 4), intent(OUT) :: edges_r

  integer, dimension(nv, lnij, nk, 4) :: edges_f   ! local space if an edge must be flipped before sending
  TYPE(MPI_Request), dimension(4) :: request_r, request_s
  TYPE(MPI_Status), dimension(4)  :: status_r, status_s
  integer :: count

  ! post asynchronous receives first
  count = nv * lnij * nk
  !    MPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST)
  call MPI_Irecv(edges_r(1, 1, 1, NORTH_EDGE), count, MPI_INTEGER, cube%north%rank, 0, cube%north%comm, request_r(NORTH_EDGE) )
  call MPI_Irecv(edges_r(1, 1, 1, SOUTH_EDGE), count, MPI_INTEGER, cube%south%rank, 0, cube%south%comm, request_r(SOUTH_EDGE) )
  call MPI_Irecv(edges_r(1, 1, 1, EAST_EDGE ), count, MPI_INTEGER, cube%east%rank,  0, cube%east%comm,  request_r(EAST_EDGE ) )
  call MPI_Irecv(edges_r(1, 1, 1, WEST_EDGE ), count, MPI_INTEGER, cube%west%rank,  0, cube%west%comm,  request_r(WEST_EDGE ) )

  ! post asynchronous sends
  !    MPI_ISEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST)
  if(cube%north%flip) then
    call flip_edge( edges_s(1, 1, 1, NORTH_EDGE), edges_f(1, 1, 1, NORTH_EDGE), nv, lnij, nk)
    call MPI_Isend(edges_f(1, 1, 1, NORTH_EDGE), count, MPI_INTEGER, cube%north%rank, 0, cube%north%comm, request_s(NORTH_EDGE) )
  else
    call MPI_Isend(edges_s(1, 1, 1, NORTH_EDGE), count, MPI_INTEGER, cube%north%rank, 0, cube%north%comm, request_s(NORTH_EDGE) )
  endif
  if(cube%south%flip) then
    call flip_edge( edges_s(1, 1, 1, SOUTH_EDGE), edges_f(1, 1, 1, SOUTH_EDGE), nv, lnij, nk)
    call MPI_Isend(edges_f(1, 1, 1, SOUTH_EDGE), count, MPI_INTEGER, cube%south%rank, 0, cube%south%comm, request_s(SOUTH_EDGE) )
  else
    call MPI_Isend(edges_s(1, 1, 1, SOUTH_EDGE), count, MPI_INTEGER, cube%south%rank, 0, cube%south%comm, request_s(SOUTH_EDGE) )
  endif
  if(cube%east%flip) then
    call flip_edge( edges_s(1, 1, 1, EAST_EDGE), edges_f(1, 1, 1, EAST_EDGE), nv, lnij, nk)
    call MPI_Isend(edges_f(1, 1, 1, EAST_EDGE ), count, MPI_INTEGER, cube%east%rank,  0, cube%east%comm,  request_s(EAST_EDGE ) )
  else
    call MPI_Isend(edges_s(1, 1, 1, EAST_EDGE ), count, MPI_INTEGER, cube%east%rank,  0, cube%east%comm,  request_s(EAST_EDGE ) )
  endif
  if(cube%west%flip) then
    call flip_edge( edges_s(1, 1, 1, WEST_EDGE), edges_f(1, 1, 1, WEST_EDGE), nv, lnij, nk)
    call MPI_Isend(edges_f(1, 1, 1, WEST_EDGE ), count, MPI_INTEGER, cube%west%rank,  0, cube%west%comm,  request_s(WEST_EDGE ) )
  else
    call MPI_Isend(edges_s(1, 1, 1, WEST_EDGE ), count, MPI_INTEGER, cube%west%rank,  0, cube%west%comm,  request_s(WEST_EDGE ) )
  endif

  ! MPI_Waitall(count, array_of_requests, array_of_statuses)
  call MPI_Waitall(4, request_r, status_r)
  call MPI_Waitall(4, request_s, status_s)
  count = compare_edges(edges_s, edges_r, nv, lnij, nk)
  print *,'Edges exchanged, errors =', count

  contains
  function compare_edges(edge_1, edge_2, nv, lnij, nk) result(diff)
    implicit none
    integer, intent(IN) :: nv, lnij, nk
    integer, dimension(nv, lnij, nk, 4), intent(IN)   :: edge_1, edge_2
    integer :: diff
    integer :: i, j, k, l

    diff = 0
    do l = 1, 4
    do k = 1, nk
    do j = 1, lnij
    do i = 1, nv
      if(edge_1(i,j,k,l) .ne. edge_2(i,j,k,l)) diff = diff + 1
    enddo
    enddo
    enddo
    enddo
    print *, diff,' differences out of',nv*lnij*nk*4
  end function compare_edges
  subroutine flip_edge(edge_i, edge_o, nv, lnij, nk)
    implicit none
    integer, intent(IN) :: nv, lnij, nk
    integer, dimension(nv, lnij, nk), intent(IN)   :: edge_i
    integer, dimension(nv, lnij, nk), intent(OUT)  :: edge_o

    edge_o(:, lnij:1:-1, :) = edge_i(:, 1:lnij:1, :)
  end subroutine flip_edge
end subroutine rpn_mpi_cube_exchange

subroutine rpn_mpi_set_cube_topology(csr, cube, npe)
  use ISO_C_BINDING
  implicit none
  type(comm_size_rank), intent(IN) :: csr     ! communicator large enough for the npe x npe x 6 cube
  type(cube_topology), intent(OUT) :: cube    ! cube topology
  integer, intent(IN) :: npe                  ! each side of the cube will use npe x npe PEs

  type(MPI_Comm) cube_comm
  integer :: ncube, nside
  logical :: north, south, east, west
  integer :: my_row, my_col, flip_row, flip_col, last, incube

  nside = npe * npe
  ncube = 6 * nside

  if(csr%size < ncube) return    ! not enough processes
  incube = 1
  if(csr%rank >= ncube) incube = 0   ! not part of the cube
  call MPI_Comm_split(csr%comm, incube, csr%rank, cube%cube_comm%comm)     ! communicator for the cube
  if(incube == 0) then
    cube%cube_comm%comm = MPI_COMM_NULL   ! not in cube
    return
  endif
  call MPI_Comm_size(cube%cube_comm%comm, cube%cube_comm%size)
  call MPI_Comm_rank(cube%cube_comm%comm, cube%cube_comm%rank)

  cube_comm = cube%cube_comm%comm
  cube%side = csr%rank / (npe*npe)                       ! side number for this PE
  call MPI_Comm_split(cube_comm, cube%side, csr%rank, cube%side_comm%comm)  ! split cube into sides
  call MPI_Comm_size(cube%side_comm%comm, cube%side_comm%size)
  call MPI_Comm_rank(cube%side_comm%comm, cube%side_comm%rank)     ! rank in side

  cube%column = mod(cube%side_comm%rank, npe)    ! column number in side
  my_col = cube%column
  flip_col = npe -1 - my_col
  call MPI_Comm_split(cube%side_comm%comm, my_col, csr%rank, cube%col_comm%comm) ! split side into columns
  call MPI_Comm_size(cube%col_comm%comm, cube%col_comm%size)
  call MPI_Comm_rank(cube%col_comm%comm, cube%col_comm%rank)

  cube%row = cube%side_comm%rank / npe           ! row number in side
  my_row = cube%row
  flip_row = npe -1 - my_row
  call MPI_Comm_split(cube%side_comm%comm, my_row, csr%rank, cube%row_comm%comm) ! split side into rows
  call MPI_Comm_size(cube%row_comm%comm, cube%row_comm%size)
  call MPI_Comm_rank(cube%row_comm%comm, cube%row_comm%rank)

  last = npe - 1

  north = cube%row == last
  south = cube%row == 0
  east = cube%column == last
  west = cube%column == 0

  ! not on an edge, same for all sides of the cube
  if(.not. west)  cube%west  = neighbor(cube%col_comm%comm, my_col - 1, .false.)        ! column to the left
  if(.not. east)  cube%east  = neighbor(cube%col_comm%comm, my_col + 1, .false.)        ! column to the right
  if(.not. north) cube%north = neighbor(cube%row_comm%comm, my_row + 1, .false.)        ! row above
  if(.not. south) cube%south = neighbor(cube%row_comm%comm, my_row - 1, .false.)        ! row below
  if( (.not. west) .and. (.not. east) .and. (.not. north) .and. (.not. south)) return   ! job done

  ! take care of edges for all sides of the cube
  select case(cube%side)
    case(0)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(3, last,  my_row), .false.)
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(1, 0    , my_row), .false.)
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(4, my_col, 0    ), .false.)
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(5, my_col, last ), .false.)
      endif

    case(1)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(0, last, my_row  ), .false.)
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(2, 0,    my_row  ), .false.)
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(4, last, my_col  ), .false.)
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(5, last, flip_col), .true. )
      endif

    case(2)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(1, last,     my_row), .false.)
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(3, 0,        my_row), .false.)
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(4, flip_col, last  ), .true. )
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(5, flip_col, 0     ), .true. )
      endif

    case(3)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(2, last, my_row  ), .false.)
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(0, 0,    my_row  ), .false.)
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(4, 0,    flip_col), .true. )
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(5, 0,    my_col  ), .false.)
      endif

    case(4)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(3, flip_row, last), .true. )
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(1, my_row,   last), .false.)
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(2, flip_col, last), .true. )
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(0, my_col,   last), .false.)
      endif

    case(5)
      if(west) then
        cube%west  = neighbor(cube_comm, cube_rank(3, my_row,   0), .false.)
      endif
      if(east) then
        cube%east  = neighbor(cube_comm, cube_rank(1, flip_row, 0), .true. )
      endif
      if(north) then
        cube%north = neighbor(cube_comm, cube_rank(0, my_col,   0), .false.)
      endif
      if(south) then
        cube%south = neighbor(cube_comm, cube_rank(2, flip_col, 0), .true.)
      endif

    case DEFAULT

  end select

  contains
  function cube_rank(side, column, row) result(rank)
    implicit none
    integer, intent(IN) :: side, column, row
    integer :: rank
    rank = (side * npe * npe) + row * npe + column
  end function cube_rank
end subroutine rpn_mpi_set_cube_topology

  function edge_compare(edge0, edge1, nv, lnij, flip) result(same)
    use ISO_C_BINDING
    implicit none
    integer, intent(IN) :: lnij
    integer, intent(IN) :: nv
    integer, dimension(nv, lnij), intent(IN) :: edge0
    integer, dimension(nv, lnij), intent(IN) :: edge1
    logical, intent(IN) :: flip
    logical :: same

    if(flip) then
      same = all(edge0(1:nv, 1:lnij) .eq. edge1(1:nv, lnij:1:-1))
    else
      same = all(edge0(1:nv, 1:lnij) .eq. edge1(1:nv, 1:lnij))
    endif
    if(.not. same) then
      print 1,edge0(1:nv, 1:lnij)
      print 1,edge1(1:nv, 1:lnij)
    endif
1 format(10(3I3,2X))
  end function edge_compare

end module

#if defined(TEST_EDGES)
program test
  call test_edges
end program
#endif

#if defined(TEST_EXCHANGE)
program test
  use ISO_C_BINDING
  use rpn_mpi_types_mod
  implicit none
  integer :: npe
  type(cube_topology) :: cube
  type(comm_size_rank) :: csr
  integer, parameter :: LNIJ = 9
  integer, parameter :: XYZ = 3
  integer, parameter :: NK = 3
  integer,  dimension(XYZ, LNIJ, NK, 4) :: edges_s  ! edges to send
  integer,  dimension(XYZ, LNIJ, NK, 4) :: edges_r  ! edges to receive

  call mpi_init()
  csr%comm = MPI_COMM_WORLD
  call mpi_comm_size(csr%comm, csr%size)
  call mpi_comm_rank(csr%comm, csr%rank)
  npe = nint(sqrt(csr%size/6*1.0))               ! each side will use npe x npe PEs

  call rpn_mpi_set_cube_topology(csr, cube, npe) ! setup of cube topology data
  print *,'PE',csr%rank+1,' of',csr%size, ', side/row/col =',cube%side,cube%column,cube%row

!      test_edges_create(side,      nv,  lnij, nk, pei,         pej,      npe, edges)
  call test_edges_create(cube%side, XYZ, LNIJ, NK, cube%column, cube%row, npe, edges_s) ! create edges to send
  call rpn_mpi_cube_exchange(cube, edges_s, edges_r, XYZ, LNIJ, NK)

  call mpi_finalize
end program
#endif

subroutine test_edges
  use ISO_C_BINDING
  use rpn_mpi_types_mod
  implicit none
  integer, parameter :: lnij = 5
  integer, parameter :: nv = 3
  integer, parameter :: nk = 1
  integer, parameter :: np = 3
!   integer, parameter :: nm = np-1

  character(len=2), dimension(4) :: edgename
  integer :: side, i, j, ip, pei, pej
  integer,  dimension(nv, lnij, nk, 4, 0:5, 0:np-1, 0:np-1) :: edges

  edges = 888
  edgename = [ 'N ', 'S ', 'E ', 'W']
  do j = 0, 5
    side = j
    print 2,'========== (',j,') =========='
    do pej = np-1, 0, -1
      do pei = 0, np-1
        call test_edges_create(j, nv, lnij, nk, pei, pej, np, edges(1, 1, 1, 1, j, pei, pej))
      enddo
      do i = 1, 4
!         print 1, edgename(i),edges(1:3,1:lnij,1,i,j,:,pej)
      enddo
!     print *,''
    enddo
!     do i = 1, 4
!       print 1, edgename(i),edges(1:3,1:lnij,1,i,j,0,0)
!     enddo
!     print *,''
  enddo
  do ip = 0, np-1
  if(edge_compare(edges(1,1,1,EAST_EDGE,0,np-1,ip) , edges(1,1,1,WEST_EDGE,1,0,ip), nv, lnij, .false.))  print *,'EAST 0  == WEST  1, p =',ip
  if(edge_compare(edges(1,1,1,EAST_EDGE,1,np-1,ip) , edges(1,1,1,WEST_EDGE,2,0,ip), nv, lnij, .false.))  print *,'EAST 1  == WEST  2, p =',ip
  if(edge_compare(edges(1,1,1,EAST_EDGE,2,np-1,ip) , edges(1,1,1,WEST_EDGE,3,0,ip), nv, lnij, .false.))  print *,'EAST 2  == WEST  3, p =',ip
  if(edge_compare(edges(1,1,1,EAST_EDGE,3,np-1,ip) , edges(1,1,1,WEST_EDGE,0,0,ip), nv, lnij, .false.))  print *,'EAST 3  == WEST  0, p =',ip
! 
  if(edge_compare(edges(1,1,1,NORTH_EDGE,0,ip,np-1) , edges(1,1,1,SOUTH_EDGE,4,ip,     0),       nv, lnij, .false.)) print *,'NORTH 0 == SOUTH 4, p =',ip
  if(edge_compare(edges(1,1,1,NORTH_EDGE,1,ip,np-1) , edges(1,1,1,EAST_EDGE, 4,np-1,   ip),      nv, lnij, .false.)) print *,'NORTH 1 == EAST  4, p =',ip
  if(edge_compare(edges(1,1,1,NORTH_EDGE,2,ip,np-1) , edges(1,1,1,NORTH_EDGE,4,np-1-ip,np-1),    nv, lnij, .true. )) print *,'NORTH 2 <> NORTH 4, p =',ip
  if(edge_compare(edges(1,1,1,NORTH_EDGE,3,ip,np-1) , edges(1,1,1,WEST_EDGE, 4,0,      np-1-ip), nv, lnij, .true. )) print *,'NORTH 3 <> WEST  4, p =',ip
! 
  if(edge_compare(edges(1,1,1,SOUTH_EDGE,0,ip,0) , edges(1,1,1,NORTH_EDGE,5,ip,     np-1),    nv, lnij, .false.)) print *,'SOUTH 0 == NORTH 5, p =',ip
  if(edge_compare(edges(1,1,1,SOUTH_EDGE,1,ip,0) , edges(1,1,1,EAST_EDGE, 5,np-1,   np-1-ip), nv, lnij, .true. )) print *,'SOUTH 1 <> EAST  5, p =',ip
  if(edge_compare(edges(1,1,1,SOUTH_EDGE,2,ip,0) , edges(1,1,1,SOUTH_EDGE,5,np-1-ip,0),       nv, lnij, .true. )) print *,'SOUTH 2 <> SOUTH 5, p =',ip
  if(edge_compare(edges(1,1,1,SOUTH_EDGE,3,ip,0) , edges(1,1,1,WEST_EDGE, 5,0,     ip),       nv, lnij, .false.)) print *,'SOUTH 3 == WEST  5, p =',ip
  enddo

1 format(A3,'(X,Y,Z) : ',20(' ',3(I3,1X),' ',2X))
2 format(A,I1,A)
end subroutine test_edges

! create N/S/E/W edges for PE(pei,pej) from side side
! there are npe PEs along each axis of the cube
! nv is expected to be 3 for this test
subroutine test_edges_create(side, nv, lnij, nk, pei, pej, npe, edges)
  use ISO_C_BINDING
  use rpn_mpi_types_mod
  implicit none
  integer, intent(IN) :: side
  integer, intent(IN) :: lnij
  integer, intent(IN) :: nv
  integer, intent(IN) :: nk
  integer, intent(IN) :: pei
  integer, intent(IN) :: pej
  integer, intent(IN) :: npe
  integer, dimension(nv, lnij, nk, 4), intent(OUT) :: edges

  integer :: i
  integer, parameter :: X = 1
  integer, parameter :: Y = 2
  integer, parameter :: Z = 3
  integer :: gnij, i20, j20

  gnij = lnij * npe       ! global number of points along i and j
  i20  = lnij * pei * 2   ! offset for tile along i
  j20  = lnij * pej * 2   ! offset for tile along j
! print *,'edge =', side, ',  lnij =',lnij
! print *,'pei, pej =',pei,pej
edges = 999
  select case(side)  ! fill interface edges accorfing to side
    case(0)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) =  gnij                   ! x = gnij
        edges(X, i, 1, SOUTH_EDGE) =  gnij 
        edges(X, i, 1, EAST_EDGE)  =  gnij 
        edges(X, i, 1, WEST_EDGE)  =  gnij 
        edges(Y, i, 1, NORTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, SOUTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, EAST_EDGE)  = -gnij + i20 + 2*lnij    ! y max
        edges(Y, i, 1, WEST_EDGE)  = -gnij + i20             ! y min
        edges(Z, i, 1, NORTH_EDGE) = -gnij + j20 + 2*lnij    ! z max
        edges(Z, i, 1, SOUTH_EDGE) = -gnij + j20             ! z min
        edges(Z, i, 1, EAST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
        edges(Z, i, 1, WEST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
      enddo
  
    case(1)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) =  gnij - i20 - 2*i + 1   ! -x axis
        edges(X, i, 1, SOUTH_EDGE) =  gnij - i20 - 2*i + 1   ! -x axis
        edges(X, i, 1, EAST_EDGE)  =  gnij - i20 - 2*lnij    ! x min
        edges(X, i, 1, WEST_EDGE)  =  gnij - i20             ! x max
        edges(Y, i, 1, NORTH_EDGE) =  gnij                   ! y = gnij
        edges(Y, i, 1, SOUTH_EDGE) =  gnij 
        edges(Y, i, 1, EAST_EDGE)  =  gnij 
        edges(Y, i, 1, WEST_EDGE)  =  gnij 
        edges(Z, i, 1, NORTH_EDGE) = -gnij + j20 + 2*lnij    ! z max
        edges(Z, i, 1, SOUTH_EDGE) = -gnij + j20             ! z min
        edges(Z, i, 1, EAST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
        edges(Z, i, 1, WEST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
      enddo
  
    case(2)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) = -gnij                   ! x = -lnij
        edges(X, i, 1, SOUTH_EDGE) = -gnij 
        edges(X, i, 1, EAST_EDGE)  = -gnij 
        edges(X, i, 1, WEST_EDGE)  = -gnij 
        edges(Y, i, 1, NORTH_EDGE) =  gnij - i20 - 2*i + 1   ! -y axis
        edges(Y, i, 1, SOUTH_EDGE) =  gnij - i20 - 2*i + 1   ! -y axis
        edges(Y, i, 1, EAST_EDGE)  =  gnij - i20 - 2*lnij    ! y min
        edges(Y, i, 1, WEST_EDGE)  =  gnij - i20             ! y max
        edges(Z, i, 1, NORTH_EDGE) = -gnij + j20 + 2*lnij    ! z max
        edges(Z, i, 1, SOUTH_EDGE) = -gnij + j20             ! z min
        edges(Z, i, 1, EAST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
        edges(Z, i, 1, WEST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
      enddo
  
    case(3)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) = -gnij + i20 + 2*i - 1   ! +x axis
        edges(X, i, 1, SOUTH_EDGE) = -gnij + i20 + 2*i - 1   ! +x axis
        edges(X, i, 1, EAST_EDGE)  = -gnij + i20 + 2*lnij    ! x max
        edges(X, i, 1, WEST_EDGE)  = -gnij + i20             ! x min
        edges(Y, i, 1, NORTH_EDGE) = -gnij                   ! y = -lnij
        edges(Y, i, 1, SOUTH_EDGE) = -gnij 
        edges(Y, i, 1, EAST_EDGE)  = -gnij 
        edges(Y, i, 1, WEST_EDGE)  = -gnij 
        edges(Z, i, 1, NORTH_EDGE) = -gnij + j20 + 2*lnij    ! z max
        edges(Z, i, 1, SOUTH_EDGE) = -gnij + j20             ! z min
        edges(Z, i, 1, EAST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
        edges(Z, i, 1, WEST_EDGE)  = -gnij + j20 + 2*i - 1   ! +z axis
      enddo
  
    case(4)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) =  gnij - j20 - 2*lnij    ! x min
        edges(X, i, 1, SOUTH_EDGE) =  gnij - j20             ! x max
        edges(X, i, 1, EAST_EDGE)  =  gnij - j20 -2*i + 1    ! -x axis
        edges(X, i, 1, WEST_EDGE)  =  gnij - j20 -2*i + 1    ! -x axis
        edges(Y, i, 1, NORTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, SOUTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, EAST_EDGE)  = -gnij + i20 + 2*lnij    ! y max
        edges(Y, i, 1, WEST_EDGE)  = -gnij + i20             ! y min
        edges(Z, i, 1, NORTH_EDGE) = gnij                    ! z = lnij
        edges(Z, i, 1, SOUTH_EDGE) = gnij 
        edges(Z, i, 1, EAST_EDGE)  = gnij 
        edges(Z, i, 1, WEST_EDGE)  = gnij 
      enddo
  
    case(5)
      do i = 1, lnij
        edges(X, i, 1, NORTH_EDGE) = -gnij + j20 + 2*lnij    ! x max
        edges(X, i, 1, SOUTH_EDGE) = -gnij + j20             ! x min
        edges(X, i, 1, EAST_EDGE)  = -gnij + j20 + 2*i - 1   ! +x axis
        edges(X, i, 1, WEST_EDGE)  = -gnij + j20 + 2*i - 1   ! +x axis
        edges(Y, i, 1, NORTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, SOUTH_EDGE) = -gnij + i20 + 2*i - 1   ! +y axis
        edges(Y, i, 1, EAST_EDGE)  = -gnij + i20 + 2*lnij    ! y max
        edges(Y, i, 1, WEST_EDGE)  = -gnij + i20             ! y min
        edges(Z, i, 1, NORTH_EDGE) = -gnij                   ! z = -lnij
        edges(Z, i, 1, SOUTH_EDGE) = -gnij 
        edges(Z, i, 1, EAST_EDGE)  = -gnij 
        edges(Z, i, 1, WEST_EDGE)  = -gnij 
      enddo
  
    case DEFAULT
      edges = -999999

  end select

  return
end subroutine test_edges_create

