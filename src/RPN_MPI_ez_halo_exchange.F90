! * Copyright (C) 2020  Recherche en Prevision Numerique
! *
! * This software is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This software is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
!====================================================================================================
module RPN_MPI_halo_cache
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
!   include 'mpif.h'
  include 'RPN_MPI_system_interfaces.inc'
  integer, parameter :: NTS = 9
  integer(kind=8), dimension(NTS) :: ts = [0,0,0,0,0,0,0,0,0]   ! timing points
!   logical, save :: redblack = .false.  ! if true, use fully synchronous (send, recv, no sendrecv) method
!   logical, save :: async    = .false.  ! if true, use fully asynchronous (isend, irecv)  method
  logical, save :: barrier  = .false.  ! if true, use a barrier between E-W and N-S exchanges
  integer, save :: nhalo = 0
  integer, save :: rowcom  = MPI_COMM_NULL
  integer, save :: colcom  = MPI_COMM_NULL
  integer, save :: rankx, ranky, sizex, sizey
contains
! an interface will be needed to use RPN_MPI_halo
! the published interface will use the "void *" approach
!
 subroutine RPN_MPI_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,row,col) BIND(C,name='RPN_MPI_halo') !InTf!
  implicit none
!! import :: RPN_MPI_Loc, RPN_MPI_Comm, C_INT                                          !InTf!
  integer(C_INT), intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Comm), intent(IN) :: row,col                                           !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                           !InTf!
! white lie in published interface, g is published as an address passed by value
  integer, intent(INOUT), dimension(minx:maxx,miny:maxy,nk) :: g
  integer, intent(IN)    :: row,col

  integer :: i, j, k, nw, ier
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(halox,lnj,nk) :: halo_from_west, halo_to_west
  integer, dimension(halox,lnj,nk) :: halo_from_east, halo_to_east
  integer, dimension(halox+lni+halox,haloy,nk) :: halo_from_north, halo_to_north
  integer, dimension(halox+lni+halox,haloy,nk) :: halo_from_south, halo_to_south
  integer(kind=8),dimension(10) :: t
!   integer(kind=8), external :: cpu_real_time_ticks

  if(row .ne. rowcom .or. col .ne. colcom) then  ! different set of row / column
    rowcom = row
    call MPI_Comm_rank(rowcom, rankx, ier)
    call MPI_Comm_size(rowcom, sizex, ier)
    colcom = col
    call MPI_Comm_rank(colcom, ranky, ier)
    call MPI_Comm_size(colcom, sizey, ier)
  endif
  t = 0
  nhalo = nhalo + 1
  t(1) = cpu_real_time_ticks()
  t(2) = t(1)
  t(3) = t(1)
  t(4) = t(1)
  if(sizex .gt. 1) then   ! is there an exchange along x ?
    nw = halox * lnj * nk       ! message length
    do k = 1, nk          ! peel west and east side of array simultaneously
    do j = 1, lnj
      halo_to_west(:,j,k) = g(1:halox        ,j,k)    ! extract west side inner halo
      halo_to_east(:,j,k) = g(lni+1-halox:lni,j,k)    ! extract east side inner halo
    enddo
    enddo
    t(2) = cpu_real_time_ticks()
    if(rankx .eq. 0)then             ! west PE, send to east, get from east
      call MPI_SENDRECV(halo_to_east  , nw, MPI_INTEGER, rankx+1, rankx, &
			halo_from_east, nw, MPI_INTEGER, rankx+1, rankx+1, &
			rowcom, STATUS, ier)
    elseif(rankx .lt. sizex-1) then  ! middle PE, (get from west, send to east) (get from east, send to west)
      call MPI_SENDRECV(halo_to_east  , nw, MPI_INTEGER, rankx+1, rankx, &
			halo_from_west, nw, MPI_INTEGER, rankx-1, rankx-1, &
			rowcom, STATUS, ier)
      call MPI_SENDRECV(halo_to_west  , nw, MPI_INTEGER, rankx-1, rankx, &
			halo_from_east, nw, MPI_INTEGER, rankx+1, rankx+1, &
			rowcom, STATUS, ier)
    else                             ! east PE, send to west, get from west
      call MPI_SENDRECV(halo_to_west  , nw, MPI_INTEGER, rankx-1, rankx, &
			halo_from_west, nw, MPI_INTEGER, rankx-1, rankx-1, &
			rowcom, STATUS, ier)
    endif
    t(3) = cpu_real_time_ticks()
    do k = 1, nk          ! insert west and east side of array simultaneously
    do j = 1, lnj
      if(rankx .gt. 0)       g(1-halox:0,j,k)       = halo_from_west(:,j,k)     ! insert into west side outer halo
      if(rankx .lt. sizex-1) g(lni+1:lni+halox,j,k) = halo_from_east(:,j,k)     ! insert into east side outer halo
    enddo
    enddo
    t(4) = cpu_real_time_ticks()
  endif
  if(barrier) call MPI_Barrier(colcom, ier)
  t(5) = cpu_real_time_ticks()
  t(6) = t(5)
  t(7) = t(5)
  t(8) = t(5)
  if(sizey .gt. 1) then   ! is there an exchange along y ?
    nw = (lni + 2*halox) * haloy * nk
    do k = 1, nk          ! peel north and south side of array
    do j = 1, haloy
      halo_to_north(:,j,k) = g(1-halox:lni+halox,lnj-haloy+j,k)    ! extract north side inner halo
      halo_to_south(:,j,k) = g(1-halox:lni+halox,j          ,k)    ! extract south side inner halo
    enddo
    enddo
    t(6) = cpu_real_time_ticks()
    if(ranky .eq. 0) then                ! south PE, send to north, receive from north
      call MPI_SENDRECV(halo_to_north  , nw, MPI_INTEGER, ranky+1, ranky, &
                        halo_from_north, nw, MPI_INTEGER, ranky+1, ranky+1, &
                        colcom, STATUS, ier)
    elseif(ranky .lt. sizey-1) then      ! middle PE,(send to north, get from south) (send to south, get from north)
      call MPI_SENDRECV(halo_to_north  , nw, MPI_INTEGER, ranky+1, ranky, &
                        halo_from_south, nw, MPI_INTEGER, ranky-1, ranky-1, &
                        colcom, STATUS, ier)
      call MPI_SENDRECV(halo_to_south  , nw, MPI_INTEGER, ranky-1, ranky, &
                        halo_from_north, nw, MPI_INTEGER, ranky+1, ranky+1, &
                        colcom, STATUS, ier)
    else                                 ! north PE, send to south, receive from south
      call MPI_SENDRECV(halo_to_south  , nw, MPI_INTEGER, ranky-1, ranky, &
                        halo_from_south, nw, MPI_INTEGER, ranky-1, ranky-1, &
                        colcom, STATUS, ier)
    endif
    t(7) = cpu_real_time_ticks()
    do k = 1, nk          ! peel north and south side of array
    do j = 1, haloy
      if(ranky .lt. sizey-1) g(1-halox:lni+halox,lnj+j  ,k) = halo_from_north(:,j,k)  ! not north row insert into north outer halo
      if(ranky .gt. 0)       g(1-halox:lni+halox,j-haloy,k) = halo_from_south(:,j,k)  ! not south row insert into south outer halo
    enddo
    enddo
    t(8) = cpu_real_time_ticks()
  endif
  ts(1:3) = ts(1:3) + t(2:4)-t(1:3)
  ts(4:6) = ts(4:6) + t(6:8)-t(5:7)
  ts(7) = ts(7) + t(4) - t(1)
  ts(8) = ts(8) + t(8) - t(5)
  ts(9) = t(5) - t(4)
  return
 end subroutine RPN_MPI_halo                                               !InTf!

 subroutine RPN_MPI_halo8(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,row,col) BIND(C,name='RPN_MPI_halo8')  !InTf!
  use ISO_C_BINDING
  implicit none
!! import :: RPN_MPI_Loc, RPN_MPI_Comm, C_INT                                          !InTf!
  integer(C_INT), intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Comm), intent(IN) :: row,col                                           !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                           !InTf!
! white lie in published interface, g is published as an address passed by value
  integer, intent(INOUT), dimension(2*minx-1:2*maxx,miny:maxy,nk) :: g
  integer, intent(IN)    :: row,col

  call RPN_MPI_halo(g,2*minx-1,2*maxx,miny,maxy,2*lni,lnj,nk,2*halox,haloy,row,col)
  return
 end subroutine RPN_MPI_halo8                                                           !InTf!

 subroutine RPN_MPI_ez_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy) bind(C,name='RPN_MPI_ez_halo')       !InTf!
  use ISO_C_BINDING
  implicit none
!!  import :: RPN_MPI_Loc                                                       !InTf!
  integer, intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
! white lie in published interface, g is published as an address passed by value
!!  type(RPN_MPI_Loc), intent(IN), value :: g                                   !InTf!
  integer, intent(INOUT), dimension(minx:maxx,miny:maxy,nk) :: g

  if(rowcom == MPI_COMM_NULL .or. colcom == MPI_COMM_NULL) then
    ! get the appropriate row and column communicators from the internal RPN_MPI data
  endif

  call RPN_MPI_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,rowcom,colcom)
  return
 end subroutine RPN_MPI_ez_halo                                                  !InTf!

 subroutine RPN_MPI_ez_halo_8(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy) bind(C,name='RPN_MPI_ez_halo_8')     !InTf!
  use ISO_C_BINDING
  implicit none
!!  import :: RPN_MPI_Loc                                                       !InTf!
  integer, intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
! white lie in published interface, g is published as an address passed by value
!!  type(RPN_MPI_Loc), intent(IN), value :: g                                   !InTf!
  integer, intent(INOUT), dimension(2*minx-1:2*maxx,miny:maxy,nk) :: g
  call RPN_MPI_ez_halo(g,2*minx-1,2*maxx,miny,maxy,2*lni,lnj,nk,2*halox,haloy)
 end subroutine RPN_MPI_ez_halo_8                                                !InTf!
end module RPN_MPI_halo_cache
!====================================================================================================

 function RPN_MPI_get_halo_timings(t,n) result(nt) BIND(C,name='RPN_MPI_get_halo_timings')   !InTf!
  use RPN_MPI_halo_cache
  implicit none
  integer, intent(IN) :: n                                                                  !InTf!
  integer(kind=8), dimension(n), intent(OUT) :: t                                           !InTf!
  integer :: nt                                                                             !InTf!
  t(1:n) = -1
  nt = -1
  if(n < NTS) return
  nt = nhalo
  t(1:NTS) = ts
  return
 end function RPN_MPI_get_halo_timings                                                       !InTf!

 subroutine RPN_MPI_reset_halo_timings() BIND(C,name='RPN_MPI_reset_halo_timings')   !InTf!
  use RPN_MPI_halo_cache
  implicit none
  nhalo = 0
  ts = 0
  return
 end subroutine RPN_MPI_reset_halo_timings                                           !InTf!

 subroutine RPN_MPI_print_halo_timings() BIND(C,name='RPN_MPI_print_halo_timings')   !InTf!
  use RPN_MPI_halo_cache
  implicit none
  integer(kind=8), dimension(:,:), allocatable :: alltsx, allts
  integer :: ier, ntot, i
  integer(kind=8), dimension(NTS) :: tsavg, tsmax, tsmin

  if(nhalo == 0) return   ! no stats collected

  ntot = sizex*sizey
  if(rankx == 0) then
    allocate(alltsx(NTS,sizex))   ! rank 0 in row will collect along x
    if(ranky == 0) then
      allocate(allts(NTS,ntot))   ! rank 0 in x and y will collect the full grid
    else
      allocate(allts(NTS,1))
    endif
  else
    allocate(alltsx(NTS,1))
    allocate(allts(NTS,1))
  endif
  call MPI_Gather(ts,     NTS,       MPI_INTEGER8, alltsx, NTS,       MPI_INTEGER8, 0, rowcom, ier)  ! collect along row
  if(rankx == 0)call MPI_Gather(alltsx, NTS*sizex, MPI_INTEGER8, allts,  NTS*sizex, MPI_INTEGER8, 0, colcom, ier)  ! collect all rows
  do i = 1, NTS
    tsavg(i) = sum(allts(i,:))
    tsmax(i) = maxval(allts(i,:))
    tsmin(i) = minval(allts(i,:))
  enddo
  tsavg = tsavg / ntot
  if(rankx == 0 .and. ranky == 0) then
    write(6,1)'nhalo, avg times =     nexch    peel-J    mesg-x    plug-J    peel-I    mesg-y    plug-I    exch-x    exch-y    barrier'
    write(6,1)'nhalo, avg times =',nhalo,tsavg/nhalo
    write(6,1)'min times        =          ',tsmin/nhalo
    write(6,1)'max times        =          ',tsmax/nhalo
    call flush(6)
  endif
  deallocate(allts)
  deallocate(alltsx)
  return
  
1 format(a,10I10)
 end subroutine RPN_MPI_print_halo_timings                                           !InTf!

! set information for halo exchange timings
 subroutine RPN_MPI_ez_halo_parms(row, col, mode) bind(C,name='RPN_MPI_ez_halo_parms')     !InTf!
  use RPN_MPI_halo_cache
  implicit none
!!  import :: RPN_MPI_Comm                                  !InTf!
!!  type(RPN_MPI_Comm), intent(IN) :: row, col              !InTf!
  integer, intent(IN) :: row, col
  character(len=1), dimension(*), intent(IN) :: mode        !InTf!
  integer :: ier
  rowcom  = row   ! row communicator
  colcom  = col   ! column communicator
  call MPI_comm_size(rowcom, sizex, ier)  ! size of row
  call MPI_comm_rank(rowcom, rankx, ier)  ! rank in row
  call MPI_comm_size(colcom, sizey, ier)  ! size of column
  call MPI_comm_rank(colcom, ranky, ier)  ! rank in column
!   redblack = trim(mode) .eq. 'REDBLACK'
!   async    = trim(mode) .eq. 'ASYNC'
  barrier  = mode(1) .eq. 'B'  ! for 'BARRIER'
  if(rankx+ranky == 0 .and. barrier) write(6,*) 'MODE = BARRIER'
  return
 end subroutine RPN_MPI_ez_halo_parms                        !InTf!

! RED/BLACK method uses only synchronous send and recv
!   West -> East move
!     EVEN PEs
!       if not East PE send to odd PE at East 
!       if not West PE recv from odd PE at West 
!     ODD PEs
!                      recv from even PE at West
!       if not East PE send to even PE at East
!   West <- East move
!     EVEN PEs
!       if not East PE recv from odd PE at East
!       if not West PE send to odd PE at West
!     ODD PEs
!                      send to even PE at West
!       if not East PE recv from even PE at East
