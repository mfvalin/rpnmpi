!/! RPN_MPI - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
!====================================================================================================
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
!****P* rpn_mpi/halo  halo exchange routines
! DESCRIPTION
! these routines are used to perform a halo exchange between members of a "grid"
! an internal timing package keeps timing statistics about the various phases of the exchnage
!
!  +=======================================================================================+___maxy
!  I       <halox>                                                           <halox>       I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I :     :                                    |                                  :     : I
!  I :  C  :                                  haloy                                :  C  : I
!  I :     :                                    |                                  :     : I
!  I +-----+=====+===========================================================+=====+-----: I___nj
!  I :     I     :                              |                            :     I     : I
!  I :     I     :                            haloy                          :     I     : I
!  I :     I     :                              |                            :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                              |                            :     I     : I
!  I :     I     :                            haloy                          :     I     : I
!  I :     I     :                              |                            :     I     : I
!  I +-----+=====+===========================================================+=====+-----: I___1
!  I :     :                                    |                                  :     : I
!  I :  C  :                                  haloy                                :  C  : I
!  I :     :                                    |                                  :     : I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I <halox>                                                                       <halox> I
!  +=======================================================================================+___miny
!  |       |                                                                       |       |
!  |       |                                                                       |       |
!  |       1                                                                       ni      |
!  minx                                                                                 maxx
!
!  the "useful" section of the array is array(i,j) where 1 <= i <= ni, 1 <= j < nj
!  the "inner" halo is the halo part inside the "useful" section of the array (1:ni,1:nj)
!  the "outer" halo is the halo part outside the "useful" section of the array (1:ni,1:nj)
!  part of the array may be outside of the "outer" halo. that part will be left as is
!
!  a East-West (x direction) exchange is performed first, for rows 1 thru nj
!  a North-South (y direction) exchange will follow, for columns 1-halox thru ni+halox
!
!  the corner parts (C) get implicitly exchanged with the corner neighbors during tne North-South phase
!  (this achieves the 2D exchange with 8 neighbors using only 4 messages)
!
! example of usage
!
!   include 'RPN_MPI.inc'
! ! macro names are CASE SENSITIVE
! #define LoC(what) rpn_mpi_loc(loc(what))
! real, dimension(:,:,:), allocatable :: z
! type(RPN_MPI_Comm) :: col_comm, row_comm
! integer :: halox, haloy, NI,NJ,NK
!
! allocate (z(1-halox:NI+halox,1-haloy:NJ+haloy,NK))
!
! call RPN_MPI_ez_halo_parms(row_comm, col_comm, 'BARRIER')
! call RPN_MPI_ez_halo(LoC(z),1-halox,NI+halox,1-haloy,NJ+haloy,NI,NJ,NK,halox,haloy)
! call RPN_MPI_halo(LoC(z),1-halox,NI+halox,1-haloy,NJ+haloy,NI,NJ,NK,halox,haloy,row_comm,col_comm)
!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
module RPN_MPI_halo_cache
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
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
!****f* rpn_mpi/RPN_MPI_halo halo exchange for 4 byte items
! DESCRIPTION
! this routine will perform a horizontal halo exchange with NO PERIODIC BOUNDARY conditions 
! for all "planes" of a "grid"
!
! row         : RPN_MPI communicator for the grid row this PE belongs to
! col         : RPN_MPI communicator for the grid column this PE belongs to
! minx, maxx, miny, maxy : horizontal(plane) dimensions of array pointed to by g
! nk          : vertical dimension of array pointed to by g
! lni, lnj    : "useful" horizontal dimensions of array pointed to byy g
! halox       : number of "halo" points along x
! haloy       : number of "halo" points along y
! g           : address (wrapped) of array
!
! SYNOPSIS
 subroutine RPN_MPI_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,row,col) BIND(C,name='RPN_MPI_halo') !InTf!
! IGNORE
  implicit none
!! import :: RPN_MPI_Loc, RPN_MPI_Comm, C_INT                                          !InTf!
! ARGUMENTS
  integer(C_INT), intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Comm), intent(IN) :: row,col                                           !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                           !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
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

!****f* rpn_mpi/RPN_MPI_halo8 8 byte version of RPN_MPI_halo
! DESCRIPTION
! see RPN_MPI_halo. (array subject to halo exchange is made of 8 byte items)
!
! SYNOPSIS
 subroutine RPN_MPI_halo8(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,row,col) BIND(C,name='RPN_MPI_halo8')  !InTf!
! IGNORE
  use ISO_C_BINDING
  implicit none
!! import :: RPN_MPI_Loc, RPN_MPI_Comm, C_INT                                          !InTf!
! ARGUMENTS
  integer(C_INT), intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Comm), intent(IN) :: row,col                                           !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                           !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
! white lie in published interface, g is published as an address passed by value
  integer, intent(INOUT), dimension(2*minx-1:2*maxx,miny:maxy,nk) :: g
  integer, intent(IN)    :: row,col

  call RPN_MPI_halo(g,2*minx-1,2*maxx,miny,maxy,2*lni,lnj,nk,2*halox,haloy,row,col)
  return
 end subroutine RPN_MPI_halo8                                                           !InTf!

!****f* rpn_mpi/RPN_MPI_ez_halo grid halo exchange for 4 byte items
! DESCRIPTION
! same functionality as RPN_MPI_halo, communicators are implicit
! (a previous call to RPN_MPI_ez_halo_parms may be needed)
!
! SYNOPSIS
 subroutine RPN_MPI_ez_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy) bind(C,name='RPN_MPI_ez_halo')       !InTf!
! IGNORE
  use ISO_C_BINDING
  implicit none
!!  import :: RPN_MPI_Loc                                                       !InTf!
! ARGUMENTS
  integer, intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
! white lie in published interface, g is published as an address passed by value
  integer, intent(INOUT), dimension(minx:maxx,miny:maxy,nk) :: g

  if(rowcom == MPI_COMM_NULL .or. colcom == MPI_COMM_NULL) then
    ! get the appropriate row and column communicators from the internal RPN_MPI data
  endif

  call RPN_MPI_halo(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy,rowcom,colcom)
  return
 end subroutine RPN_MPI_ez_halo                                                  !InTf!

!****f* rpn_mpi/RPN_MPI_ez_halo_8 8 byte version of RPN_MPI_ez_halo
! DESCRIPTION
! same functionality as RPN_MPI_halo8, communicators are implicit
! (a previous call to RPN_MPI_ez_halo_parms may be needed)
!
! SYNOPSIS
 subroutine RPN_MPI_ez_halo8(g,minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy) bind(C,name='RPN_MPI_ez_halo8')     !InTf!
! IGNORE
  use ISO_C_BINDING
  implicit none
!!  import :: RPN_MPI_Loc                                                       !InTf!
! ARGUMENTS
  integer, intent(IN)    :: minx,maxx,miny,maxy,lni,lnj,nk,halox,haloy          !InTf!
!! type(RPN_MPI_Loc), intent(IN), value :: g                                   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
! white lie in published interface, g is published as an address passed by value
  integer, intent(INOUT), dimension(2*minx-1:2*maxx,miny:maxy,nk) :: g
  call RPN_MPI_ez_halo(g,2*minx-1,2*maxx,miny,maxy,2*lni,lnj,nk,2*halox,haloy)
 end subroutine RPN_MPI_ez_halo8                                                !InTf!
end module RPN_MPI_halo_cache
!====================================================================================================

!****f* rpn_mpi/RPN_MPI_get_halo_timings get a copy of the current timing stats
! DESCRIPTION
! t      : array of 64 bit integers to receive timing stats
! n      : dimension of array t
! function value : number of halo exchanges performed since last stats reset
! SYNOPSIS
 function RPN_MPI_get_halo_timings(t,n) result(nt) BIND(C,name='RPN_MPI_get_halo_timings')   !InTf!
! IGNORE
  use RPN_MPI_halo_cache
  implicit none
! ARGUMENTS
  integer, intent(IN) :: n                                                                  !InTf!
  integer(kind=8), dimension(n), intent(OUT) :: t                                           !InTf!
  integer :: nt                                                                             !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
  t(1:n) = -1
  nt = -1
  if(n < NTS) return
  nt = nhalo
  t(1:NTS) = ts
  return
 end function RPN_MPI_get_halo_timings                                                       !InTf!

!****f* rpn_mpi/RPN_MPI_reset_halo_timings reset timing stats
! SYNOPSIS
 subroutine RPN_MPI_reset_halo_timings() BIND(C,name='RPN_MPI_reset_halo_timings')   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
  use RPN_MPI_halo_cache
  implicit none
  nhalo = 0
  ts = 0
  return
 end subroutine RPN_MPI_reset_halo_timings                                           !InTf!

!****f* rpn_mpi/RPN_MPI_print_halo_timings print accumulated timing stats
! SYNOPSIS
 subroutine RPN_MPI_print_halo_timings() BIND(C,name='RPN_MPI_print_halo_timings')   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
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
!****f* rpn_mpi/RPN_MPI_ez_halo_parms setup call for ez_halo routines
! DESCRIPTION
! row         : RPN_MPI communicator for the grid row this PE belongs to
! col         : RPN_MPI communicator for the grid column this PE belongs to
! mode == 'B' : add a MPI_Barrier call between the x and y phases of the exchange
!               (used for timing purposes)
! SYNOPSIS
 subroutine RPN_MPI_ez_halo_parms(row, col, mode) bind(C,name='RPN_MPI_ez_halo_parms')     !InTf!
! IGNORE
  use RPN_MPI_halo_cache
  implicit none
!!  import :: RPN_MPI_Comm                                  !InTf!
! ARGUMENTS
!! type(RPN_MPI_Comm), intent(IN) :: row, col               !InTf!
  integer, intent(IN) :: row, col
  character(len=1), dimension(*), intent(IN) :: mode        !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
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
