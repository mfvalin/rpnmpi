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
!                         GEM solver with FFTW memory layouts
!
!
!                         ===    with x and y transposes  ===
!
!                         dimension(lnix,      lnjy,      gnk)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose X                     inverse transpose X
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lnix,  lnjy,  lnkx , npex)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
! forward FFT(x) <------> dimension(gni,      lnjy,      lnkx)  <------> inverse FFT(x)
! normalisation           |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lniy,  lnjy,  lnkx,  npey)  [ all i indices on same PE ]
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose Y                     inverse transpose Y
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lniy,  lnjy,  lnkx,  npey)  [ all j indices on same PE ]
!                                            ^
!                                            |
!                                            v
!                                 tridiagonal solver along j
!                         
!
!                     ===  with x transpose and ripple solver  ===
!
!                         dimension(mini:maxi, minj:maxj, gnk) [ or (lnix,lnjy,gnk) ]
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lnix,      lnjy,      gnk)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose X                     inverse transpose X
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lnix,  lnjy,  lnkx , npex)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
! forward FFT(x) <------> dimension(gni,      lnjy,      lnkx)  <------> inverse FFT(x)
!                             [ SHARED MEMORY on NUMA node]
!                                            ^
!                                            |
!                                            v
!                     ripple distributed tridiagonal solver along j
module RPN_MPI_fast_transpose_cache
  use rpn_mpi_mpif
  implicit none
  integer, parameter :: MAX_PEX = 128
  integer, dimension(MAX_PEX), save :: sendcount, senddispl, recvcount, recvdispl, countnk
  integer, save :: nk = 0
  integer, save :: npex = 0
  integer, save :: npey = 0
  integer, save :: rankx = -1
  integer, save :: rowcom = MPI_COMM_NULL
  integer, save :: colcom = MPI_COMM_NULL
end module RPN_MPI_fast_transpose_cache
!
! transpose along x
!
! in the original array a given PE has 
! all slices(consecutive in memory) along z (all k indices, gnk values)
! one slice along y (length lnjy)
! one slice along x (length lnix)
!
! in the transposed array, a given PE has
! one slice along z (length lnkx)
! one slice along y (length lnjy)
! all slices along x (all i indices, gni values)
!
! forward transpose (forward == .true.)
! z         : original array, dimension(lnix,lnjy,gnk)
! zt        : transposed arrray, dimension(lnix,lnjy,lnkx,npex)
! lnix      : number of local points along x (assumed to be IDENTICAL on ALL PEs in row)
! lnjy      : number of local points along y (assumed to be IDENTICAL on ALL PEs in row)
! gnk       : total number of levels
! lnkx      : local number of levels in transposed array zt (the SUM of lnkx MUST be equal to gnk)
! row_comm  : communicator for this transpose
! ierr      : error flag, same as MPI
!
 subroutine RPN_MPI_fast_transpose_x(z, zt, forward, lnix, lnjy, gnk, lnkx, row_comm, ierr) !InTf!
  use ISO_C_BINDING
  use RPN_MPI_fast_transpose_cache
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lnix, lnjy, gnk, lnkx                  !InTf!
  integer, intent(IN) :: row_comm                               !InTf!
! little white lie in interface, z, zt are advertised as addresses passed by value
!! type(RPN_MPI_Loc), intent(IN), value :: z, zt                !InTf!
  integer, dimension(lnix,lnjy,gnk), intent(INOUT)    :: z      ! NO HALO in arrays 
  integer, dimension(lnix,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npex
  integer, intent(OUT) :: ierr                                  !InTf!

  integer :: n
  integer :: i, j, k, m

  ierr = MPI_ERROR
  if(npex == 0 .or. row_comm .ne. rowcom .or. nk .ne. gnk) then

    rowcom = row_comm
    call MPI_Comm_rank(rowcom, rankx, ierr)
    call MPI_Comm_size(rowcom, npex, ierr)
    if(npex > MAX_PEX) goto 111                        ! ERROR, tables are too small
    call MPI_Allgather(lnkx, 1, MPI_INTEGER, countnk, 1, MPI_INTEGER, rowcom, ierr) ! get lnk counts
    if(sum(countnk(1:npex)) .ne. gnk) goto 222         ! ERROR, missing or extra levels
  endif
! write(0,*) 'sendcount = ',sendcount(1:npex)
! write(0,*) 'lnix, lnjy, gnk, lnkx =',lnix, lnjy, gnk, lnkx
  if(forward) then                                     ! z -> zt
    sendcount(1:npex) = countnk(1:npex) * lnix * lnjy    ! may be zero for some PEs
    recvcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    senddispl(1) = 0
    recvdispl(1) = 0
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
! write(0,*) 'sendcount = ',sendcount(1:npex)
! write(0,*) 'senddispl = ',senddispl(1:npex)
! write(0,*) 'recvcount = ',recvcount(1:npex)
! write(0,*) 'recvdispl = ',recvdispl(1:npex)
! write(0,*) 'z'
! do j = lnjy,1,-1
!   write(0,2) (z(1:lnix,j,k),k=1,gnk)
! enddo
    call MPI_Alltoallv(z , sendcount, senddispl, MPI_INTEGER,               &
                       zt, recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
! write(0,*) 'zt'
! do j = lnjy,1,-1
!   write(0,2) (zt(1:lnix,j,1,m),m=1,npex)
! enddo
  else                                                 ! zt -> z
    sendcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    recvcount(1:npex) = countnk(1:npex) * lnix * lnjy    ! may be zero for some PEs
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    call MPI_Alltoallv(zt, sendcount, senddispl, MPI_INTEGER,               &
                       z , recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
  endif

1 return
111 continue
  write(0,*) 'ERROR: more than',MAX_PEX,' levels'
  goto 1
222 continue
  write(0,*) 'ERROR: the sum of lnk() does not match gnk'
  goto 1
2 format(30I10.9)
 end subroutine RPN_MPI_fast_transpose_x !InTf!
!
! transpose along y, the original and transposed arrays have the same dimensions
!
! in the original array a given PE has
!
! all slices along x (all i indices, gni values)
! one slice along y (length lnjy)
! one slice along z (length lnkx)
!
! in the transposed array, a given PE has
! one slice along x (length lniy)
! all slices along y (all j indices, gnj values)
! one slice along z (length lnkx)
!
! forward transpose (forward == .true.)
! z         : original array, dimension(lniy,lnjy,lnkx,npey)
! zt        : transposed arrray, dimension(lniy,lnjy,lnkx,npey)
! lniy      : number of local points along x (assumed to be IDENTICAL on ALL PEs in column)
! lnjy      : number of local points along y (assumed to be IDENTICAL on ALL PEs in column)
! lnkx      : local number of levels in transposed array zt (the SUM of lnkx MUST be equal to gnk)
! col_comm  : communicator for this transpose
! ierr      : error flag, same as MPI
!
 subroutine RPN_MPI_fast_transpose_y(z, zt, forward, lniy, lnjy, lnkx, col_comm, ierr) !InTf!
  use ISO_C_BINDING
  use RPN_MPI_fast_transpose_cache
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lniy, lnjy, lnkx                  !InTf!
  integer, intent(IN) :: col_comm                               !InTf!
! little white lie in interface, z, zt are advertised as addresses passed by value
!! type(RPN_MPI_Loc), intent(IN), value :: z, zt                !InTf!
  integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: z      ! NO HALO in arrays 
  integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npey
  integer, intent(OUT) :: ierr                                  !InTf!

  integer :: n
  integer :: i, j, k, m

  ierr = MPI_ERROR
  if(npey == 0 .or. col_comm .ne. colcom) then
    colcom = col_comm
    call MPI_Comm_size(colcom, npey, ierr)
  endif

  if(forward) then                                     ! z -> zt
    call MPI_Alltoall(z , lniy*lnjy*lnkx, MPI_INTEGER,               &
                      zt, lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
  else                                                 ! zt -> z
    call MPI_Alltoall(zt, lniy*lnjy*lnkx, MPI_INTEGER,               &
                      z , lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
  endif

1 return
111 continue
  write(0,*) 'ERROR: more than',MAX_PEX,' levels'
  goto 1
222 continue
  write(0,*) 'ERROR: the sum of lnk() does not match gnk'
  goto 1
2 format(30I10.9)
 end subroutine RPN_MPI_fast_transpose_y !InTf!



