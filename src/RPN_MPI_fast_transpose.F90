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
  integer, save :: myrank = -1
  integer, save :: rowcom = MPI_COMM_NULL
end module RPN_MPI_fast_transpose_cache
!
! forward transpose (forward == .true.)
! z         : original array, dimension(lnix,lnjy,gnk)
! zt        : transposed arrray, dimension(lnix,lnjy,lnkx,npex)
! lnix      : number of local points along x (assumed to be IDENTICAL on ALL PEs)
! lniy      : number of local points along y (assumed to be IDENTICAL on ALL PEs)
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
  integer, dimension(lnix,lnjy,gnk), intent(INOUT)    :: z      ! NO HALO in local arrays 
  integer, dimension(lnix,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npex
  integer, intent(OUT) :: ierr                                  !InTf!

  integer :: n

  ierr = MPI_ERROR
  if(npex == 0 .or. row_comm .ne. rowcom .or. nk .ne. gnk) then
    rowcom = row_comm
    call MPI_Comm_rank(rowcom, myrank, ierr)
    call MPI_Comm_size(rowcom, npex, ierr)
    if(npex > MAX_PEX) goto 111                        ! ERROR, tables are too small
    call MPI_Allgather(lnkx, 1, MPI_INTEGER, countnk, 1, MPI_INTEGER, rowcom, ierr) ! get lnk counts
    if(sum(countnk(1:npex)) .ne. gnk) goto 222         ! ERROR, missing or extra levels
  endif

  if(forward) then                                     ! z -> zt
    sendcount(1:npex) = countnk(1:nk) * lnix * lnjy    ! may be zero for some PEs
    recvcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    senddispl(1) = 0
    recvdispl(1) = 0
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    call MPI_Alltoallv(z , sendcount, senddispl, MPI_INTEGER,               &
                       zt, recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
  else                                                 ! zt -> z
    sendcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    recvcount(1:npex) = countnk(1:nk) * lnix * lnjy    ! may be zero for some PEs
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
 end subroutine RPN_MPI_fast_transpose_x !InTf!



