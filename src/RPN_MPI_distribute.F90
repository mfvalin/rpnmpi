!/* RPN_MPI - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! */
!**function  RPN_data_distribute global domain decomposition function (along one dimension)
!InTf!
!****f* RPN_COMM/RPN_data_distribute
! SYNOPSIS
 subroutine RPN_data_distribute(myrank, npe, gmin, gmax, lmin, lmax, count, offset, mode)    !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
      implicit none                                                !InTf!
!
! ARGUMENTS
!  IN   myrank "tile" ordinal along decomposition axis (0 -> npe-1)
!  IN   npe    number of "tiles" (PEs) along this dimension
!  IN   gmin,gmax
!              global index space along this dimension is gmin:gmax (usually 1:n)
!  OUT  lmin,lmax
!              this "tile" will cover index range lmin:lmax in global space
!              lmax << lmin signals a distribution error
!              lmax = lmin -1 is OK, it happens for ZERO length tiles at end
!  OUT  count(1:npe)
!              count(i) = number of points along this dimension for PE with ordinal I-1
!  OUT  offset(1:npe)
!              offset(i) = offset from gmin for PE with ordinal I-1
!  IN   mode   decomposition mode, OPTIONAL argument (0 by default)
!          0 : strict mode, all tiles but last one must have same dimension, 
!              last tile may be shorter but may not have zero dimension
!          1 : some tiles at end may be 1 shorter than tiles at beginning, 
!              zero size not allowed for these shorter tiles
!          2 : same as mode=1 but zero dimension tiles at end are allowed 
!              (useful only if more PEs than points)
!          3 : tiles with same length possibly followed by ONE shorter tile 
!              possibly followed by ONE or MORE zero size tiles
!          4 : spread as evenly as possible, consecutive tiles may not have same length
!
! Notes
!     this function is totally stand alone and could eventually be moved into the rmnlib library
!     mode 2 only makes sense when one has more PEs than points
!******
!*
  integer, intent(IN) ::  myrank, npe, gmin, gmax               !InTf!
  integer, intent(OUT) :: lmin, lmax                            !InTf!
  integer, intent(OUT) :: count(npe), offset(npe)               !InTf!
  integer, intent(IN), optional ::  mode                        !InTf!

  integer :: the_mode, total, run, i, npts

  count  = 0
  total  = 0
  offset = 0
  lmin   = 0
  lmax   = lmin - 999999                    ! error flag, lmax WAY smaller than lmin
  npts   = gmax - gmin + 1                  ! number of data points to distribute
  if(myrank >= npe .or. myrank < 0) return  ! ERROR, impossible rank value
  if(npts <= 0) return                      ! ERROR, invalid number of points
  if(npe <= 0) return                       ! ERROR, invalid number of "tiles" (PEs)

  the_mode = 0
  if(present(mode)) the_mode = mode
  SELECT CASE(the_mode)

  CASE(0,3)                                 ! strict mode, all tiles but last must have same dimension (0)
                                            ! tile(s) with zero length allowed at end  (3)
    run = (npts + npe - 1) / npe            ! ceiling(npts / npe)
    do i = 1 , npe
      offset(i) = total
      count(i)  = min(npts - total , run)   ! 
      total     = min(npts, total + run)
    enddo
    if(count(npe) <= 0 .and. the_mode == 0) return  ! ERROR, ZERO length tile(s) at end not allowed in mode 0

  CASE(1,2)                                 ! some tiles at end may be 1 shorter than tiles at beginning (1)
                                            ! zero dimension tile(s) at end allowed (2) (only happens when npe > npts)
    run = npts / npe                        ! floor(npts / npe)  run will be zero if npe > npts (OK for mode 2)
    count = run
    do i = 1 , npts - run * npe             ! some tiles at beginning have to be ONE longer
      count(i) = count(i) + 1
    enddo
    if(count(npe) <= 0 .and. the_mode == 1) return  ! ERROR, ZERO length tile(s) at end not allowed in mode 1
    do i = 2 , npe
      offset(i) = offset(i-1) + count(i-1)
    enddo

  CASE(4)                                   ! spread as evenly as possible, consecutive tiles may have different length
    run = -npe/2
    do i = 1, npe
      run = run + npts
      offset(i) = total
      do while (run > 0)
	count(i) = count(i) + 1
	total = total + 1
	run = run - npe
      enddo
    enddo

  CASE DEFAULT
    return                                  ! ERROR, invalid mode
  END SELECT

  if(offset(npe) + count(npe) .ne. npts) return  ! ERROR, count mismatch

! for tiles with ZERO length at end, lmin would be gmax + 1, lmax would be gmax
  lmin = gmin + offset(myrank + 1)          ! myrank is in origin 0, offset indexing is in origin 1
  lmax = lmin + count(myrank + 1) -1        ! count indexing is in origin 1
  lmax = min(gmax, lmax)                    ! this should NEVER happen
  return
 end subroutine RPN_data_distribute                               !InTf!
