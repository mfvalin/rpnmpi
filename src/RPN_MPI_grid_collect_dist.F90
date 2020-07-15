!/! RPN_MPI - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2020  Recherche en Prevision Numerique
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
! !/
#include <RPN_MPI_macros.hf>
 subroutine RPN_MPI_grid_collect(topo, &                                               !InTf!
                            zloc, mini, maxi, minj, maxj, &                            !InTf!
                            zglb, gni, gnj, ierr) bind(C,name='RPN_MPI_grid_collect')  !InTf!
  use rpn_mpi_mpif
  implicit none
  include 'RPN_MPI_mpi_symbols.inc'
  include 'RPN_MPI_mpi_layout.inc'
!! import :: RPN_MPI_Loc, RPN_MPI_Ftopo, C_INT                   !InTf!
  type(RPN_MPI_Ftopo), intent(IN) :: topo                        !InTf!
  integer(C_INT), intent(IN) :: mini,maxi,minj,maxj,gni,gnj      !InTf!
  integer, intent(OUT) :: ierr                                   !InTf!
! white lie in published interface, zloc and zglb are published as addresses passed by value
!! type(RPN_MPI_Loc), intent(IN), value :: zloc, zglb            !InTf!
  integer, intent(IN),  dimension(mini:maxi,minj:maxj) :: zloc
  integer, intent(OUT), dimension(gni,gnj)             :: zglb

  integer :: npex, npey, lni, lnj, nw, i, i0, j, np, lnimax, lnjmax
  integer, dimension(:,:,:), allocatable :: t1
  integer, dimension(:,:),   allocatable :: t2
  integer, dimension(:),     allocatable :: counts, displs
  integer, dimension(:), pointer :: lis, ljs

  ierr = MPI_ERROR
  if(VaL(topo%grd%comm) == MPI_COMM_NULL .or. &
     VaL(topo%row%comm) == MPI_COMM_NULL .or. &
     VaL(topo%col%comm) == MPI_COMM_NULL) return

  npex = topo%row%size
  npey = topo%col%size
  lnimax  = (gni + npex - 1) / npex
  lnjmax  = (gnj + npey - 1) / npey
  lnj = topo%lj
  if(topo%row%rank == 0) then         ! only needed on column 0
    allocate(t1(mini:maxi,lnj,npex))
  else
    allocate(t1(1,1,1))
  endif
  nw = (maxi - mini + 1) * lnj

! gather along row
  call MPI_Gather(zloc, nw, MPI_INTEGER, t1, nw,  MPI_INTEGER, 0, topo%row%comm, ierr)
  if(topo%row%rank .ne. 0) goto 1  ! not on column 0, job done

! reshape and eliminate row halo
  allocate(t2(gni,lnj), counts(npey), displs(npey))
  call C_F_POINTER(topo%lis, lis, [npex])
  i0 = 0
  do np = 1, npex
    lni = lis(np)   ! length along i for tile np
    do j = 1 , lnj
      t2(i0+1:i0+lni,j) = t1(1:lni,j,np)
    enddo
    i0 = i0 + lni
  enddo

! now fill zglb array from t2
  call C_F_POINTER(topo%ljs, ljs, [npey])
  counts(1:npey) = ljs(1:npey) * gni
  displs(1) = 0
  do i = 2 , npey
    displs(i) = displs(i-1) + counts(i)
  enddo
  call MPI_Gatherv(t2, gni*lnj, MPI_INTEGER, zglb, counts,  displs, MPI_INTEGER, 0, topo%row%comm, ierr)
  deallocate(t2, counts, displs)         ! no longer needed

1 continue
  deallocate(t1)         ! no longer needed
  ierr = MPI_SUCCESS
 end subroutine RPN_MPI_grid_collect                             !InTf!

 subroutine RPN_MPI_grid_dist(topo, &                                               !InTf!
                            zloc, mini, maxi, minj, maxj, &                         !InTf!
                            zglb, gni, gnj, ierr) bind(C,name='RPN_MPI_grid_dist')  !InTf!
  use rpn_mpi_mpif
  implicit none
  include 'RPN_MPI_mpi_symbols.inc'
  include 'RPN_MPI_mpi_layout.inc'
!! import :: RPN_MPI_Loc, RPN_MPI_Ftopo, C_INT                   !InTf!
  type(RPN_MPI_Ftopo), intent(IN) :: topo                        !InTf!
  integer(C_INT), intent(IN) :: mini,maxi,minj,maxj,gni,gnj      !InTf!
  integer, intent(OUT) :: ierr                                   !InTf!
! white lie in published interface, zloc and zglb are published as addresses passed by value
!! type(RPN_MPI_Loc), intent(IN), value :: zloc, zglb            !InTf!
  integer, intent(OUT), dimension(mini:maxi,minj:maxj) :: zloc
  integer, intent(IN), dimension(gni,gnj)              :: zglb

  ierr = MPI_ERROR
  if(VaL(topo%grd%comm) == MPI_COMM_NULL .or. &
     VaL(topo%row%comm) == MPI_COMM_NULL .or. &
     VaL(topo%col%comm) == MPI_COMM_NULL) return

  ierr = MPI_SUCCESS
 end subroutine RPN_MPI_grid_dist                                !InTf!
