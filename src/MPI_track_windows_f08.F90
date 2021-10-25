! functions for FORTRAN MPI
! Copyright (C) 2021  Environnement et Changement climatique Canada
!
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! Author:
!     M. Valin,   Recherche en Prevision Numerique, 2021

#if ! defined(DEBUG)
#define DEBUG 0
#endif

! interceptors for the Fortran 2008 interface ( use mpi_f08 )

subroutine MPI_Win_allocate_f08(size, disp_unit, info, comm, baseptr, win, ierror)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_ALLOCATE, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
  INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
  INTEGER, INTENT(IN) :: disp_unit
  TYPE(MPI_Info), INTENT(IN) :: info
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(C_PTR), INTENT(OUT) :: baseptr
  TYPE(MPI_Win), INTENT(OUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_ALLOCATE_f08'
  call PMPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, my_rank + 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_Win_allocate_f08

subroutine MPI_WIN_ALLOCATE_SHARED_f08(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_Win_allocate_shared, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
  INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
  INTEGER, INTENT(IN) :: disp_unit
  TYPE(MPI_Info), INTENT(IN) :: info
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(C_PTR), INTENT(OUT) :: baseptr
  TYPE(MPI_Win), INTENT(OUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_ALLOCATE_SHARED_f08'
  call PMPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, my_rank + 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_ALLOCATE_SHARED_f08

subroutine MPI_WIN_CREATE_f08(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)     ! use mpi_f08 OpenMPI
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_CREATE, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
!     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: base
  integer, dimension(*), ASYNCHRONOUS :: base
  INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
  INTEGER, INTENT(IN) :: disp_unit
  TYPE(MPI_Info), INTENT(IN) :: info
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(MPI_Win), INTENT(OUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_CREATE_f08'
  call PMPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, my_rank + 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_CREATE_f08

#if defined(USE_MPICH)
! subroutine MPI_WIN_CREATE_f08ts(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)     ! use mpi_f08 Mpich
!   use :: MPI_tracked_windows_mod
!   USE :: mpi_f08, ONLY : PMPI_WIN_CREATE, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
! !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: base
!   integer, dimension(*), ASYNCHRONOUS :: base
!   INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
!   INTEGER, INTENT(IN) :: disp_unit
!   TYPE(MPI_Info), INTENT(IN) :: info
!   TYPE(MPI_Comm), INTENT(IN) :: comm
!   TYPE(MPI_Win), INTENT(OUT) :: win
!   INTEGER, OPTIONAL, INTENT(OUT) :: ierror
!   integer :: ierr
!   if(DEBUG > 0) print *,'MPI_WIN_CREATE_f08ts'
!   call PMPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, ierr)  ! seems to call the C function MPI_Win_create
! !   call MPI_Track_Window(win%MPI_VAL, 1)  ! insert
!   if(present(IERROR)) IERROR = ierr
! end subroutine MPI_WIN_CREATE_f08ts
#endif

subroutine MPI_WIN_CREATE_DYNAMIC_f08(INFO, COMM, WIN, IERROR)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_CREATE_DYNAMIC, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
  TYPE(MPI_Info), INTENT(IN) :: info
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(MPI_Win), INTENT(OUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_CREATE_DYNAMIC_f08'
  call PMPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, my_rank + 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_CREATE_DYNAMIC_f08

subroutine MPI_WIN_FREE_f08(WIN, IERROR)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_FREE, MPI_Win
  TYPE(MPI_Win), INTENT(INOUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: IERROR
  integer :: ierr
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_FREE_f08',win
  call MPI_Track_Window(win%MPI_VAL, 0)  ! remove from tables
  call PMPI_WIN_FREE(WIN, ierr)
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_FREE_f08

#if 0
subroutine MPI_Sendrecv_f08(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
  USE :: mpi_f08, ONLY : MPI_Datatype, MPI_Comm, MPI_Status
!     TYPE(*), DIMENSION(..), INTENT(IN) :: sendbuf
  integer, dimension(*), INTENT(IN) :: sendbuf
!     TYPE(*), DIMENSION(..) :: recvbuf
  integer, dimension(*) :: recvbuf
  INTEGER, INTENT(IN) :: sendcount, dest, sendtag, recvcount, source, recvtag
  TYPE(MPI_Datatype), INTENT(IN) :: sendtype, recvtype
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(MPI_Status) :: status
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  if(DEBUG > 0) print *,'DEBUG: MPI_Sendrecv_f08'
  call PMPI_Sendrecv(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, ierr)
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_Sendrecv_f08
#endif
