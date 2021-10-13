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

! Mpich Fortran semble appeler la routine C directement, contrairement a OpenMpi
! teste avec Cray , Intel OneAPI, gcc+mpich
#if ! defined(USE_MPICH)
subroutine MPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER(KIND=C_INTPTR_T) SIZE, BASEPTR
  INTEGER :: DISP_UNIT, INFO, COMM, WIN, IERROR
  if(DEBUG > 0) print *,'MPI_WIN_ALLOCATE'
  call PMPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)
  call MPI_Track_Window(WIN, 1)  ! insert
end subroutine MPI_WIN_ALLOCATE

subroutine MPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER(KIND=C_INTPTR_T) SIZE, BASEPTR
  INTEGER DISP_UNIT, INFO, COMM, WIN, IERROR
  if(DEBUG > 0) print *,'MPI_WIN_ALLOCATE_SHARED'
  call PMPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)
  call MPI_Track_Window(WIN, 1)  ! insert
end subroutine MPI_WIN_ALLOCATE_SHARED

subroutine MPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER BASE(*)
  INTEGER(KIND=C_INTPTR_T) SIZE
  INTEGER DISP_UNIT, INFO, COMM, WIN, IERROR
  if(DEBUG > 0) print *,'MPI_WIN_CREATE'
  call PMPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)
  call MPI_Track_Window(WIN, 1)  ! insert
end subroutine MPI_WIN_CREATE

subroutine MPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, IERROR)   ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  INTEGER INFO, COMM, WIN, IERROR
  call PMPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, IERROR)
  call MPI_Track_Window(WIN, 1)  ! insert
end subroutine MPI_WIN_CREATE_DYNAMIC

subroutine MPI_WIN_FREE(WIN, IERROR)
  use :: MPI_tracked_windows_mod
  INTEGER :: WIN, IERROR
  if(DEBUG > 0) print *,'MPI_WIN_FREE',win
  call MPI_Track_Window(WIN, 0)  ! remove from tables
  call PMPI_WIN_FREE(WIN, IERROR)
end subroutine MPI_WIN_FREE
#endif

#if defined(MPIF08)

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
  if(DEBUG > 0) print *,'MPI_WIN_ALLOCATE_f08'
  call PMPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, 1)  ! insert
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
  if(DEBUG > 0) print *,'MPI_WIN_ALLOCATE_SHARED_f08'
  call PMPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, 1)  ! insert
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
  if(DEBUG > 0) print *,'MPI_WIN_CREATE_f08'
  call PMPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_CREATE_f08

subroutine MPI_WIN_CREATE_DYNAMIC_f08(INFO, COMM, WIN, IERROR)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_CREATE_DYNAMIC, MPI_Info, MPI_Comm, MPI_Win, MPI_ADDRESS_KIND
  TYPE(MPI_Info), INTENT(IN) :: info
  TYPE(MPI_Comm), INTENT(IN) :: comm
  TYPE(MPI_Win), INTENT(OUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  if(DEBUG > 0) print *,'MPI_WIN_CREATE_DYNAMIC_f08'
  call PMPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, ierr)
  call MPI_Track_Window(win%MPI_VAL, 1)  ! insert
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_CREATE_DYNAMIC_f08

subroutine MPI_WIN_FREE_f08(WIN, IERROR)     ! use mpi_f08 OpenMPI + mpich
  use :: MPI_tracked_windows_mod
  USE :: mpi_f08, ONLY : PMPI_WIN_FREE, MPI_Win
  TYPE(MPI_Win), INTENT(INOUT) :: win
  INTEGER, OPTIONAL, INTENT(OUT) :: IERROR
  integer :: ierr
  if(DEBUG > 0) print *,'MPI_WIN_FREE_f08',win
  call MPI_Track_Window(win%MPI_VAL, 0)  ! remove from tables
  call PMPI_WIN_FREE(WIN, ierr)
  if(present(IERROR)) IERROR = ierr
end subroutine MPI_WIN_FREE_f08

#endif
