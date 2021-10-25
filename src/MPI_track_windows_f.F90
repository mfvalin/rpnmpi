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

! interceptors for the traditional Fortran interface ( use mpi or include 'mpif.h' )
! apparently not needed (nor wanted) when using mpich

#if ! defined(USE_MPICH)
subroutine MPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER(KIND=C_INTPTR_T) SIZE, BASEPTR
  INTEGER :: DISP_UNIT, INFO, COMM, WIN, IERROR
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_ALLOCATE Fortran'
  call PMPI_WIN_ALLOCATE(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)
  call MPI_Track_Window(WIN, my_rank + 1)  ! insert
end subroutine MPI_WIN_ALLOCATE

subroutine MPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER(KIND=C_INTPTR_T) SIZE, BASEPTR
  INTEGER DISP_UNIT, INFO, COMM, WIN, IERROR
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_ALLOCATE_SHARED Fortran'
  call PMPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)
  call MPI_Track_Window(WIN, my_rank + 1)  ! insert
end subroutine MPI_WIN_ALLOCATE_SHARED

subroutine MPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)    ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
  INTEGER BASE(*)
  INTEGER(KIND=C_INTPTR_T) SIZE
  INTEGER DISP_UNIT, INFO, COMM, WIN, IERROR
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_CREATE Fortran'
  call PMPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)
  call MPI_Track_Window(WIN, my_rank + 1)  ! insert
end subroutine MPI_WIN_CREATE

subroutine MPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, IERROR)   ! include 'mpif.h'  |  use mpi
  use :: MPI_tracked_windows_mod
  INTEGER INFO, COMM, WIN, IERROR
  integer :: my_rank
  call PMPI_COMM_RANK(COMM, my_rank, IERROR)
  my_rank = my_rank + my_rank
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_CREATE_DYNAMIC Fortran'
  call PMPI_WIN_CREATE_DYNAMIC(INFO, COMM, WIN, IERROR)
  call MPI_Track_Window(WIN, my_rank + 1)  ! insert
end subroutine MPI_WIN_CREATE_DYNAMIC

subroutine MPI_WIN_FREE(WIN, IERROR)
  use :: MPI_tracked_windows_mod
  INTEGER :: WIN, IERROR
  if(DEBUG > 0) print *,'DEBUG: MPI_WIN_FREE Fortran',win
  call MPI_Track_Window(WIN, 0)  ! remove from tables
  call PMPI_WIN_FREE(WIN, IERROR)
end subroutine MPI_WIN_FREE

! subroutine MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
!     INTEGER    BUF(*)
!     INTEGER    COUNT, DATATYPE, DEST, TAG, COMM, IERROR
!     if(DEBUG > 0) print *,'DEBUG: MPI_SEND Fortran'
!     call PMPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
! end subroutine MPI_SEND

subroutine MPI_SENDRECV(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
  use :: mpi, ONLY : MPI_STATUS_SIZE
  INTEGER    SENDBUF(*), RECVBUF(*)
  INTEGER    SENDCOUNT, SENDTYPE, DEST, SENDTAG
  INTEGER    RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM
  INTEGER    STATUS(MPI_STATUS_SIZE), IERROR
  if(DEBUG > 0) print *,'DEBUG: MPI_SENDRECV Fortran'
    call PMPI_SENDRECV(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
end subroutine MPI_SENDRECV
#endif
