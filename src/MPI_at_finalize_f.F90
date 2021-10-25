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

! interfaces to C functions for 1 sided window tracking
module MPI_tracked_windows_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
  interface

    subroutine MPI_Track_Window(win, insert) bind(C,name='MPI_Track_Window_f')     ! window tracker with Fortran MPI window
      import :: C_INT
      integer(C_INT), intent(IN), value :: win                                     ! Fortran 1 sided MPI window descriptor
      integer(C_INT), intent(IN), value :: insert                                  ! 1 if created window, 0 if window to free
    end subroutine

    subroutine MPI_Free_Tracked_Windows() bind(C,name='MPI_Free_Tracked_Windows')  ! frees all remaining not freed MPI windows
    end subroutine

  end interface
end module MPI_tracked_windows_mod

! interfaces to C functions for at finalize registration
module MPI_At_Finalize_mod
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_FUNPTR
  interface

    subroutine MPI_At_Finalize_exec() bind(C,name='MPI_At_Finalize_exec')   ! execute functions registered with MPI_At_Finalize
    end subroutine

    subroutine MPI_At_Finalize_c(fn) bind(C,name='MPI_At_Finalize')         ! register a function for execution before MPI_Finalize
      import :: C_FUNPTR
      type(C_FUNPTR), intent(IN), value :: fn                               ! pointer to function  ( C_FUNLOC(external_function) )
    end subroutine

  end interface
end module MPI_At_Finalize_mod

! intercept user call to Fortran MPI_init  ( include 'mpif.h'  |  use mpi )
#if ! defined(USE_MPICH)
subroutine MPI_init(ierror)
  use :: MPI_tracked_windows_mod
  use :: MPI_At_Finalize_mod
  integer :: ierror
  if(DEBUG > 0) print *,'DEBUG: MPI_init Fortran'
  call PMPI_init(ierror)
#if ! defined(USE_MPICH)
  call MPI_At_Finalize(MPI_Free_Tracked_Windows)  ! will be registered by C MPI_Init with Mpich
#endif
end
#endif

! intercept user call to Fortran MPI_finalize  ( include 'mpif.h'  |  use mpi )
#if ! defined(USE_MPICH)
subroutine MPI_finalize(ierror)
  use :: MPI_tracked_windows_mod
  use :: MPI_At_Finalize_mod
  integer :: ierror
  if(DEBUG > 0) print *,'DEBUG: MPI_finalize Fortran'
#if ! defined(USE_MPICH)
  call MPI_At_Finalize_exec()                        ! redundant for Mpich, necessary for OpenMPI
#endif
  call PMPI_finalize(ierror)
end
#endif

! register a function for execution before MPI_Finalize
subroutine MPI_At_Finalize(fn)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_FUNLOC
  use MPI_tracked_windows_mod
  use :: MPI_At_Finalize_mod
  external :: fn                                     ! Fortran function to register
  call MPI_At_Finalize_c(C_FUNLOC(fn))               ! call C function
end subroutine

! things are slightly different with Fortran 2008 interface
#if defined(MPIF08)

! intercept user call to Fortran MPI_init
subroutine MPI_init_f08(ierror)          ! use mpi_f08 OpenMPI + mpich
  USE :: mpi_f08, ONLY : PMPI_init
  use :: MPI_tracked_windows_mod
  use :: MPI_At_Finalize_mod
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  if(DEBUG > 0) print *,'DEBUG: MPI_init_f08 Fortran 2008'
  call PMPI_init(ierr)
  if(present(IERROR)) IERROR = ierr
  call MPI_At_Finalize(MPI_Free_Tracked_Windows)
return
end

! intercept user call to Fortran MPI_finalize
subroutine MPI_finalize_f08(ierror)      ! use mpi_f08 OpenMPI + mpich
  USE :: mpi_f08, ONLY : PMPI_finalize
  use :: MPI_tracked_windows_mod
  use :: MPI_At_Finalize_mod
  INTEGER, OPTIONAL, INTENT(OUT) :: ierror
  integer :: ierr
  if(DEBUG > 0) print *,'DEBUG: MPI_finalize_f08 Fortran 2008'
  call MPI_At_Finalize_exec()                       ! get functions registered with MPI_At_Finalize executed
  call PMPI_finalize(ierr)
  if(present(IERROR)) IERROR = ierr
return
end

#endif
