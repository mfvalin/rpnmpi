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
subroutine RPN_MPI_get_mpi_definitions_raw(what, ierr)          !InTf! ! get a copy of MPI definitions
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI_mpi_symbols.inc'
!! import :: RPN_MPI_mpi_definitions_raw                        !InTf! 
  type(RPN_MPI_mpi_definitions_raw), intent(INOUT) :: what      !InTf! 
  integer, intent(OUT) :: ierr                                  !InTf! 
  call RPN_MPI_get_mpi_definitions(what, ierr)
  return
end subroutine RPN_MPI_get_mpi_definitions_raw          !InTf! ! get a copy of MPI definitions

subroutine RPN_MPI_get_mpi_definitions(what, ierr)          !InTf! ! get a copy of MPI definitions
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_MPI_mpi_symbols.inc'
! white lie in published interface, what is treated as the "raw" type rather than the "wrapped" type
!! import :: RPN_MPI_mpi_definitions                        !InTf! 
!!  type(RPN_MPI_mpi_definitions), intent(INOUT) :: what    !InTf! 
  type(RPN_MPI_mpi_definitions_raw), intent(INOUT) :: what
  integer, intent(OUT) :: ierr                              !InTf! 

  ierr = MPI_ERROR
  if(what%version .ne. mpi_symbols_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_definitions) version mismatch, expected :',mpi_symbols_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS
  what%MPI_GROUP_NULL = MPI_GROUP_NULL
  what%MPI_COMM_NULL = MPI_COMM_NULL
  what%MPI_DATATYPE_NULL = MPI_DATATYPE_NULL
  what%MPI_REQUEST_NULL = MPI_REQUEST_NULL
  what%MPI_OP_NULL = MPI_OP_NULL
  what%MPI_ERRHANDLER_NULL = MPI_ERRHANDLER_NULL
  what%MPI_INFO_NULL = MPI_INFO_NULL
  what%MPI_WIN_NULL = MPI_WIN_NULL
  what%MPI_STATUS_SIZE = MPI_STATUS_SIZE
  what%MPI_ANY_SOURCE = MPI_ANY_SOURCE
  what%MPI_ANY_TAG = MPI_ANY_TAG
  what%MPI_SUCCESS = MPI_SUCCESS
  what%MPI_ERROR = MPI_ERROR
  what%MPI_COMM_WORLD = MPI_COMM_WORLD
  what%MPI_COMM_SELF = MPI_COMM_SELF
  what%MPI_GROUP_EMPTY = MPI_GROUP_EMPTY
  what%MPI_COMM_TYPE_SHARED = MPI_COMM_TYPE_SHARED
  what%MPI_ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL
  what%MPI_ERRORS_RETURN = MPI_ERRORS_RETURN
  what%MPI_BYTE = MPI_BYTE
  what%MPI_PACKED = MPI_PACKED
  what%MPI_UB = MPI_UB
  what%MPI_LB = MPI_LB
  what%MPI_CHARACTER = MPI_CHARACTER
  what%MPI_LOGICAL = MPI_LOGICAL
  what%MPI_INTEGER = MPI_INTEGER
  what%MPI_INTEGER1 = MPI_INTEGER1
  what%MPI_INTEGER2 = MPI_INTEGER2
  what%MPI_INTEGER4 = MPI_INTEGER4
  what%MPI_INTEGER8 = MPI_INTEGER8
  what%MPI_INTEGER16 = MPI_INTEGER16
  what%MPI_REAL = MPI_REAL
!   what%MPI_REAL2 = MPI_REAL2
  what%MPI_REAL4 = MPI_REAL4
  what%MPI_REAL8 = MPI_REAL8
  what%MPI_REAL16 = MPI_REAL16
  what%MPI_DOUBLE_PRECISION = MPI_DOUBLE_PRECISION
  what%MPI_COMPLEX = MPI_COMPLEX
  what%MPI_COMPLEX8 = MPI_COMPLEX8
  what%MPI_COMPLEX16 = MPI_COMPLEX16
  what%MPI_COMPLEX32 = MPI_COMPLEX32
  what%MPI_DOUBLE_COMPLEX = MPI_DOUBLE_COMPLEX
  what%MPI_2REAL = MPI_2REAL
  what%MPI_2DOUBLE_PRECISION = MPI_2DOUBLE_PRECISION
  what%MPI_2INTEGER = MPI_2INTEGER
!   what%MPI_2COMPLEX = MPI_2COMPLEX
!   what%MPI_2DOUBLE_COMPLEX = MPI_2DOUBLE_COMPLEX
!   what%MPI_LOGICAL1 = MPI_LOGICAL1
!   what%MPI_LOGICAL2 = MPI_LOGICAL2
!   what%MPI_LOGICAL4 = MPI_LOGICAL4
!   what%MPI_LOGICAL8 = MPI_LOGICAL8
  what%MPI_MAX = MPI_MAX
  what%MPI_MIN = MPI_MIN
  what%MPI_SUM = MPI_SUM
  what%MPI_PROD = MPI_PROD
  what%MPI_LAND = MPI_LAND
  what%MPI_BAND = MPI_BAND
  what%MPI_LOR = MPI_LOR
  what%MPI_BOR = MPI_BOR
  what%MPI_LXOR = MPI_LXOR
  what%MPI_BXOR = MPI_BXOR
  what%MPI_MAXLOC = MPI_MAXLOC
  what%MPI_MINLOC = MPI_MINLOC
  what%MPI_REPLACE = MPI_REPLACE
  what%MPI_THREAD_SINGLE = MPI_THREAD_SINGLE
  what%MPI_THREAD_FUNNELED = MPI_THREAD_FUNNELED
  what%MPI_THREAD_SERIALIZED = MPI_THREAD_SERIALIZED
  what%MPI_THREAD_MULTIPLE = MPI_THREAD_MULTIPLE

  return
end subroutine RPN_MPI_get_mpi_definitions                  !InTf! 

module RPN_MPI_mpi_layout
  use ISO_C_BINDING
  implicit none
! cannot include RPN_MPI_mpi_definitions.inc because it contains interfaces 
! for RPN_MPI_get_mpi_definitions and RPN_MPI_get_mpi_layout
  include 'RPN_MPI_mpi_symbols.inc'
  include 'RPN_MPI_mpi_layout.inc'
  type(mpi_layout_internal), save :: ml
end module RPN_MPI_mpi_layout

! initialize RPN_MPI internal mpi layout (communicators, ranks, sizes)
subroutine RPN_MPI_reset_mpi_layout()                       !InTf! 
  use RPN_MPI_mpi_layout
  implicit none
  include 'mpif.h'
  include 'RPN_MPI_system_interfaces.inc'
  integer :: cpu

  ml%host = get_host_id()              ! get linux host id
  cpu = sched_get_my_cpu()             ! get logical cpu number
  ml%numa = numa_node(cpu)             ! get numa space number for this cpu
  ml%colors = [-1, -1, -1]

  ml%comm%wrld = application(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL)
  ml%rank%wrld = application(-1, -1, -1)
  ml%size%wrld = application(-1, -1, -1)

  ml%comm%appl = ml%comm%wrld
  ml%rank%appl = ml%rank%wrld
  ml%size%appl = ml%size%wrld

  ml%comm%sgrd = mpigrid(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
			 MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
			 MPI_COMM_NULL, MPI_COMM_NULL)
  ml%rank%sgrd = mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
  ml%size%sgrd = mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

  ml%comm%grid = ml%comm%sgrd
  ml%rank%grid = ml%rank%sgrd
  ml%size%grid = ml%size%sgrd

  ml%comm%blck = subgrid(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL)
  ml%rank%blck = subgrid(-1, -1, -1)
  ml%size%blck = subgrid(-1, -1, -1)

  return
end subroutine RPN_MPI_reset_mpi_layout                     !InTf! 

! same as RPN_MPI_get_mpi_layout but advertized with type mpi_layout_internal rather than mpi_layout
! get a copy on RPN_MPI internal mpi layout (communicators, ranks, sizes)
subroutine RPN_MPI_get_mpi_layout_raw(what, ierr)            !InTf! 
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
!! import :: mpi_layout_internal, C_INT                      !InTf! 
  type(mpi_layout_internal), intent(INOUT) :: what           !InTf! 
  integer(C_INT), intent(OUT) :: ierr                        !InTf! 
  call RPN_MPI_get_mpi_layout(what, ierr)
  return
end subroutine RPN_MPI_get_mpi_layout_raw                    !InTf! 

! get a copy on RPN_MPI internal mpi layout (communicators, ranks, sizes)
subroutine RPN_MPI_get_mpi_layout(what, ierr)               !InTf! 
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
  include 'mpif.h'
! white lie in published interface, what is treated as the "internal" type rather than the "wrapped" type
!! import :: mpi_layout, C_INT                              !InTf! 
!! type(mpi_layout), intent(INOUT) :: what                  !InTf! 
  type(mpi_layout_internal), intent(INOUT) :: what
  integer(C_INT), intent(OUT) :: ierr                       !InTf! 
  integer :: i

  ierr = MPI_ERROR
  if(what%version .ne. layout_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_layout) version mismatch, expected :',layout_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS

  what = ml
  what%comm%sgrd%row       = MPI_COMM_NULL   ! row is not defined for supergrids
  what%rank%sgrd%row       = -1
  what%size%sgrd%row       = -1
  what%comm%sgrd%column    = MPI_COMM_NULL   ! and neither is column
  what%rank%sgrd%column    = -1
  what%size%sgrd%column    = -1
  return
end subroutine RPN_MPI_get_mpi_layout                       !InTf! 
