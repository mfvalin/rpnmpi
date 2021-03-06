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
      integer function RPN_MPI_get_a_free_unit()                   !InTf!
      implicit none
      integer :: i
      character (len=16) :: access_mode
      RPN_MPI_get_a_free_unit=-1
      do i = 99,1,-1  ! find an available unit number
        inquire(UNIT=i,ACCESS=access_mode)
        if(trim(access_mode) == 'UNDEFINED')then ! found
          RPN_MPI_get_a_free_unit = i
          exit
        endif
      enddo
      return
      end function RPN_MPI_get_a_free_unit                         !InTf!
!     legacy routine for rpn_comm routines
      integer function RPN_COMM_get_a_free_unit()                   !InTf!
      implicit none
      integer, external :: RPN_MPI_get_a_free_unit
      RPN_COMM_get_a_free_unit = RPN_MPI_get_a_free_unit()
      return
      end function RPN_COMM_get_a_free_unit                         !InTf!
