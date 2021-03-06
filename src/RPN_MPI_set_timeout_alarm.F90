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
      function RPN_MPI_set_timeout_alarm(seconds) result(seconds_since)  !InTf!
      use ISO_C_BINDING
      implicit none
      integer, intent(IN) :: seconds  !InTf!
      integer :: seconds_since  !InTf!
#include <RPN_MPI_system_interfaces.hf>

      seconds_since = c_alarm(seconds)
      return
      end function RPN_MPI_set_timeout_alarm                             !InTf!
!     legacy function for rpn_comm
      function RPN_COMM_set_timeout_alarm(seconds) result(seconds_since)  !InTf!
      use ISO_C_BINDING
      implicit none
      integer, intent(IN) :: seconds  !InTf!
      integer :: seconds_since  !InTf!
#include <RPN_MPI_system_interfaces.hf>

      seconds_since = c_alarm(seconds)
      return
      end function RPN_COMM_set_timeout_alarm                             !InTf!
