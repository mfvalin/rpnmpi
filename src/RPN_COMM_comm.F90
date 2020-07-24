
!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
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
module rpncomm_com
  use rpn_comm
  implicit none
  integer, parameter :: MAX_COMM_TAB=128
  type(symtab), dimension(:), pointer, save :: com_tab => NULL()  ! communicator translation table
  integer, save :: defcom_index=-1                                ! index for RPN_COMM_DEFAULT
  integer, save :: max_com_index=0                    
 contains
!
! allocate and initialize com_tab, the communicator table
!
  subroutine init_com_tab
    implicit none
    integer :: i

    if( associated(com_tab) ) return  ! job already done

    allocate(com_tab(MAX_COMM_TAB+1))   ! allocate table and initialize "the usual suspects"
    com_tab( 1) = symtab(pe_grid_peers,RPN_COMM_GRIDPEERS)
    com_tab( 2) = symtab(pe_indomm,RPN_COMM_GRID)
    com_tab( 3) = symtab(pe_indomm,RPN_COMM_DOMM)
    com_tab( 4) = symtab(WORLD_COMM_MPI,RPN_COMM_WORLD)
    com_tab( 5) = symtab(WORLD_COMM_MPI,RPN_COMM_ALLDOMAINS)
    com_tab( 6) = symtab(pe_a_domain,RPN_COMM_ALLGRIDS)
    com_tab( 7) = symtab(pe_indomms,RPN_COMM_MULTIGRID)
    com_tab( 8) = symtab(pe_wcomm,RPN_COMM_ALL)
    com_tab( 9) = symtab(pe_defcomm,RPN_COMM_DEFAULT)  ; defcom_index = 9  ! this item might get updated by rpn_comm_defo
    com_tab(10) = symtab(pe_myrow,RPN_COMM_EW)
    com_tab(11) = symtab(pe_mycol,RPN_COMM_NS)
    com_tab(12) = symtab(pe_blocmaster,RPN_COMM_BLOCMASTER)
    com_tab(13) = symtab(pe_bloc,RPN_COMM_BLOCK)
    com_tab(14) = symtab(MPI_COMM_WORLD,RPN_COMM_UNIVERSE)
    com_tab(15) = symtab(MPI_COMM_NULL,RPN_COMM_NULL)
    max_com_index = 15
    do i = 16,MAX_COMM_TAB+1
      com_tab(i) = symtab(MPI_COMM_NULL,"")
    enddo
!print *,'DEBUG: com_tab initialized'
!print *,com_tab(1:15)
  end subroutine init_com_tab
!
! add new communicator to table
!
  function new_com(com,mpicom) result(indx)
    implicit none
    character(len=*), intent(IN) :: com      ! rpn_comm character style communicator
    integer, intent(IN) :: mpicom            ! MPI integer communicator
    integer :: indx

    if(.not. associated(com_tab)) call init_com_tab
    if(max_com_index < MAX_COMM_TAB) then
      max_com_index = max_com_index + 1
      com_tab(max_com_index)%number = mpicom
      com_tab(max_com_index)%string = trim(com)
      indx = max_com_index
    else
      indx = -1
    endif
    return
  end function new_com
!
! get index of string com in com_tab
!
  function indx_com_tab(com) result(indx)
    implicit none
    character(len=*), intent(IN) :: com
    integer :: indx
    integer :: i

    if(.not. associated(com_tab)) call init_com_tab
    indx = MAX_COMM_TAB+1  ! will be returned if com is not found
    do i = 1,max_com_index
      if( trim(com_tab(i)%string) == trim(com) ) then
        indx = i
!print *,'DEBUG: index of "'//trim(com)//'" is',indx
        return
      endif
    enddo
  end function indx_com_tab
end module rpncomm_com
!InTf!
      integer function RPN_COMM_comm(com)                    !InTf!
!	Luc Corbeil, 2000-11-21
!
!	lien entre chaine de caractere de communicateur
!	GRID, EW et NS et leur numero assigne par
!	MPI.
!
      use rpncomm_com
      implicit none                                 !InTf!
      character(len=*), intent(IN) :: com           !InTf!
      character(len=32) comm
      integer :: indx

      if(.not. associated(com_tab)) call init_com_tab
      call rpn_comm_low2up(com,comm)
!print *,'DEBUG: avant indx_com_tab'
      indx = indx_com_tab(comm)
!print *,'DEBUG: indx_com_tab=',indx
      RPN_COMM_comm = com_tab(indx_com_tab(comm))%number
!print *,'DEBUG: number=',RPN_COMM_comm

      return
      end function RPN_COMM_comm                                  !InTf!
