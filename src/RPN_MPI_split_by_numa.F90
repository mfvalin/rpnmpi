! RPN_MPI - Library of useful routines for C and FORTRAN programming
! Copyright (C) 2020  Division de Recherche en Prevision Numerique
!                     Environnement Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
#define IN_RPN_MPI_split_by
!****m* rpn_mpi/NumaSplit
! DESCRIPTION
!  manage communicators across NUMA spaces and/or SMP nodes
!  - intra node communicators (1 SMP node is likely to include >1 NUMA domains)
!  - communicators between Pes in different SMP nodes but with same intra-node communicator rank
!  - intra NUMA space communicators
!  - communicators between Pes in different NUMA domains with same rank in NUMA space communicators
!
! SEE ALSO
!   RPN_MPI_split_by_memory_domain
!
! AUTHOR
!  M.Valin, Recherche En Prevision Numerique 2020
! IGNORE
!
module RPN_MPI_memspace
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
  save
  integer, parameter :: MAX_CACHE=16
! VARIABLES
! internally kept variables
! cache for communicator split
  integer, dimension(MAX_CACHE) :: corig     ! original communicator
  integer, dimension(MAX_CACHE) :: cnode     ! communicator for PEs on same node
  integer, dimension(MAX_CACHE) :: cnuma     ! communicator for PEs in same NUMA (subset of cnode)
  integer, dimension(MAX_CACHE) :: cpeer     ! communicator for "same rank on node" peers across corig
  integer, dimension(MAX_CACHE) :: speer     ! communicator for "same rank in NUMA" peers across corig
  integer       :: ncached = 0               ! number of cached entries
  integer       :: my_numa = -1              ! NUMA domain for this PE
!******

end module RPN_MPI_memspace
!****P* rpn_mpi/SMP_Commnicators
! DESCRIPTION
!   intra SMP node communicator management package (split by node or by NUMA space)
!******
!InTf!
!****f* rpn_mpi/RPN_MPI_split_by_memory_domain
! DESCRIPTION
!   split a communicator on a node(host), socket, or NUMA space basis
!   THIS version treata split by socket as split by numa
! SYNOPSIS
subroutine RPN_MPI_split_by_memory_domain(   &       !InTf!
           origcomm, split_by,               &       !InTf!
           nodecomm, numacomm, peercomm,     &       !InTf!
           noderank, numarank, peerrank, isiz, ierr) !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! IGNORE
  use RPN_MPI_memspace
  implicit none
#include <RPN_MPI_system_interfaces.hf>
#include <RPN_MPI_mpi_symbols.hf>
! ARGUMENTS
  integer, intent(IN)  :: origcomm  ! MPI communicator to split on a host basis        !InTf!
  integer, intent(IN)  :: split_by  ! type of split desired                            !InTf!
!                         RPN_MPI_BY_NUMA, RPN_MPI_BY_SOCKET, RPN_MPI_BY_HOST, RPN_MPI_BY_NODE
  integer, intent(OUT) :: nodecomm  ! new communicator to be used py PEs on same host  !InTf!
  integer, intent(OUT) :: numacomm  ! new communicator to be used py PEs in same NUMA  !InTf!
  integer, intent(OUT) :: peercomm  ! communicator for node/numa peers                 !InTf!
  integer, intent(OUT) :: noderank  ! rank in PEs on same host communicator            !InTf!
  integer, intent(OUT) :: numarank  ! rank in PEs in same NUMA communicator            !InTf!
  integer, intent(OUT) :: peerrank  ! rank in node/numa peers                          !InTf!
  integer, intent(OUT) :: isiz      ! size of PEs on same host communicator            !InTf!
  integer, intent(OUT) :: ierr      ! error code                                       !InTf!
!******
  integer :: i, rank, lstring, status, nnuma, numapop
  character(len=32) :: force_numa

  noderank = -1
  peerrank = -1
  isiz = 0
  nodecomm = MPI_COMM_NULL
  numacomm = MPI_COMM_NULL
  peercomm = MPI_COMM_NULL

  do i = 1, ncached
    if(corig(i) == origcomm) then  ! cached entry found
      nodecomm = cnode(i)
      numacomm = cnuma(i)
      if(split_by == RPN_MPI_BY_NUMA .or. split_by == RPN_MPI_BY_SOCKET) then
        peercomm = speer(i)
      else
        peercomm = cpeer(i)
      endif
      exit
    endif
  enddo

  if(nodecomm == MPI_COMM_NULL) then          ! nothing useful found in cache
    ncached = min(MAX_CACHE, ncached + 1)
    call mpi_comm_rank(origcomm, rank, ierr)  ! rank in original communicator
    if(ierr .ne. MPI_SUCCESS) return
    corig(ncached)   = origcomm

    ! split on shared memory node basis
    call mpi_comm_split_type(origcomm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, nodecomm, ierr)
    if(ierr .ne. MPI_SUCCESS) return
    cnode(ncached)   = nodecomm
    ! rank of this PE on this SMP node (belonging to origcomm)
    call MPI_Comm_rank(nodecomm, noderank,ierr)
    if(ierr .ne. MPI_SUCCESS) return
    ! split origcomm using noderank as the color to make "same rank on node" (node peers) communicator
    call MPI_Comm_split(origcomm, noderank, rank, peercomm, ierr)
    if(ierr .ne. MPI_SUCCESS) return
    cpeer(ncached)   = peercomm

    if(my_numa == -1) then          ! initialize my_numa
      ! can be used to force number of numa spaces (usually sockets) on node
      ! export RPN_MPI_FORCE_NUMA=number_of_numa_spaces
      call GET_ENVIRONMENT_VARIABLE('RPN_MPI_FORCE_NUMA', force_numa, lstring, status)
      if(status == 0) then
        read(force_numa,*) nnuma
        call MPI_comm_size(nodecomm, numapop, ierr)       ! node population
        numapop = (numapop + nnuma - 1) / nnuma           ! numa space population
        my_numa = noderank / numapop                      ! overriding c_numa_node_of_cpu()
      else
        my_numa = c_numa_node_of_cpu(c_sched_getcpu())    ! get numa space associated with the core this PE is using
      endif
    endif

    ! re split nodecomm by socket (NUMA domain) number
    call mpi_comm_split(nodecomm, my_numa, noderank, numacomm, ierr)
    if(ierr .ne. MPI_SUCCESS) return
    cnuma(ncached)   = numacomm
    ! rank in NUMA domain communicator
    call MPI_Comm_rank(numacomm, numarank,ierr)
    if(ierr .ne. MPI_SUCCESS) return
    ! create socket peers communicator (split original communicator)
    call MPI_Comm_split(origcomm, numarank, rank, peercomm, ierr)      
    if(ierr .ne. MPI_SUCCESS) return
    speer(ncached)   = peercomm
    if(split_by == RPN_MPI_BY_NUMA .or. split_by == RPN_MPI_BY_SOCKET) then
      peercomm = speer(ncached)
    else
      peercomm = cpeer(ncached) 
    endif
  endif

  ! rank of this PE in the peers communicator
  call MPI_Comm_rank(peercomm, peerrank,ierr)
  if(ierr .ne. MPI_SUCCESS) return
  ! number of PEs on this SMP node (belonging to origcomm)
  call MPI_Comm_size(nodecomm, isiz, ierr)
  if(ierr .ne. MPI_SUCCESS) return

  ierr = MPI_SUCCESS
  return
  
end subroutine RPN_MPI_split_by_memory_domain  !InTf!
#if 0
!****f* rpn_mpi/RPN_MPI_split_by_socket
! DESCRIPTION
!   split a communicator on a socket basis (same as NUMA space for this version)
! SYNOPSIS
subroutine RPN_MPI_split_by_socket(origcomm, sockcomm, sockrank, peercomm, peerrank, ierr)   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! IGNORE
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
#include <RPN_MPI_system_interfaces.hf>
#include <RPN_MPI_mpi_symbols.hf>
! ARGUMENTS
  integer, intent(IN)  :: origcomm  ! MPI communicator to split on a socket basis        !InTf!
  integer, intent(OUT) :: sockcomm  ! new communicator to be used py PEs on same socket  !InTf!
  integer, intent(OUT) :: sockrank  ! rank in socket communicator                        !InTf!
  integer, intent(OUT) :: peercomm  ! new communicator for socket peers                  !InTf!
  integer, intent(OUT) :: peerrank  ! rank in socket peers                               !InTf!
  integer, intent(OUT) :: ierr      ! error code                                         !InTf!
!******
  integer :: socket, rank, indx
  integer :: nnuma, my_numa, my_core, lstring, status, numapop
  character(len=32) :: force_numa
  integer, external :: RPN_MPI_split_by_memory_domain

  indx = RPN_MPI_split_by_memory_domain(origcomm, nodecomm, peercomm, noderank, peerrank, isiz, ierr) ! split by node first
  if(ierr .ne. MPI_SUCCESS) return


  return
end subroutine RPN_MPI_split_by_socket  !InTf!
#endif
#if defined(SELF_TEST)
program test_numa
!
! s.f90 -DSELF_TEST -mpi -o a.out RPN_MPI_split_by_memory_domain.F90
! export export RPN_MPI_FORCE_NUMA=4
! r.run_in_parallel -npex 8 -npey 1 -pgm ./a.out -inorder -maxcores
!
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: nodecomm, sockcomm, peercomm, noderank, sockrank, peerrank, isiz
  call mpi_init(ierr)
  call RPN_MPI_split_by_socket(MPI_COMM_WORLD, nodecomm, sockcomm, peercomm, noderank, sockrank, peerrank, isiz, ierr)
  print 1,'noderank, sockrank, peerrank =',noderank, sockrank, peerrank
1 format(A,10I10)
  call mpi_finalize(ierr)
end program
#endif
