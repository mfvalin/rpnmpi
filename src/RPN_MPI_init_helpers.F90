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
!===================================================================
!     legacy support for rpn_comm entry points has been added
!===================================================================
!InTf!
      ! must be called before RPN_MPI_init
      ! provide call_back routine to get domain information
      subroutine RPN_MPI_mydomain (call_back, mydomain)             !InTf!
      use rpn_mpi
      implicit none
      external :: call_back                                          !InTf!
      integer, intent(OUT) :: mydomain                               !InTf!
!
      logical mpi_started
      integer ierr, err, err_all, pe_tot2, pe_me2, npe_per_domain, offset, ndomains
!
      if(RPN_MPI_IS_INITIALIZED .or. RPN_COMM_IS_INITIALIZED) then  ! check legacy RPN_COMM flags
         if (pe_me.eq.0) write (6,1002)
         stop
      endif
      call MPI_INITIALIZED(mpi_started,ierr)
      if (.not. mpi_started ) call MPI_init(ierr)

      ! use WORLD_COMM_MPI instead of MPI_COMM_WORLD
      call mpi_comm_size (WORLD_COMM_MPI,pe_tot2,ierr)
      call mpi_comm_rank (WORLD_COMM_MPI,pe_me2 ,ierr)

      call call_back (ndomains, offset, err)  ! get number of domains and offset from callback

      call mpi_allreduce(err,err_all,1,MPI_INTEGER,MPI_MIN,WORLD_COMM_MPI,ierr)

      if (err_all.lt.0) then
         if (pe_me2.eq.0) write (6,1001) 
         call RPN_MPI_FINALIZE(ierr)
         stop
      endif

      npe_per_domain = pe_tot2 / ndomains
      mydomain       = pe_me2  / npe_per_domain + offset
!
 1001 format(/'===> RPN_MPI_mydomain: FATAL error in call_back    --- ABORTING ---'/)
 1002 format(/'===> RPN_MPI_mydomain: RPN_MPI already initialized --- ABORTING ---'/)
      return
      end subroutine RPN_MPI_mydomain                                !InTf!
!     legacy routine for rpn_comm
      subroutine RPN_COMM_mydomain (call_back, mydomain)             !InTf!
      use rpn_comm
      implicit none
      external :: call_back                                          !InTf!
      integer, intent(OUT) :: mydomain                               !InTf!
      call RPN_MPI_mydomain (call_back, mydomain)
      return
      end subroutine RPN_COMM_mydomain                               !InTf!
!===================================================================
      ! must be called before RPN_MPI_init
      ! set world_comm as "WORLD" communicator instead of MPI_COMM_WORLD
!InTf!
      subroutine RPN_MPI_world_set(world_comm)                    !InTf!
      use rpn_mpi
      implicit none
      integer, intent(IN) ::  world_comm                           !InTf!
!
      if(RPN_MPI_IS_INITIALIZED .or. RPN_COMM_IS_INITIALIZED) then  ! check legacy RPN_COMM flags
         if (pe_me.eq.0) write (6,1002)
         stop
      endif
      WORLD_COMM_MPI=world_comm
      WORLD_COMM_MPI_INIT=.true.
 1002 format(/'===> RPN_MPI_mydomain: RPN_MPI already initialized --- ABORTING ---'/)
      return
      end subroutine RPN_MPI_world_set                            !InTf!
!     legacy routine for rpn_comm
      subroutine RPN_COMM_world_set(world_comm)                    !InTf!
      use rpn_comm
      implicit none
      integer, intent(IN) ::  world_comm                           !InTf!
      call RPN_MPI_world_set(world_comm)
      return
      end subroutine RPN_COMM_world_set                            !InTf!
!===================================================================
      ! get the core number of the CPU this process is running on
      ! get the NUMA space number for the CPU this process is running on
      subroutine RPN_MPI_get_core_and_numa(core, numa)            !InTf!
      use ISO_C_BINDING
      implicit none
      integer, intent(OUT) ::  core, numa                          !InTf!
#include <RPN_MPI_system_interfaces.hf>

      core = sched_get_my_cpu()
      numa = numa_node(core)
      return
      end subroutine RPN_MPI_get_core_and_numa                    !InTf!
!===================================================================
      ! find if RPN_MPI_init has been called
      function RPN_MPI_initialized() result(status)             !InTf!
      use rpn_mpi
      implicit none
      logical ::  status                                        !InTf!
      status = RPN_MPI_IS_INITIALIZED .or. RPN_COMM_IS_INITIALIZED
      end function RPN_MPI_initialized                          !InTf!
      