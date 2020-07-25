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
!     legacy variables from module rpn_comm are initialized to provide
!     compatibility with rpn_comm_init family
!===================================================================
!InTf!
!
! multiple levels initialization function application/supergrid/grid/block
!
      INTEGER FUNCTION RPN_MPI_init &                                !InTf!
      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,AppID,Io)   !InTf!
      use RPN_MPI
      use RPN_MPI_mpi_layout
      implicit none                                                  !InTf!
      external :: Userinit                                           !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      integer, intent(in)    :: MultiGrids                           !InTf!
      integer, intent(in)    :: Grids                                !InTf!
      character(len=*)       :: AppID                                !InTf!
      integer, intent(in)    :: Io                                   !InTf!
!arguments
!  I  Userinit  User routine that will be called by PE 0 to
!               get the processor grid topology if it is not supplied
!               (Pex .eq. 0 ) or (Pey .eq. 0)
!  O  Pelocal    PE rank (number) of local PE in its GRID
!  O  Petotal    Number of PEs in this PEs GRID
! I/O Pex        Number of PEs along the X axis. If Pex=0 upon entry
!                it will be set to the proper value upon exit
! I/O Pey        Number of PEs along the Y axis. If Pey=0 upon entry
!                it will be set to the proper value upon exit
!  I  MultiGrids number of simultaneous MultiGrids (usually 1)
!  I  Grids      number of grids in MultiGrids
!  I  AppID      application ID, character string, 5 chars max
! I/O Io         Io mode (there will be mod(io,10) service (IO) PEs per io/10 compute PEs)
!                Io = 182 : 2 service PEs for 18 compute PEs (2 out of 20 PEs are IO PEs)
!
!notes
!      processor topology common /pe/ will be filled here
!      positions are calculated from 0 (ORIGIN 0)
!
!       this is also intended to be a cleanup/refactoring of RPN_MPI_init_multi_level
!
#include <RPN_MPI_system_interfaces.hf>
      integer ierr, i, j, count, reste, status, core
      integer :: pe_type    ! 0 = compute, 1 = service
      integer :: service, compute
      logical mpi_started
      logical ok, allok
      integer RPN_MPI_petopo, RPN_MPI_version
      character *4096 SYSTEM_COMMAND
      integer my_color,iun
      integer version_marker, version_max
      integer pe_my_location(10)
!       external RPN_MPI_unbind_process
      integer, external :: RPN_MPI_get_a_free_unit, RPN_MPI_set_timeout_alarm
      integer :: ApplID
      character(len=5) :: appid5
!
      if(RPN_COMM_IS_INITIALIZED .or. RPN_MPI_IS_INITIALIZED) then ! ignore with warning message or abort ?
        if(lr%grid%all%rank == 0) &
          write(rpn_u,*) 'ERROR: RPN_MPI/RPN_COMM already initialized, ABORTING execution'
        call mpi_finalize(ierr)
        stop
      endif
!
!     build application (domain) "color" from first 5 (at most) characters of AppID
      ApplID = 0
      do i = 1, min(5 , len(AppID))    ! 6 bits per character ASCII 32-96, case insensitive
        ApplID = ApplID * 64 + and(63, 32 + ichar(AppID(i:i)))  ! i th character
      enddo
      call RPN_MPI_init_mpi_layout   ! initialize NEW style layout structure
      call RPN_MPI_get_core_and_numa(core, lr%numa)  ! get numa space for this PE
      lr%host = get_host_id()   ! get host id for this PE
      compute = Io / 10         ! number of compute PEs in a PE block (compute + service PEs)
      service = mod(Io, 10)     ! number of service(IO) PEs in a PE block
!
!     initialize OLD style variables
      pe_indomm=-1
      pe_indomms=-1
      pe_dommtot=-1
      pe_all_domains=MPI_COMM_NULL
      pe_me_all_domains=-1
      pe_a_domain=MPI_COMM_NULL
      pe_me_a_domain=-1
      pe_multi_grid=MPI_COMM_NULL
      pe_me_multi_grid=-1
      pe_grid=MPI_COMM_NULL
      pe_me_grid=-1
!
      RPN_MPI_IS_INITIALIZED  =.true.
      RPN_COMM_IS_INITIALIZED =.true.      ! legacy to maintain compatibility with RPN_COMM_init family
      if( .not. WORLD_COMM_MPI_INIT ) then ! world not set before, use MPI_COMM_WORLD
        WORLD_COMM_MPI_INIT=.true.
        WORLD_COMM_MPI=MPI_COMM_WORLD
      endif
!
      call MPI_INITIALIZED(mpi_started,ierr)      ! is the MPI library already initialized ?
      status = RPN_MPI_set_timeout_alarm(60)      ! maximum of 60 seconds for MPI_init
      if (.not. mpi_started ) call MPI_init(ierr)
      status = RPN_MPI_set_timeout_alarm(0)       ! timeout reset to infinity (no timeout)
      call RPN_MPI_init_mpi_layout                ! make sure layout and constants structures are properly initialized
!
!     --------------------------------------------------------------------------
!     all applications (domains)
      pe_wcomm=WORLD_COMM_MPI                     ! the UNIVERSE
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)     ! rank in UNIVERSE
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)    ! size of UNIVERSE
!     NEW style
      lr%wrld%all%comm = pe_wcomm
      lr%wrld%all%rank = pe_me
      lr%wrld%all%size = pe_tot
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, lr%wrld%same_node%comm, ierr) ! same node
      call MPI_COMM_RANK(lr%wrld%same_node%comm, lr%wrld%same_node%rank, ierr)
      call MPI_COMM_SIZE(lr%wrld%same_node%comm, lr%wrld%same_node%size, ierr)
      call MPI_COMM_SPLIT(lr%wrld%same_node%comm, lr%numa, pe_me, lr%wrld%same_numa%comm, ierr)                    ! same numa space
      call MPI_COMM_RANK(lr%wrld%same_numa%comm, lr%wrld%same_numa%rank, ierr)
      call MPI_COMM_SIZE(lr%wrld%same_numa%comm, lr%wrld%same_numa%size, ierr)
!     OLD style
      pe_all_domains = pe_wcomm
      pe_me_all_domains = pe_me
      pe_tot_all_domains = pe_tot
      call MPI_COMM_GROUP(pe_all_domains,pe_gr_all_domains,ierr)
!
!     --------------------------------------------------------------------------
      if(pe_me == 0)then                           ! if requested produce a "status" file
        call get_environment_variable("RPN_MPI_MONITOR",SYSTEM_COMMAND)  ! export RPN_MPI_MONITOR=filename
        if(SYSTEM_COMMAND .ne. " ") then
          iun = RPN_MPI_get_a_free_unit()
          open(iun,file=trim(SYSTEM_COMMAND),form='FORMATTED')
          write(iun,'(A)')'mpi_init successful'
          close(iun)
        endif
      endif
!       call RPN_MPI_unbind_process ! unbind processes if requested (FULL_UNBIND environment variable, linux only)
!
!     --------------------------------------------------------------------------
!     set diagnostic mode. if >= 2, some diagnostics are printed (pe topo)
      diag_mode=1
      SYSTEM_COMMAND=" "
      call get_environment_variable("RPN_MPI_DIAG",SYSTEM_COMMAND)  ! export RPN_MPI_DIAG=n
      if( SYSTEM_COMMAND .ne. " " ) read(SYSTEM_COMMAND,*) diag_mode
!
!     --------------------------------------------------------------------------
!     check that all Processes use the same version of RPN_MPI
      version_marker=RPN_MPI_version()
      call mpi_allreduce(version_marker, version_max, 1, MPI_INTEGER, MPI_BOR, WORLD_COMM_MPI, ierr)
      if(version_max .ne. version_marker)then
        write(rpn_u,*) 'ERROR: RPN_MPI version mismatch, PLS recompile, ABORTING execution'
        call mpi_finalize(ierr)
        stop
      endif
!
!     even if environment variable RPN_MPI_DOM is set, ignore it (deprecated feature)
!
!     --------------------------------------------------------------------------
!     split WORLD_COMM_MPI by application(domain) using applid as "color"
!
      my_color = ApplID
      call MPI_COMM_SPLIT(WORLD_COMM_MPI, my_color, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my appl
      if(pe_me .eq. 0 .and. diag_mode .ge. 1) then
        write(rpn_u,1000)trim(AppID),pe_tot,MultiGrids,Grids
1000    format('domain=',A,I6,' PEs,',I3,' SuperGrids with',I3,' Grids each')
      endif
      ok = .true.                             
      i=((pe_tot/MultiGrids)/Grids)           ! number of PEs in a grid
      j=pe_tot - i*MultiGrids*Grids           ! remainder
      ok = (j .eq. 0)                         ! check that pe_tot is a multiple of MultiGrids and Grids
      if(.not. ok .and. pe_me .eq. 0) &
         write(rpn_u,*)'ERROR: number of PEs in domain is not a multiple of MultiGrids*Grids, ABORTING execution'
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok .and. pe_me .eq. 0) then
           call RPN_MPI_finalize(ierr)
           stop
      endif
      RPN_MPI_init = my_color
!     NEW style
      lr%appl%all%comm = pe_wcomm
      lr%appl%all%rank = pe_me
      lr%appl%all%size = pe_tot
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, lr%appl%same_node%comm, ierr) ! same node
      call MPI_COMM_RANK(lr%appl%same_node%comm, lr%appl%same_node%rank, ierr)
      call MPI_COMM_SIZE(lr%appl%same_node%comm, lr%appl%same_node%size, ierr)
      call MPI_COMM_SPLIT(lr%appl%same_node%comm, lr%numa, pe_me, lr%appl%same_numa%comm, ierr)                    ! same numa space
      call MPI_COMM_RANK(lr%appl%same_numa%comm, lr%appl%same_numa%rank, ierr)
      call MPI_COMM_SIZE(lr%appl%same_numa%comm, lr%appl%same_numa%size, ierr)
      lr%colors(1)     = my_color
!     OLD style
      my_colors(1) = my_color
      pe_a_domain = pe_wcomm
      pe_me_a_domain = pe_me
      pe_tot_a_domain = pe_tot
      call MPI_COMM_GROUP(pe_a_domain,pe_gr_a_domain,ierr)
if(pe_me == 0) print *,'application split done'
!
!     --------------------------------------------------------------------------
!     if more than 1 supergrid in application, split application into supergrids
!     (compute and service PEs at this point)
!
      my_color = 0
      pe_wcomm = lr%appl%all%comm
      if(MultiGrids .gt. 1) then
        my_color=lr%appl%all%rank / (lr%appl%all%size / MultiGrids)
        RPN_MPI_init = my_color
        call MPI_COMM_SPLIT(lr%appl%all%comm, my_color, lr%appl%all%rank, pe_wcomm, ierr)
      endif
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                 ! rank in supergrid
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                ! size of supergrid
!     NEW style
      lr%sgrd%all%comm = pe_wcomm                             ! supergrid communicator
      lr%sgrd%all%rank = pe_me                                ! ordinal in supergrid communicator
      lr%sgrd%all%size = pe_tot                               ! population in supergrid communicator
      lr%colors(2)     = my_color
      lr%sgrd%row%comm = MPI_COMM_NULL                        ! row and column are not defined for supergrids
      lr%sgrd%row%rank = -1
      lr%sgrd%row%size = -1
      lr%sgrd%column%comm = MPI_COMM_NULL
      lr%sgrd%column%rank = -1
      lr%sgrd%column%size = -1
!
!     supergrid to supergrid peers in application (PEs with same rank in supergrid)
!
      call MPI_COMM_SPLIT(lr%appl%all%comm, pe_me, lr%appl%all%rank, lr%sgrd%grid_peer%comm, ierr)  ! supergrid peers
      call MPI_COMM_RANK(lr%sgrd%grid_peer%comm, lr%sgrd%grid_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%sgrd%grid_peer%comm, lr%sgrd%grid_peer%size, ierr)
!
!     split supergrid communicator into same node and node peers
!
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, lr%sgrd%same_node%comm, ierr)  ! same node
      call MPI_COMM_RANK(lr%sgrd%same_node%comm, lr%sgrd%same_node%rank, ierr)
      call MPI_COMM_SIZE(lr%sgrd%same_node%comm, lr%sgrd%same_node%size, ierr)
      call MPI_COMM_SPLIT(pe_wcomm, lr%sgrd%same_node%rank, lr%sgrd%all%rank, lr%sgrd%node_peer%comm, ierr)         ! node peers
      call MPI_COMM_RANK(lr%sgrd%node_peer%comm, lr%sgrd%node_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%sgrd%node_peer%comm, lr%sgrd%node_peer%size, ierr)
!
!     split same node into same NUMA space and supergrid into NUMA space peers
!
      call MPI_COMM_SPLIT(lr%sgrd%same_node%comm, lr%numa, pe_me, lr%sgrd%same_numa%comm, ierr)             ! same numa space
      call MPI_COMM_RANK(lr%sgrd%same_numa%comm, lr%sgrd%same_numa%rank, ierr)
      call MPI_COMM_SIZE(lr%sgrd%same_numa%comm, lr%sgrd%same_numa%size, ierr)
      call MPI_COMM_SPLIT(pe_wcomm, lr%sgrd%same_numa%rank, lr%sgrd%all%rank, lr%sgrd%numa_peer%comm, ierr) ! numa peers
      call MPI_COMM_RANK(lr%sgrd%numa_peer%comm, lr%sgrd%numa_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%sgrd%numa_peer%comm, lr%sgrd%numa_peer%size, ierr)
!     --------------------------------------------------------------------------
!     TODO : take care of IO processors  (at the grid level)
!            "grid" will have to be split into compute and IO processes
!            will need to look at RPN_MPI_IO_CONFIG environment variable
!            may have to redefine pe_wcomm to be compute PEs
!            do we instead want IO at the supergrid level ?
!            grids will have to cooperate at the supergrid level
!            my_colors for compute PEs only ?
!     options:
!          return to the caller with special code
!          call user supplied subroutine
!          call own IO server code
!     --------------------------------------------------------------------------
!
!     split supergrid into compute and service (IO) PEs
!
      pe_type = 0   ! compute by default, service PEs out of every (service+compute) are service PEs
      if( mod(pe_me,compute+service) < service ) pe_type = 1    ! determine if compute or service PE
! print *,'PE type =',pe_type
!
      call MPI_COMM_SPLIT(lr%sgrd%all%comm, pe_type, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                   ! rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                  ! size
      if(pe_type == 0) then                                     ! compute PEs
!           NEW style
        lr%sgrd%compute%comm = pe_wcomm
        lr%sgrd%compute%rank = pe_me
        lr%sgrd%compute%size = pe_tot
        lr%sgrd%service%comm = MPI_COMM_NULL
        lr%sgrd%service%rank = -1
        lr%sgrd%service%size = -1
      else                                                      ! service PEs
        lr%sgrd%compute%comm = MPI_COMM_NULL
        lr%sgrd%compute%rank = -1
        lr%sgrd%compute%size = -1
        lr%sgrd%service%comm = pe_wcomm
        lr%sgrd%service%rank = pe_me
        lr%sgrd%service%size = pe_tot
      endif
!     OLD style, compute PEs only
      my_colors(2) = my_color                                 ! my multigrid number
      pe_multi_grid = pe_wcomm                                ! multigrid communicator
      pe_me_multi_grid = pe_me                                ! my rank in multigrid
      pe_tot_multi_grid = pe_tot                              ! multigrid population
      call MPI_COMM_GROUP(pe_multi_grid,pe_gr_multi_grid,ierr)
      pe_indomms = pe_multi_grid                              ! multigrid communicator
      call MPI_COMM_GROUP(pe_indomms,pe_gr_indomms,ierr)
! print *,'supergrid split done'
!     --------------------------------------------------------------------------
!     if more than 1 grid per supergrid, split supergrid into grids 
!     (compute and service PEs at this point)
!
      my_color = 0
      pe_wcomm = lr%sgrd%all%comm
      if(Grids .gt. 1) then
        my_color=lr%sgrd%all%rank / (lr%sgrd%all%size / Grids)
        RPN_MPI_init = my_color
        call MPI_COMM_SPLIT(lr%sgrd%all%comm, my_color, lr%sgrd%all%rank, pe_wcomm, ierr)
      endif
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
!     NEW style
      lr%grid%all%comm = pe_wcomm
      lr%grid%all%rank = pe_me
      lr%grid%all%size = pe_tot
      lr%colors(3)     = my_color
!
!     grid to grid peers in supergrid (PEs with same rank in grid)
!
! print *,'splitting for grid peers'
      call MPI_COMM_SPLIT(lr%sgrd%all%comm, pe_me, lr%sgrd%all%rank, lr%grid%grid_peer%comm, ierr)  ! grid peers
      call MPI_COMM_RANK(lr%grid%grid_peer%comm, lr%grid%grid_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%grid_peer%comm, lr%grid%grid_peer%size, ierr)
!     OLD style
      pe_grid_peers = lr%grid%grid_peer%comm
      call MPI_COMM_SIZE(pe_grid_peers,pe_tot_peer,ierr)  ! This must be equal to the number of grids
      call MPI_COMM_RANK(pe_grid_peers,pe_me_peer,ierr)   ! in a supergrid (or else...)
      call MPI_COMM_GROUP(pe_grid_peers,pe_gr_grid_peers,ierr)
!
!     split grid into nodes and node peers
!
! print *,'splitting for same node'
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, lr%grid%same_node%comm, ierr)  ! same node
      call MPI_COMM_RANK(lr%grid%same_node%comm, lr%grid%same_node%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%same_node%comm, lr%grid%same_node%size, ierr)
! print *,'splitting for node peers'
      call MPI_COMM_SPLIT(pe_wcomm, lr%grid%same_node%rank, lr%sgrd%all%rank, lr%grid%node_peer%comm, ierr) ! node peers
      call MPI_COMM_RANK(lr%grid%node_peer%comm, lr%grid%node_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%node_peer%comm, lr%grid%node_peer%size, ierr)
!     OLD style
      pe_grid_host = lr%grid%same_node%comm
      call MPI_COMM_RANK(pe_grid_host,pe_me_grid_host,ierr)    ! my rank on this host
      call MPI_COMM_SIZE(pe_grid_host,pe_tot_grid_host,ierr)   ! population of this host
      call MPI_COMM_GROUP(pe_grid_host,pe_gr_grid_host,ierr)   ! group communicator
!
!     split grid into NUMA spaces and NUMA peers
!
! print *,'splitting for same numa'
      call MPI_COMM_SPLIT(lr%grid%same_node%comm, lr%numa, pe_me, lr%grid%same_numa%comm, ierr)                     ! same numa
      call MPI_COMM_RANK(lr%grid%same_numa%comm, lr%grid%same_numa%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%same_numa%comm, lr%grid%same_numa%size, ierr)
! print *,'splitting for numa peers'
      call MPI_COMM_SPLIT(lr%grid%all%comm, lr%grid%same_numa%rank, lr%grid%all%rank, lr%grid%numa_peer%comm, ierr) ! numa peers
      call MPI_COMM_RANK(lr%grid%numa_peer%comm, lr%grid%numa_peer%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%numa_peer%comm, lr%grid%numa_peer%size, ierr)
!
!     split grid into compute and service (IO) PEs
!
!     pe_type already set above when dealing with supergrids
! print *,'splitting for service'
      call MPI_COMM_SPLIT(lr%grid%all%comm, pe_type, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                   ! rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                  ! size
      if(pe_type == 0) then                                     ! compute PEs
!           NEW style
        lr%grid%compute%comm = pe_wcomm
        lr%grid%compute%rank = pe_me
        lr%grid%compute%size = pe_tot
        lr%grid%service%comm = MPI_COMM_NULL
        lr%grid%service%rank = -1
        lr%grid%service%size = -1
      else                                                      ! service PEs
        lr%grid%compute%comm = MPI_COMM_NULL
        lr%grid%compute%rank = -1
        lr%grid%compute%size = -1
        lr%grid%service%comm = pe_wcomm
        lr%grid%service%rank = pe_me
        lr%grid%service%size = pe_tot
      endif
!     OLD style, compute PEs only
      my_colors(3) = my_color                                 ! my grid number in my multigrid
      pe_grid = pe_wcomm                                      ! grid communicator
      pe_me_grid = pe_me                                      ! my rank in grid
      pe_tot_grid = pe_tot                                    ! grid population
      call MPI_COMM_GROUP(pe_grid,pe_gr_grid,ierr)

      Pelocal = pe_me                                          ! my rank in my grid
      Petotal = pe_tot                                         ! number of pes in my grid

      pe_indomm = pe_grid
      pe_gr_indomm = pe_gr_grid
      pe_pe0 = 0
      pe_dommtot = pe_tot
      call MPI_COMM_GROUP(pe_wcomm,pe_gr_wcomm,ierr)
! if(pe_me == 0) print *,'grid split done'
!     --------------------------------------------------------------------------
!     Grid initialization (compute PEs) 
!     get PEs along X and Y
!     call Userinit if appropriate
!
      ok = .true.
      if(pe_me .eq. pe_pe0) then
        if ( Pex.eq.0 .or. Pey.eq.0  ) then ! get processor topology from Userinit
          WORLD_pe(1)=pe_tot
          WORLD_pe(2)=1
          call Userinit(WORLD_pe(1),WORLD_pe(2))
          if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
            ok = .false.
            write(rpn_u,*) 'ERROR: (RPN_MPI_init) Inconsistency between'
            write(rpn_u,*) '       userinit Subroutine and total number of PEs'
            write(rpn_u,*) '       please double check topology'
          endif
          if(diag_mode.ge.1) then
            write(rpn_u,*)'INFO: Requested topology = ',WORLD_pe(1),' by ',WORLD_pe(2)
            write(rpn_u,*)'      Grid will use ',pe_dommtot,' processes'
          endif
        else ! ( Pex.ne.0 .and. Pey.ne.0  )
          write(rpn_u,*) 'INFO: (RPN_MPI_init) Forced topology :',Pex,' by',Pey
          WORLD_pe(1) = Pex
          WORLD_pe(2) = Pey
          if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
            ok = .false.
            write(rpn_u,*) 'ERROR: (RPN_MPI_init) Inconsistency between Pex, Pey and total number of PEs'
            write(rpn_u,*) '       please double check topology'
          endif
          if(diag_mode.ge.1) then
            write(rpn_u,*)'Requested topology =',WORLD_pe(1),' by ',WORLD_pe(2)
          endif
        endif ! ( Pex.eq.0 .or. Pey.eq.0  )
!
        if(WORLD_pe(1)*WORLD_pe(2) .gt. pe_dommtot) then
          write(rpn_u,*)' ERROR: not enough PEs for requested decomposition '
          write(rpn_u,*)'        REQUESTED=',WORLD_pe(1)*WORLD_pe(2)
          write(rpn_u,*)'        AVAILABLE=',pe_dommtot
          ok = .false.
        endif
      endif  ! (pe_me .eq. pe_pe0)
!
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL, MPI_LAND, WORLD_COMM_MPI, ierr)
      if(.not.allok) then
        if(.not. ok .and. pe_me .eq. pe_pe0 ) write(rpn_u,*)'ERROR: problem in grid initialization'
        call RPN_MPI_finalize(ierr)
        stop
      endif
!
!      send WORLD topology to all PEs. That will allow all PEs
!      to compute other PE topology parameters locally.
!       for doing this, we need to define some basic domains
!       communicators.

      call MPI_COMM_rank(pe_indomm,pe_medomm,ierr)
      pe_defcomm = pe_indomm
      pe_defgroup = pe_gr_indomm
!
!       broadcast number of PEs along X and Y axis
!       broadcast PE block sizes (deltai and deltaj)
!
      WORLD_pe(3) = deltai
      WORLD_pe(4) = deltaj
      call MPI_BCAST(WORLD_pe, size(WORLD_pe), MPI_INTEGER, 0, pe_indomm, ierr)
      pe_nx  = WORLD_pe(1)
      pe_ny  = WORLD_pe(2)
      deltai = WORLD_pe(3)
      deltaj = WORLD_pe(4)
!
      if ( Pex.eq.0 .or. Pey.eq.0  ) then ! return processor topology
        Pex = WORLD_pe(1)
        Pey = WORLD_pe(2)
      endif
!
!      pe_pe0 is not equal to 0 if there are more than one domain
!      computational grid
!
      count = pe_pe0
!
!      fill tables containing the position along the X axis (pe_xtab)
!      and along the Y axis (pe_ytab) for all processors
!     --------------------------------------------------------------------------
!     PE topology
      ierr = RPN_MPI_petopo(WORLD_pe(1),WORLD_pe(2))
      lr%grid%row%comm    = pe_myrow
      call MPI_COMM_RANK(lr%grid%compute%comm, lr%grid%row%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%compute%comm, lr%grid%row%size, ierr)
      lr%grid%column%comm = pe_mycol
      call MPI_COMM_RANK(lr%grid%compute%comm, lr%grid%column%rank, ierr)
      call MPI_COMM_SIZE(lr%grid%compute%comm, lr%grid%column%size, ierr)

      pe_my_location(1) = pe_mex
      pe_my_location(2) = pe_mey
      pe_my_location(3) = pe_me_grid
      pe_my_location(4) = my_colors(3)
      pe_my_location(5) = pe_me_multi_grid
      pe_my_location(6) = my_colors(2)
      pe_my_location(7) = pe_me_a_domain
      pe_my_location(8) = my_colors(1)
      pe_my_location(9) = lr%host
      pe_my_location(10)= lr%numa

      pe_tot = lr%wrld%all%size
print *,'pe_tot =',pe_tot
! diag_mode = 3
      allocate(pe_location(10,0:pe_tot-1))
      call MPI_allgather( pe_my_location,10,MPI_INTEGER, &
                          pe_location,   10,MPI_INTEGER, &
                          WORLD_COMM_MPI, ierr)
      if( pe_me_all_domains .eq. 0 .and. diag_mode .ge.3) then
        write(rpn_u,*)'                         FULL PE MAP'
        write(rpn_u,*)'    mex     mey   me(g)    grid  me(sg)   sgrid   me(d)  application appid     host   numa'
        do j=0,pe_tot_all_domains-1
          appid5 = '     '
          reste = pe_location(8,j)
          do i = 5, 1, -1
            appid5(i:i) = char(32 + mod(reste,64))
            reste = reste / 64
          enddo
          write(rpn_u,1001)pe_location(1:8,j),appid5,pe_location(9:10,j)
1001      format(7I8,I12,3X,A5,Z9,I3)
        enddo
      endif
      deallocate(pe_location)
!     --------------------------------------------------------------------------
!     PE blocks, initialized to 1 x 1 blocks
      BLOC_EXIST   =.false.
      BLOC_SIZEX   = 1
      BLOC_SIZEY   = 1
      BLOC_mybloc  = 0
      BLOC_myblocx = 0
      BLOC_myblocy = 0
      BLOC_me      = pe_me
      BLOC_comm_world = pe_indomm
      BLOC_comm_row = pe_myrow
      BLOC_comm_col = pe_mycol
      BLOC_corner = pe_pe0
      BLOC_master = 0
      if(pe_me.eq.pe_pe0) then
        BLOC_master=1
      endif
      pe_bloc = pe_indomm

      call MPI_Group_incl(pe_gr_indomm, 1, 0, pe_gr_blocmaster, ierr) 
      call MPI_Comm_create(pe_indomm,pe_gr_blocmaster, pe_blocmaster, ierr)

!       for each communicator, store the corresponding group
!
      call MPI_COMM_GROUP(pe_myrow,pe_gr_myrow,ierr)
      call MPI_COMM_GROUP(pe_mycol,pe_gr_mycol,ierr)
      call MPI_COMM_GROUP(pe_bloc,pe_gr_bloc,ierr)
      return
      end FUNCTION RPN_MPI_init                        !InTf!
