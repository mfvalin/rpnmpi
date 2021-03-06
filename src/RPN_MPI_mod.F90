!/* rpn_mpi - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
module rpn_mpi_comm_0
  use iso_c_binding
  use rpn_mpi_mpif
  implicit none
  save
!       ---------------------------------------------------------------------
!       this module is used by rpn_comm and rpn_mpi 
!       ---------------------------------------------------------------------
!
! control variables
!
  integer, target :: WORLD_COMM_MPI=MPI_COMM_WORLD                ! plays the role of MPI_COMM_WORLD for rpn_mpi
  integer :: diag_mode
  logical :: WORLD_COMM_MPI_INIT=.false.
  logical :: RPN_MPI_IS_INITIALIZED=.false.
  logical :: RPN_COMM_IS_INITIALIZED=.false.
  integer :: deltai=1   ! deltai and deltaj are used by rpn_mpi_petopo to distribute grid over PEs
  integer :: deltaj=1   ! default PE distribution is X increasing, then Y increasing
!
! Output file unit
!
  integer :: rpn_u = 6
!
!        GLOBAL information, will be BROADCAST
!
!        WORLD_pe(1) number of PEs along x in grid
!        WORLD_pe(2) number of PEs along y in grid
!        WORLD_pe(3) deltai, size of PE blocks along x
!        WORLD_pe(4) deltaj, size of PE blocks along Y
!        WORLD_pe(5:10) provision for future expansion
!
  integer WORLD_pe(10)
!
!   logical :: async_exch=.true.        ! asynchronous halo exchange (level 1)
  logical :: full_async_exch=.false.  ! fully asynchronous halo exchange (level 2)
  logical :: rpn_ew_ext_l = .false.   ! extended halo option (haloy extra rows on North and South tiles)
!
! domain boundary flags, LOGICAL
! .true. if a PE(TILE) is on a domain edge
! all 4 values .true. if single tile
! normally set by routine rpn_mpi_init
!
  logical bnd_east,bnd_west,bnd_north,bnd_south
end module rpn_mpi_comm_0

module rpn_mpi_comm_1
  use iso_c_binding
  implicit none
!
!       ---------------------------------------------------------------------
!       the following block is used mainly by rpn_comm,
!       but is filled by rpn_mpi for legacy support
!       ---------------------------------------------------------------------
!        characteristics of the local PE
!        plus some PE grid topology information
!        normally set by routine rpn_mpi_init
!
!        pe_me    PE number of this PE
!        pe_mex   x coordinate of this PE in domain (origin 0)
!        pe_mey   y coordinate of this PE in domain (origin 0)
!        pe_myrow communicator for PEs in same ROW(same pe_mey)
!        pe_mycol communicator for PEs in same COLUMN(same pe_mex)
!        pe_tot   total number of PEs involved in a grid
!        pe_nx    number of PEs along x axis in a grid
!        pe_ny    number of PEs along y axis in a grid
!        pe_pe0   PE number of first PE in a grid (normally zero)
!        pe_extra flag =1 if pe in compute grid, =0 otherwise
!        pe_grid   communicator for PEs in grid and out of grid
!       pe_ingrid  communicator for PEs in grid
!       pe_outgrid communicator for PEs out of grid
!       pe_bloc    communicator for internal blocs in grid
!       pe_blocmaster communicator for bloc_corner PEs in grid
!       pe_id       matrix of PE numbers in grid
!
!                 pe_me=pe_id(pe_mex,pe_mey)
!
!
  integer :: pe_me,pe_mex,pe_mey,pe_myrow,pe_mycol
  integer :: pe_tot,pe_nx,pe_ny,pe_pe0,pe_extra
  integer :: pe_gr_extra,pe_gr_myrow,pe_gr_mycol
  integer :: pe_bloc, pe_blocmaster, pe_defgroup
  integer :: pe_gr_bloc, pe_gr_blocmaster, pe_defcomm
  integer :: pe_gr_indomm, pe_gr_outdomm, pe_indomm, pe_outdomm
  integer :: pe_gr_indomms, pe_indomms ! multigrid countepart to pe_gr_indomm and pe_indomm
  integer :: pe_wcomm, pe_gr_wcomm
  integer :: pe_wcomms, pe_gr_wcomms  ! multigrid countepart to pe_wcomm and pe_gr_wcomm
  integer :: pe_dommtot, pe_medomm
  integer :: pe_all_domains, pe_gr_all_domains      ! all the domains
  integer :: pe_me_all_domains, pe_tot_all_domains  ! all the domains
  integer :: pe_a_domain, pe_gr_a_domain            ! all multigrids in a domain
  integer :: pe_me_a_domain, pe_tot_a_domain        ! all multigrids in a domain
  integer :: pe_multi_grid, pe_gr_multi_grid        ! all the grids in a multigrid
  integer :: pe_me_multi_grid, pe_tot_multi_grid    ! all the grids in a multigrid
  integer :: pe_grid, pe_gr_grid                    ! a single grid
  integer :: pe_me_grid, pe_tot_grid                ! a single grid
  integer :: pe_grid_peers, pe_gr_grid_peers        ! PEs with same rank on different grids of same multigrid
  integer :: pe_me_peer, pe_tot_peer                ! PEs with same rank on different grids of same multigrid
  integer :: pe_grid_host, pe_gr_grid_host          ! PEs on same host node (belonging to same "grid")
  integer :: pe_me_grid_host, pe_tot_grid_host      ! PEs on same host node (belonging to same "grid")
  integer :: my_colors(3)                           ! domain ordinal/multigrid ordinal/grid ordinal (color)
  integer, allocatable, dimension(:,:) :: pe_id     !  O( pe_id(pe_nx,pe_ny) )
  integer, allocatable, dimension(:)   :: pe_xtab,pe_ytab  ! O(pe_xtab(pe_nx)) O(pe_ytab(pe_ny)) 
  integer, allocatable, dimension(:,:) :: pe_location   ! pe_x,pe_x,pe_ingrid,pe_insgrid,pe_indomain pe_location(8,total_nb_of_pes)
!
  integer ord_max
  integer, allocatable, dimension(:) :: ord_tab
end module rpn_mpi_comm_1

module rpn_mpi_comm_2
  use iso_c_binding
  implicit none
!       ---------------------------------------------------------------------
!       options
!       ---------------------------------------------------------------------
!
!        pe_optn  table of options ( option name)
!        pe_opiv  values of integer options
!        pe_oprv  values of real options
!        pe_opcv  values of character options
  integer, parameter :: MAX_OPTN=10
  character *4 pe_optn(MAX_OPTN)
  integer :: pe_opiv(MAX_OPTN)
  real *4 pe_oprv(MAX_OPTN)
  character *4 pe_opcv(MAX_OPTN)
  integer, private :: RESTE
  parameter (RESTE=MAX_OPTN-1)
  data pe_optn/'FILL',RESTE*'    '/
  data pe_opcv/MAX_OPTN*'    '/
  data pe_oprv/MAX_OPTN*0.0/
  data pe_opiv/MAX_OPTN*0/
end module rpn_mpi_comm_2

module rpn_mpi_comm_3
  use iso_c_binding
  implicit none
!
!       Subgrid information (blocks)
!
  logical BLOC_EXIST
  integer BLOC_SIZEX,BLOC_SIZEY, BLOC_MASTER
  integer BLOC_mybloc,BLOC_myblocx,BLOC_myblocy
  integer BLOC_me,BLOC_corner
  integer BLOC_comm_world, bloc_comm_row, bloc_comm_col
end module rpn_mpi_comm_3

!       ---------------------------------------------------------------------
!       this module is used by rpn_comm
!       ---------------------------------------------------------------------
module rpn_mpi_comm_4
  use iso_c_binding
  use rpn_mpi_mpif
  implicit none
  include 'RPN_COMM_constants.inc'
!
!  symbol tables
!
  type :: SYMTAB
    integer :: number
    character (len=32) :: string
  end type
!
  type(SYMTAB), dimension(6) :: &  ! miscellaneous symbols
  misc_tab=(/ &
  SYMTAB(MPI_ANY_SOURCE,RPN_COMM_ANY_SOURCE), &
  SYMTAB(MPI_ANY_TAG,RPN_COMM_ANY_TAG), &
  SYMTAB(MPI_COMM_NULL,RPN_COMM_COMM_NULL), &
  SYMTAB(MPI_GROUP_NULL,RPN_COMM_GROUP_NULL), &
  SYMTAB(MPI_COMM_WORLD,RPN_COMM_COMM_WORLD), &
  SYMTAB(MPI_SUCCESS,RPN_COMM_SUCCESS)  &
  /)
  type(SYMTAB), dimension(25) :: &  ! data types
  type_tab=(/ &
  SYMTAB(MPI_BYTE,RPN_COMM_BYTE),  &
  SYMTAB(MPI_PACKED,RPN_COMM_PACKED),  &
  SYMTAB(MPI_UB,RPN_COMM_UB),  &
  SYMTAB(MPI_LB,RPN_COMM_LB),  &
  SYMTAB(MPI_CHARACTER,RPN_COMM_CHARACTER),  &
  SYMTAB(MPI_LOGICAL,RPN_COMM_LOGICAL),  &
  SYMTAB(MPI_INTEGER,RPN_COMM_INTEGER),  &
  SYMTAB(MPI_INTEGER1,RPN_COMM_INTEGER1),  &
  SYMTAB(MPI_INTEGER2,RPN_COMM_INTEGER2),  &
  SYMTAB(MPI_INTEGER4,RPN_COMM_INTEGER4),  &
  SYMTAB(MPI_INTEGER8,RPN_COMM_INTEGER8),  &
!        SYMTAB(MPI_INTEGER16,RPN_COMM_INTEGER16),  &
  SYMTAB(MPI_DATATYPE_NULL,RPN_COMM_DATATYPE_NULL), &
  SYMTAB(MPI_REAL,RPN_COMM_REAL),  &
  SYMTAB(MPI_REAL4,RPN_COMM_REAL4),  &
  SYMTAB(MPI_REAL8,RPN_COMM_REAL8),  &
  SYMTAB(MPI_REAL16,RPN_COMM_REAL16),  &
  SYMTAB(MPI_DOUBLE_PRECISION,RPN_COMM_DOUBLE_PRECISION),  &
  SYMTAB(MPI_COMPLEX,RPN_COMM_COMPLEX),  &
  SYMTAB(MPI_COMPLEX8,RPN_COMM_COMPLEX8),  &
  SYMTAB(MPI_COMPLEX16,RPN_COMM_COMPLEX16),  &
  SYMTAB(MPI_COMPLEX32,RPN_COMM_COMPLEX32),  &
  SYMTAB(MPI_DOUBLE_COMPLEX,RPN_COMM_DOUBLE_COMPLEX),  &
  SYMTAB(MPI_2REAL,RPN_COMM_2REAL),  &
  SYMTAB(MPI_2DOUBLE_PRECISION,RPN_COMM_2DOUBLE_PRECISION),  &
  SYMTAB(MPI_2INTEGER,RPN_COMM_2INTEGER)  &
!        SYMTAB(MPI_2COMPLEX,RPN_COMM_2COMPLEX),  &
!        SYMTAB(MPI_2DOUBLE_COMPLEX,RPN_COMM_2DOUBLE_COMPLEX),  &
  /)
  type(SYMTAB), dimension(14) :: &  ! operation types
  op_tab=(/ &
  SYMTAB(MPI_MAX,RPN_COMM_MAX),  &
  SYMTAB(MPI_MIN,RPN_COMM_MIN),  &
  SYMTAB(MPI_SUM,RPN_COMM_SUM),  &
  SYMTAB(MPI_PROD,RPN_COMM_PROD),  &
  SYMTAB(MPI_LAND,RPN_COMM_LAND),  &
  SYMTAB(MPI_BAND,RPN_COMM_BAND),  &
  SYMTAB(MPI_LOR,RPN_COMM_LOR),  &
  SYMTAB(MPI_BOR,RPN_COMM_BOR),  &
  SYMTAB(MPI_LXOR,RPN_COMM_LXOR),  &
  SYMTAB(MPI_BXOR,RPN_COMM_BXOR),  &
  SYMTAB(MPI_MAXLOC,RPN_COMM_MAXLOC),  &
  SYMTAB(MPI_MINLOC,RPN_COMM_MINLOC),  &
  SYMTAB(MPI_OP_NULL,RPN_COMM_OP_NULL), &
  SYMTAB(MPI_REPLACE,RPN_COMM_REPLACE)  &
  /)
end module rpn_mpi_comm_4