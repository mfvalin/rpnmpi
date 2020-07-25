! RPN_MPI - Library of useful routines for C and FORTRAN programming
! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
!                          Environnement Canada
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
#define IN_RPN_MPI_cube_topo
#define IN_RPN_MPI_LIBRARY
#if 0
    =========  The "CUBE"  =========

    6 * (M+1) * (M+1) elements in cube

    (M+1)*(M+1)       elements in a cube face

    M+1               elements along a cube edge

    each of the 6 cube faces has M+1 elements along x, and M+1 elements along y

    cube map

    [F](x,y) : face F, element(x,y)
    (x,y,F)  : face F, element(x,y)

    [n]  : face n (0->5)

    <e>  : edge e (0->11)

    edges 0, 1, 2,  3  Yaxis <-> Yaxis
          0  :  [0](0,y) <-> [3](M,y)
          1  :  [1](0,y) <-> [0](M,y)
          2  :  [2](0,y) <-> [1](M,y)
          3  :  [3](0,y) <-> [2](M,y)

    edges 5, 7, 9, 11  Xaxis <-> Xaxis
          5  :  [0](x,M) <-> [4](x  ,0)
          7  :  [2](x,M) <-> [4](M-x,M)   (axis reversal) (possible PE remapping)
          9  :  [0](x,0) <-> [5](x  ,M)
        11  :  [2](x,0) <-> [5](M-x,0)   (axis reversal) (possible PE remapping)

    edges 4, 6, 8, 10  Xaxis <-> Yaxis (possible npex / npey remapping)
          4  :  [3](x,M) <-> [4](0,M-x)   (axis reversal)
          6  :  [1](x,M) <-> [4](M,  x)
          8  :  [3](x,0) <-> [5](0,  x)
        10  :  [1](x,0) <-> [5](M,M-x)   (axis reversal)

    along x, M+1 elements are distributed across npex PEs
    along y, M+1 elements are distributed across npey PEs
    Xaxis <-> Xaxis will need PE remapping IF 
              M+1 is not a MULTIPLE of npex 
              AND
              there is an axis reversal
    Xaxis <-> Yaxis will need PE remapping IF 
              M+1 is not a MULTIPLE of npex and npey 
              OR
              npex is not equal to npey
          
    -> and <- facing each other means axis reversal during edge exchange

                                      (M-x,M,2)<-
                                  +------<7>------+
                                  |  ( x ,M,4)->  |
                                  |               |
                       (M-y,M,3)<-|->(0,y,4)      |
                                  |               |
                                 <4>     [4]     <6>
                                  |               |
                                  |      (M,y,4)->|->(y,M,1)
                                  |               |
                      (0,M-x,4)<- |    (x,0,4)->  |    (M,x,4)->     (M-x,M,4)<-
                  +-------<4>-----+------<5>------+------<6>------+------<7>------+
                  |   (x, M ,3)-> |    (x,M,0)->  |    (x,M,1)->  |    (x,M,2)->  |
                  |               |               |               |               |
                  |      (M,y,3)->|->(0,y,0)      |      (M,y,1)->|->(0,y,2)      |
                  |               |               |               |               |
                 <3>      [3]    <0>     [0]     <1>     [1]     <2>     [2]     <3>
                  |               |               |               |               |
         (M,y,2)->|->(0,y,3)      |      (M,y,0)->|->(0,y,1)      |      (M,y,2)->|->(0,y,3)
                  |               |               |               |               |
                  |     (x,0,3)-> |    (x,0,0)->  |  (x, 0 ,1)->  |  ( x ,0,2)->  |
                  +-------<8>-----+------<9>------+----<10>-------+-----<11>------+
                        (0,x,5)-> |    (x,M,5)->  |  (M,M-x,5)<-     (M-x,0,5)<-
                                  |               |
                         (y,0,3)->|->(0,y,5)      |
                                  |               |
                                 <8>     [5]     <10>
                                  |               |
                                  |      (M,y,5)->|<-(M-y,0,1)
                                  |               |
                                  |  ( x ,0,5)->  |
                                  +-----<11>------+
                                      (M-x,0,2)<-

    element coordinates (global) in a cube face

              M   +---------------+
              ^   |(0,M)     (M,M)|
              |   |               |
           Y axis |               |
              |   |               |
              |   |               |
              0   |(0,0)     (M,0)|
                  +---------------+
                  ---------------->
                  0    X axis     M

#endif
!
! set PE cube topology
! it is assumed that topo%cube contains a VALID communicator/size/rank combo upon entry
!
 subroutine RPN_MPI_cube_topo(topo, npe, npex, npey, blkx, blky, x_first, ierr)   !InTf!
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
#include <RPN_MPI.hf>
!! import :: RPN_MPI_Fcube                             !InTf!
  type(RPN_MPI_Fcube), intent(INOUT) :: topo           !InTf!
  integer, intent(IN) :: npe                           !InTf!
  integer, intent(IN) :: npex, npey                    !InTf!
  integer, intent(IN) :: blkx, blky                    !InTf!
  logical, intent(IN) :: x_first                       !InTf!
  integer, intent(OUT) :: ierr                         !InTf!

  integer :: me_x, me_y, edge, color, key, fn, edgeflag
  logical :: north, south, east, west
  type(RPN_MPI_Fcom) :: edgecom

  ierr = MPI_ERROR
  if(topo%version .ne. mpi_symbols_version) return     ! RPN_MPI version mismatch
  if(npe .ne. npex*npey*6)                  return     ! bad number of PEs

  ! nullify all edges
  topo%west  = RPN_MPI_Fedge_NULL
  topo%east  = RPN_MPI_Fedge_NULL
  topo%north = RPN_MPI_Fedge_NULL
  topo%south = RPN_MPI_Fedge_NULL

  topo%faceno = topo%cube%size / topo%cube%rank
  fn = topo%faceno

  ! split cube into faces, color = face number, key is global rank in cube
  call MPI_Comm_split(topo%cube%comm, topo%faceno, topo%cube%rank, topo%face%grd%comm, ierr)
  call MPI_Comm_rank(topo%face%grd%comm, topo%face%grd%rank, ierr)
  call MPI_Comm_size(topo%face%grd%comm, topo%face%grd%size, ierr)

  ! split face into rows and columns
  call RPN_MPI_grid_topo(topo%face, npex, npey, blkx, blky, x_first, ierr)

  me_x = topo%face%row%rank
  me_y = topo%face%col%rank
  ! are we on a face edge ?
  west  = ( me_x == 0 )
  east  = ( me_x == npex - 1 )
  north = ( me_y == npey - 1 )
  south = ( me_y == 0 )

  ! now we can start building the edge communicators 
  ! (color will be 1 only for the selected edge)
  do edge = 0, 11
    color = 0    ! NOT AN EDGE
    key = 1000000000 * fn   ! will rank by face number first, then rank on local edge
    edgeflag = EDGE_NONE
    select case(edge)
      case( 0)   ! 0 W, 3 E
        if( (fn == 0 .and. west) .or. (fn == 3 .and.east) ) then
          color = 1
          key   = key + me_x
          edgeflag = edgeflag + EDGE_YY
        endif
      case( 1)   ! 0 E, 1 W
        if( (fn == 0 .and. east) .or. (fn == 1 .and. west) ) then
          color = 1
          key   = key + me_x
          edgeflag = edgeflag + EDGE_YY
        endif
      case( 2)   ! 1 E, 2 W
        if( (fn == 1 .and. east) .or. (fn == 2 .and. west) ) then
          color = 1
          key   = key + me_x
          edgeflag = edgeflag + EDGE_YY
        endif
      case( 3)   ! 2 E, 3 W
        if( (fn == 2 .and. east) .or. (fn == 3 .and. west) ) then
          color = 1
          key   = key + me_x
          edgeflag = edgeflag + EDGE_YY
        endif
      case( 4)   ! 3 N, 4 W, reverse
        if( (fn == 3 .and. north) .or. (fn == 4 .and. west) ) then
          color = 1
          if(fn == 3) key = key + me_y
          if(fn == 4) key = key + me_x
          edgeflag = edgeflag + EDGE_XY
          edgeflag = edgeflag + EDGE_REV
        endif
      case( 5)   ! 0 N, 4 S
        if( (fn == 0 .and. north) .or. (fn == 4 .and. south) ) then
          color = 1
          key = key + me_y
          edgeflag = edgeflag + EDGE_XX
        endif
      case( 6)   ! 1 N, 4 E
        if( (fn == 1 .and. north) .or. (fn == 4 .and. east) ) then
          color = 1
          if(fn == 1) key = key + me_y
          if(fn == 4) key = key + me_x
          edgeflag = edgeflag + EDGE_XY
        endif
      case( 7)   ! 2 N, 4 N, reverse
        if( (fn == 2 .and. north) .or. (fn == 4 .and. north) ) then
          color = 1
          key = key + me_y
          edgeflag = edgeflag + EDGE_REV
          edgeflag = edgeflag + EDGE_XX
        endif
      case( 8)   ! 3 S, 5 W
        if( (fn == 3 .and. south) .or. (fn == 5 .and. west) ) then
          color = 1
          if(fn == 3) key = key + me_y
          if(fn == 5) key = key + me_x
          edgeflag = edgeflag + EDGE_XY
        endif
      case( 9)   ! 0 S, 5 N
        if( (fn == 0 .and. south) .or. (fn == 5 .and. north) ) then
          color = 1
          key = key + me_y
          edgeflag = edgeflag + EDGE_XX
        endif
      case(10)   ! 1 S, 5 E, reverse
        if( (fn == 1 .and. south) .or. (fn == 5 .and. east) ) then
          color = 1
          if(fn == 1) key = key + me_y
          if(fn == 5) key = key + me_x
          edgeflag = edgeflag + EDGE_REV
          edgeflag = edgeflag + EDGE_XY
        endif
      case(11)   ! 2 S, 5 S, reverse
        if( (fn == 2 .and. south) .or. (fn == 5 .and. south) ) then
          color = 1
          key = key + me_y
          edgeflag = edgeflag + EDGE_REV
          edgeflag = edgeflag + EDGE_XX
        endif
    end select
    call MPI_Comm_split(topo%cube%comm, color, key, edgecom%comm, ierr)
    call MPI_Comm_rank(edgecom%comm, edgecom%rank, ierr)
    call MPI_Comm_size(edgecom%comm, edgecom%size, ierr)
    select case(edge)
      case( 0)   ! 0 W, 3 E
        if(fn == 0 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
        if(fn == 3 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
      case( 1)   ! 0 E, 1 W
        if(fn == 0 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
        if(fn == 1 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
      case( 2)   ! 1 E, 2 W
        if(fn == 1 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
        if(fn == 2 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
      case( 3)   ! 2 E, 3 W
        if(fn == 2 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
        if(fn == 3 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npey, npey)
      case( 4)   ! 3 N, 4 W, reverse
        if(fn == 3 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
        if(fn == 4 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
      case( 5)   ! 0 N, 4 S
        if(fn == 0 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
        if(fn == 4 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
      case( 6)   ! 1 N, 4 E
        if(fn == 1 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
        if(fn == 4 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
      case( 7)   ! 2 N, 4 N, reverse
        if(fn == 2 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
        if(fn == 4 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
      case( 8)   ! 3 S, 5 W
        if(fn == 3 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
        if(fn == 5 .and. west ) topo%west  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
      case( 9)   ! 0 S, 5 N
        if(fn == 0 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
        if(fn == 5 .and. north) topo%north = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
      case(10)   ! 1 S, 5 E, reverse
        if(fn == 1 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
        if(fn == 5 .and. east ) topo%east  = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npey)
      case(11)   ! 2 S, 5 S, reverse
        if(fn == 2 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
        if(fn == 5 .and. south) topo%south = RPN_MPI_Fedge(edgecom, edgeflag, edge, npex, npex)
    end select
  enddo

 end subroutine RPN_MPI_cube_topo                     !InTf!

! initialize a CUBE topology to null values
 subroutine RPN_MPI_init_cube(topo, ierr)             !InTf!
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
#include <RPN_MPI.hf>

!! import :: RPN_MPI_Fcube                            !InTf!
  type(RPN_MPI_Fcube), intent(INOUT) :: topo          !InTf!
  integer, intent(OUT) :: ierr                        !InTf!

  ierr = MPI_ERROR
  if(topo%version .ne. mpi_symbols_version) return     ! RPN_MPI version mismatch

  topo = RPN_MPI_Fcube_NULL

  ierr = MPI_SUCCESS
 end subroutine RPN_MPI_init_cube                     !InTf!
