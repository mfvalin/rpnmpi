#if defined(STAND_ALONE)
program test_106
call rpn_mpi_test_106
end
#endif

module mod_106
  implicit none
  save
  type :: edge
    real :: x, y, z
  end type
  character(len=64) fmtx
  character(len=64) fmty
  character(len=64) fmtz

end module mod_106

subroutine rpn_mpi_test_106
  use ISO_C_BINDING
  use mod_106
  implicit none
#include <RPN_MPI.hf>
  integer :: ier, npe, rank, larg
  type(RPN_MPI_mpi_definitions_raw) :: dr
  integer :: i, ln, side, stat, npex, npey, myside, siderank, li, lj, pe_i, pe_j
  character(len=128) :: argv(4)

  call MPI_Init(ier)

  call RPN_MPI_get_mpi_definitions_raw(dr, ier)       ! get "raw" definitions

  call MPI_Comm_size(dr%MPI_COMM_WORLD, npe, ier)     ! get world size
  call MPI_Comm_rank(dr%MPI_COMM_WORLD, rank, ier)    ! get world rank

  write(6,*) 'I am PE',rank+1,' of',npe

  do i=1,4
    call get_command_argument(i,argv(i),larg,stat)    ! get 3 program arguments
    if(stat .ne. 0) goto 777
  enddo
  read(argv(1),*,err=777) npex                        ! PE topology (npex, npey)
  read(argv(2),*,err=777) npey
!   read(argv(3),*,err=777) ln                          ! cube edge length
  read(argv(3),*,err=777) li                          ! PE edge length
  read(argv(4),*,err=777) lj                          ! PE edge length
  myside = rank / (npex * npey)
  siderank = rank - myside * npex * npey
  write(6,*) 'myside =',myside,', siderank =',siderank

  if(npex*npey*6 .ne. npe .and. npe .ne. 1) goto 888  ! npex x npey PEs per side, 6 sides, if a cube
  ln = li * npex
  if(ln .ne. lj * npey) goto 888


  if(npe .eq. 1) then
    write(fmtx,1) '(A,',ln,'F7.3,A,F7.3,A,F7.3)'
    write(fmty,1) '(A,F7.3,A,',ln,'F7.3,A,F7.3)'
    write(fmtz,1) '(A,F7.3,A,F7.3,A,',ln,'F7.3,A)'
    do side = 0 , 5
      call test_edges(side, ln)
    enddo
  else
    if(siderank == 0) then
      write(fmtx,1) '(A,',ln,'F7.3,A,F7.3,A,F7.3)'
      write(fmty,1) '(A,F7.3,A,',ln,'F7.3,A,F7.3)'
      write(fmtz,1) '(A,F7.3,A,F7.3,A,',ln,'F7.3,A)'
      call test_edges(myside, ln)
    endif
    if(npex * npey > 1) then      !  sides are decomposed
      write(fmtx,1) '(A,',li,'F7.3,A,F7.3,A,F7.3)'
      write(fmty,1) '(A,F7.3,A,',li,'F7.3,A,F7.3)'
      write(fmtz,1) '(A,F7.3,A,F7.3,A,',li,'F7.3,A)'
      pe_j = siderank / npex
      pe_i = mod(siderank, npex)
      write(6,*) 'I am sub PE',pe_i, pe_j
      call test_all_edges(myside, ln, li, lj, pe_i, pe_j, npex, npey)
    endif
  endif

7 continue
  call MPI_Finalize(ier)
  stop

777 continue
  write(6,*) 'MISSING or ill formatted arguments'
  goto 7

888 continue
  write(6,*) 'WRONG values in arguments'
  goto 7

1 format(A,I0,A)
  return
end subroutine rpn_mpi_test_106

subroutine test_edges(side, ln)  ! fill North, South, East, and West edges of a given side of the cube
  use ISO_C_BINDING
  use mod_106
  implicit none
  integer, intent(IN) :: side    ! cube side number (0-5)
  integer, intent(IN) :: ln      ! side dimension is (ln x ln), edge dimension is (ln)
  type(edge), dimension(ln) :: n, s, e, w  ! the 4 edges (North, South, East, and West)

  call set_edge_values(side, ln, n, s, w, e) ! cube has a range of -1.0 to 1.0 alng x, y, z
  print 2,'side =',side
  !                                                       4
  ! there are 6 sides to the cube, labelled 0 to 5     3  0  1  2
  !                                                       5
  select case(side)
  case(0)  ! x = 1.0  (front side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(1)  ! y = 1.0  (right side)
    print fmtx,'North x = [',n(:)%x,'], y = ',n(1)%y,', z =',n(1)%z
    print fmtx,'South x = [',s(:)%x,'], y = ',s(1)%y,', z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(2)  ! x = -1.0  (back side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(3)  ! y = -1.0  (left side)
    print fmtx,'North x = [',n(:)%x,'], y = ',n(1)%y,', z =',n(1)%z
    print fmtx,'South x = [',s(:)%x,'], y = ',s(1)%y,', z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(4)  ! z = 1.0  (top side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtx,'East  x = [',e(:)%x,'], y = ',e(1)%y,', z =',e(1)%z
    print fmtx,'West  x = [',w(:)%x,'], y = ',w(1)%y,', z =',w(1)%z
  case(5)  ! z = -1.0  (bottom side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtx,'East  x = [',e(:)%x,'], y = ',e(1)%y,', z =',e(1)%z
    print fmtx,'West  x = [',w(:)%x,'], y = ',w(1)%y,', z =',w(1)%z
  end select
1 format(A,F7.3,A,20F7.3)
2 format(A,I2,A,F7.3)
  return
end subroutine test_edges

subroutine test_all_edges(side, ln, li, lj, pe_i, pe_j, npe_i, npe_j)  ! fill North, South, East, and West edges
  use ISO_C_BINDING
  use mod_106
  implicit none
  integer, intent(IN) :: side    ! cube side number (0-5)
  integer, intent(IN) :: ln      ! side dimension is (ln x ln), edge dimension is (ln)
  integer, intent(IN) :: li, lj   ! side tile dimensions (li , lj)
  integer, intent(IN) :: pe_i, pe_j   ! tile coordinates on a side ( 0 : npe_i-1 , 0 : npe_j-1 )
  integer, intent(IN) :: npe_i, npe_j ! number of PEs for a cube side
  type(edge), dimension(li) :: n, s  ! North and South edges
  type(edge), dimension(lj) :: w, e  ! West and East edges

  call set_edges(side, ln, li, lj, pe_i, pe_j, npe_i, npe_j, n, s, w, e) ! cube has a range of -1.0 to 1.0 alng x, y, z
  print 2,'side =',side
  !                                                       4
  ! there are 6 sides to the cube, labelled 0 to 5     3  0  1  2
  !                                                       5
  select case(side)
  case(0)  ! x = 1.0  (front side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(1)  ! y = 1.0  (right side)
    print fmtx,'North x = [',n(:)%x,'], y = ',n(1)%y,', z =',n(1)%z
    print fmtx,'South x = [',s(:)%x,'], y = ',s(1)%y,', z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(2)  ! x = -1.0  (back side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(3)  ! y = -1.0  (left side)
    print fmtx,'North x = [',n(:)%x,'], y = ',n(1)%y,', z =',n(1)%z
    print fmtx,'South x = [',s(:)%x,'], y = ',s(1)%y,', z =',s(1)%z
    print fmtz,'East  x = ',e(1)%x,', y = ',e(1)%y,', z = [',e(:)%z,']'
    print fmtz,'West  x = ',w(1)%x,', y = ',w(1)%y,', z = [',w(:)%z,']'
  case(4)  ! z = 1.0  (top side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtx,'East  x = [',e(:)%x,'], y = ',e(1)%y,', z =',e(1)%z
    print fmtx,'West  x = [',w(:)%x,'], y = ',w(1)%y,', z =',w(1)%z
  case(5)  ! z = -1.0  (bottom side)
    print fmty,'North x = ',n(1)%x,', y = [',n(:)%y,'], z =',n(1)%z
    print fmty,'South x = ',s(1)%x,', y = [',s(:)%y,'], z =',s(1)%z
    print fmtx,'East  x = [',e(:)%x,'], y = ',e(1)%y,', z =',e(1)%z
    print fmtx,'West  x = [',w(:)%x,'], y = ',w(1)%y,', z =',w(1)%z
  end select
1 format(A,F7.3,A,20F7.3)
2 format(A,I2,A,F7.3)
  return
end subroutine test_all_edges

subroutine set_edges(side, ln, li, lj, pe_i, pe_j, npe_i, npe_j, n, s, w, e)
  use ISO_C_BINDING
  use mod_106
  implicit none
  integer, intent(IN) :: side     ! cube side number (0-5)
  integer, intent(IN) :: ln       ! cube side dimensions (ln , ln)
  integer, intent(IN) :: li, lj   ! side tile dimensions (li , lj)
  integer, intent(IN) :: pe_i, pe_j   ! tile coordinates on a side ( 0 : npe_i-1 , 0 : npe_j-1 )
  integer, intent(IN) :: npe_i, npe_j ! number of PEs for a cube side
  type(edge), intent(OUT), dimension(li) :: n, s  ! North and South edges
  type(edge), intent(OUT), dimension(lj) :: w, e  ! West and East edges

  integer :: i
  real, dimension(ln) :: x
  real :: deltai, deltaj, delta, basei, basej, topi, topj

  delta  = 1.0 / ln       ! half spacing between points along edges ( ln points : -1 -> +1 )
  deltai = 2.0 / npe_i    ! tile spacing along i
  basei  = pe_i * deltai  ! start of tile along i
  topi   = basei + deltai ! end of tile along i
  deltaj = 2.0 / npe_j    ! tile spacing along j
  basej  = pe_j * deltaj  ! start of tile along i
  topj   = basej + deltaj ! end of tile along j
  do i = 1, max(li,lj)
    x(i) = (2*i -1)
  enddo
                     !                                                       4
  select case(side)  ! there are 6 sides to the cube, labelled 0 to 5     3  0  1  2
                     !                                                       5
  case(0)  ! x = 1.0, I axis = +y, J axis = +z  (front side)
    n%x = 1.0
    n%y = -1.0 + basei + delta * x(1:li)
    n%z = -1.0 + topj
    s%x = 1.0
    s%y = -1.0 + basei + delta * x(1:li)
    s%z = -1.0 + basej
    e%x = 1.0
    e%y = -1.0 + topi
    e%z = -1.0 + basej + delta * x(1:lj)
    w%x = 1.0
    w%y = -1.0 + basei
    w%z = -1.0 + basej + delta * x(1:lj)
  case(1)  ! y = 1.0, I axis = -x, J axis = +z  (right side)
    n%x = 1.0 - basei - delta * x(1:li)
    n%y = 1.0
    n%z = -1.0 + topj
    s%x = 1.0 - basei - delta * x(1:li)
    s%y = 1.0
    s%z = -1.0 + basej
    e%x = 1.0 - topi
    e%y = 1.0
    e%z = -1.0 + basej + delta * x(1:lj)
    w%x = 1.0 - basei
    w%y = 1.0
    w%z = -1.0 + basej + delta * x(1:lj)
  case(2)  ! x = -1.0, I axis = -y, J axis = +z  (back side)
    n%x = -1.0
    n%y = 1.0 - basei - delta * x(1:li)
    n%z = -1.0 + topj
    s%x = -1.0
    s%y = 1.0 - basei - delta * x(1:li)
    s%z = -1.0 + basej
    e%x = -1.0
    e%y = 1.0 - topi
    e%z = -1.0 + basej + delta * x(1:lj)
    w%x = -1.0
    w%y = 1.0 - basei
    w%z = -1.0 + basej + delta * x(1:lj)
  case(3)  ! y = -1.0, I axis = +x, J axis = +z  (left side)
    n%x = -1.0 + basei + deltai * x(1:li)
    n%y = -1.0
    n%z = -1.0 + topj
    s%x = -1.0 + basei + deltai * x(1:li)
    s%y = -1.0
    s%z = -1.0 + basej
    e%x = -1.0 + topi
    e%y = -1.0
    e%z = -1.0 + basej + delta * x(1:lj)
    w%x = -1.0 + basei
    w%y = -1.0
    w%z = -1.0 + basej + delta * x(1:lj)
  case(4)  ! z = 1.0, I axis = +y, J axis = -x  (top side)
    n%x = 1.0 - topj
    n%y = -1.0 + basei + delta * x(1:li)
    n%z = 1.0
    s%x = 1.0 - basej
    s%y = -1.0 + basei + delta * x(1:li)
    s%z = 1.0
    e%x = 1.0 - basej - delta * x(1:lj)
    e%y = -1.0 + topi
    e%z = 1.0
    w%x = 1.0 - basej - delta * x(1:lj)
    w%y = -1.0 + basei
    w%z = 1.0
  case(5)  ! z = -1.0, I axis = +y, J axis = +x  (bottom side)
    n%x = -1.0 + topj
    n%y = -1.0 + basei + delta * x(1:li)
    n%z = -1.0
    s%x = -1.0 + basej
    s%y = -1.0 + basei + delta * x(1:li)
    s%z = -1.0
    e%x = -1.0 + basej + delta * x(1:lj)
    e%y = -1.0 + topi
    e%z = -1.0
    w%x = -1.0 + basej + delta * x(1:lj)
    w%y = -1.0 + basei
    w%z = -1.0
  end select

  return
end subroutine set_edges

subroutine set_edge_values(side, ln, n, s, w, e)
  use ISO_C_BINDING
  use mod_106
  implicit none
  integer, intent(IN) :: side    ! cube side (0-5)
  integer, intent(IN) :: ln      ! side dimension (ln x ln)
  type(edge), intent(OUT), dimension(ln) :: n, s, w, e  ! the 4 edges

  integer :: i
  real, dimension(ln) :: x
  real :: delta

  delta = 1.0 / ln
  do i = 1, ln
    x(i) = -1.0 + (2*i -1) * delta  ! -1 + delta , .... , 1 - delta
  enddo
                     !                                                       4
  select case(side)  ! there are 6 sides to the cube, labelled 0 to 5     3  0  1  2
                     !                                                       5
  case(0)  ! x = 1.0  (front side)
    n%x = 1.0
    n%y = x
    n%z = 1.0
    s%x = 1.0
    s%y = x
    s%z = -1.0
    e%x = 1.0
    e%y = 1.0
    e%z = x
    w%x = 1.0
    w%y = -1.0
    w%z = x
  case(1)  ! y = 1.0  (right side)
    n%x = - x
    n%y = 1.0
    n%z = 1.0
    s%x = - x
    s%y = 1.0
    s%z = -1.0
    e%x = -1.0
    e%y = 1.0
    e%z = x
    w%x = 1.0
    w%y = 1.0
    w%z = x
  case(2)  ! x = -1.0  (back side)
    n%x = -1.0
    n%y = -x
    n%z = 1.0
    s%x = -1.0
    s%y = -x
    s%z = -1.0
    e%x = -1.0
    e%y = -1.0
    e%z = x
    w%x = -1.0
    w%y = 1.0
    w%z = x
  case(3)  ! y = -1.0  (left side)
    n%x = x
    n%y = -1.0
    n%z = 1.0
    s%x = x
    s%y = -1.0
    s%z = -1.0
    e%x = 1.0
    e%y = -1.0
    e%z = x
    w%x = -1.0
    w%y = -1.0
    w%z = x
  case(4)  ! z = 1.0  (top side)
    n%x = -1.0
    n%y = x
    n%z = 1.0
    s%x = 1.0
    s%y = x
    s%z = 1.0
    e%x = -x
    e%y = 1.0
    e%z = 1.0
    w%x = -x
    w%y = -1.0
    w%z = 1.0
  case(5)  ! z = -1.0  (bottom side)
    n%x = 1.0
    n%y = x
    n%z = -1.0
    s%x = -1.0
    s%y = x
    s%z = -1.0
    e%x = x
    e%y = 1.0
    e%z = -1.0
    w%x = x
    w%y = -1.0
    w%z = -1.0
  end select

  return
end subroutine set_edge_values
