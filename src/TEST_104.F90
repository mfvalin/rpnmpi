subroutine rpn_mpi_test_104
!
! simple, standalone test for RPN_MPI_ez_halo_parms
! does not use the full RPN_MPI setup
!
  use ISO_C_BINDING
  implicit none
  include 'RPN_MPI.inc'
!   include 'RPN_MPI_mpif.inc'

  integer :: ier, npe, rank
  type(RPN_MPI_mpi_definitions_raw) :: dr
  integer :: i, j, k, n, larg, stat, lni, lnj, lnk, gni, gnj, gnk, npex, npey, ranki, rankj
  integer :: i0, j0, k0, i0y, lniy, lnimaxy
  character(len=128) :: argv(5), mode
  integer :: lnimax, lnjmax, lnkmax, errors
  integer(kind=8), dimension(:,:,:), allocatable   :: z, z2        ! original and restored original
  integer(kind=8), dimension(:), allocatable       :: zgi, zgj     ! full row, full column
  integer(kind=8), dimension(:,:,:,:), allocatable :: zt, zy, zty
  integer :: row_comm, col_comm
  integer :: myrow, mycol

  call MPI_Init(ier)

  call RPN_MPI_get_mpi_definitions_raw(dr, ier)  ! get "raw" definitions
  call MPI_Comm_size(dr%MPI_COMM_WORLD, npe, ier)
  call MPI_Comm_rank(dr%MPI_COMM_WORLD, rank, ier)

  do i=1,5
    call get_command_argument(i,argv(i),larg,stat)
    if(stat .ne. 0) goto 777
  enddo
  write(0,*) 'I am PE',rank+1,' of',npe
  read(argv(1),*,err=777) gni
  read(argv(2),*,err=777) gnj
  read(argv(3),*,err=777) gnk
  read(argv(4),*,err=777) npex
  read(argv(5),*,err=777) npey
  if(npex*npey .ne. npe) then
    write(0,*) 'got',npe,' PEs, expected',npex*npey
    goto 777
  endif

  ranki = mod(rank,npex)
  mycol = ranki
  call MPI_Comm_split(dr%MPI_COMM_WORLD, mycol, rank, col_comm, ier)
  rankj = rank / npex
  myrow = rankj
  call MPI_Comm_split(dr%MPI_COMM_WORLD, myrow, rank, row_comm, ier)
  lnimax = (gni + npex -1) / npex      ! gni distributed over x
  i0 = lnimax * ranki
  lnimaxy = (gni + npey -1) / npey     ! gni distributed over y
  i0y = lnimaxy * rankj
  lnjmax = (gnj + npey -1) / npey      ! gnj distributed over y
  j0 = lnjmax * rankj
  lnkmax = (gnk + npex -1) / npex      ! gnk distributed over x
  k0 = lnkmax * ranki
  lni = min(lnimax, gni - (lnimax * ranki))
  lni = max(0,lni)
  lniy = min(lnimaxy, gni - (lnimaxy * ranki))
  lniy = max(0,lniy)
  lnj = min(lnjmax, gnj - (lnjmax * rankj))
  lnj = max(0,lnj)
  lnk = min(lnkmax, gnk - (lnkmax * ranki))
  lnk = max(0,lnk)
  write(0,*) 'gni, gnj, gnk, lni, lnj, lnk =',gni, gnj, gnk, lni, lnj, lnk
  write(0,*) 'ranki, rankj, lnimax, lnimaxy =',ranki, rankj, lnimax, lnimaxy

  allocate( z(lnimax,lnjmax,gnk))                 ! original array
  allocate(z2(lnimax,lnjmax,gnk))                 ! restored original array
  allocate(zt(lnimax,lnjmax,max(1,lnk),npex))     ! after x transpose
  allocate( zgi(lnimax*npex))                         ! a full row, used to reshape zt into zy, and zy into zt
  allocate( zy(lnimaxy,lnjmax,max(1,lnk),npey))   ! before y transpose, after inverse y transpose
  allocate(zty(lnimaxy,lnjmax,max(1,lnk),npey))   ! after y transpose
  allocate(zgj(npey*lnjmax))                      ! a full column

  do k = 1, gnk
  do j = 1, lnj
  do i = 1, lni
    z(i,j,k) = (i0+i)*1000000 + (j0+j)*1000 + k
  enddo
  do i = lni+1, lnimax
    z(i,j,k) = 0
  enddo
  enddo
  enddo
  if(gni <= 16) then
    write(0,*) 'before x transpose'
    do j = lnj,1,-1
      write(0,2)(z(1:lnimax,j,k),k=1,gnk)
    enddo
  endif
  call RPN_MPI_fast_transpose_x(rpn_mpi_loc(loc(z)), rpn_mpi_loc(loc(zt)), .true., lnimax*2, lnjmax, gnk, lnk, row_comm, ier)
  if(gni <= 16) then
    write(0,*) 'after x transpose'
    do j = lnj,1,-1
      write(0,2)(zt(1:lnimax,j,1,n), n=1,npex)
    enddo
  endif

  ! reshape (lnix,lnjy,lnkx,npex) -> (lniy,lnjy,lnkx,npey)
  do k = 1,lnk
  do j = 1,lnjmax
    do n = 0,npex-1    ! (1:lnix,i,j,1:npex) -> (gni)
    do i = 1,lnimax
      zgi (n*lnimax+i) = zt(i,j,k,n+1)   ! gather full row
    enddo
    enddo
    do n = 0,npey-1    ! (gni) -> (1:lniy,i,j,1:npey) 
    do i = 1,lnimaxy
      zy(i,j,k,n+1) = zgi (n*lnimaxy+i)  ! scatter full row
    enddo
    enddo
  enddo
  enddo

  if(gni <= 16) then
    write(0,*) 'before y transpose'
    do j = lnj,1,-1
    enddo
  endif
  call RPN_MPI_fast_transpose_y(rpn_mpi_loc(loc(zy)), rpn_mpi_loc(loc(zty)), .true., lnimaxy*2, lnjmax, lnk, col_comm, ier)
  if(gni <= 16) then
    write(0,*) 'after y transpose'
    do j = lnj,1,-1
    enddo
  endif
  call RPN_MPI_fast_transpose_y(rpn_mpi_loc(loc(zy)), rpn_mpi_loc(loc(zty)), .false., lnimaxy*2, lnjmax, lnk, col_comm, ier)
  if(gni <= 16) then
    write(0,*) 'after inverse y transpose'
    do j = lnj,1,-1
    enddo
  endif

  ! reshape (lniy,lnjy,lnkx,npey) -> (lnix,lnjy,lnkx,npex)
  do k = 1,lnk
  do j = 1,lnjmax
    do n = 0,npey-1    ! (1:lniy,j,k,1:npey) -> (gni)
    do i = 1,lnimaxy
      zgi (n*lnimaxy+i) = zy(i,j,k,n+1)   ! gather full row
    enddo
    enddo
    do n = 0,npex-1    ! (gni) -> (1:lnix,j,k,1:npex) 
    do i = 1,lnimax
      zt(i,j,k,n+1) = zgi (n*lnimax+i)    ! scatter full row
    enddo
    enddo
  enddo
  enddo

  z2 = 0
  call RPN_MPI_fast_transpose_x(rpn_mpi_loc(loc(z2)), rpn_mpi_loc(loc(zt)), .false., lnimax*2, lnjmax, gnk, lnk, row_comm, ier)
  if(gni <= 16) then
    write(0,*) 'after inverse x transpose'
    do j = lnj,1,-1
      write(0,2)(z2(1:lnimax,j,k),k=1,gnk)
    enddo
  endif

  errors = 0
  do k = 1, gnk
  do j = 1, lnj
  do i = 1, lni
    if(z(i,j,k) .ne. z2(i,j,k)) errors = errors + 1
  enddo
  enddo
  enddo
  write(0,*) 'errors =',errors

1 call MPI_Finalize(ier)
  return

777 continue
  write(0,*) 'ERROR in arguments'
  goto 1
2 format(30I10.9)
end subroutine rpn_mpi_test_104
