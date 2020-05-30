!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
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
      SUBROUTINE RPN_COMM_transpose(za,min1,max1,n1g,n2,min3,max3,n3g,zb,type,size)
      use rpn_comm
      implicit none
      integer min1,max1,n1g,n2,min3,max3,n3g,type,size
      integer za(size,min1:max1,n2,n3g)
      integer zb(size,n2,min3:max3,n1g)
!
      integer n3partiel,n1partiel,npe,pecomm
!
      if(abs(type).eq.1) then   ! transpose along X
        npe=pe_nx
        pecomm=pe_myrow
      else                      ! transpose along Y
        npe=pe_ny
        pecomm=pe_mycol
      endif
      n3partiel=(n3g+npe-1)/npe
      n1partiel=(n1g+npe-1)/npe
!
!	check that min1>0, max1>=n1partiel
!	check that min3>0, max3>=n3partiel
!	check that size = 1 or 2 (integer/real, real*8)
!
      if(type.gt.0) then  ! forward transpose
        call RPN_COMM_Xpose1(n3partiel,npe,pecomm,n1partiel, za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
      else ! backward transpose
        call RPN_COMM_Xpose2(n3partiel,npe,pecomm,n1partiel,za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
      endif
      return
      end SUBROUTINE RPN_COMM_transpose
      SUBROUTINE RPN_COMM_Xpose1(n3partiel,npe,pecomm,n1partiel,za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
	use rpn_comm
!	forward transpose, za to zb
!	gather first dimension into processor,
!	distribute last dimension
!	result array has gathered index as last dimension
!	(last dimension of arrays is always in-processor)
!
	implicit none

	integer n3partiel,n1partiel,npe,pecomm,size
	integer min1,max1,n1g,n2,min3,max3,n3g
	integer za(size,min1:max1,n2,n3g)
	integer zb(size,n2,min3:max3,n1g)
	integer *8 za8(size/2,min1:max1,n2,n3g)
	integer *8 zb8(size/2,n2,min3:max3,n1g)
	pointer (za8_,za8)
	pointer (zb8_,zb8)
	integer, allocatable :: ta(:,:,:,:,:)
	integer*8 :: ta8(size/2,n2,min3:max3,n1partiel,npe)
	pointer (ta8_,ta8)
	integer i,j,k,iter,i0,ierr,isz,isz2,n3w
	logical odd
!
	za8_ = loc(za)
	zb8_ = loc(zb)
	isz2 = size/2
	odd = isz2*2 .ne. size
	if(odd) isz2=size
	n3w = max3-min3+1
!
	if(npe.eq.1)then
	  do isz=1,isz2
	  do k=1,n3g
	  do j=1,n2
	    if(odd) then
	      do i=1,n1g
	        zb(isz,j,k,i)=za(isz,i,j,k)
	      enddo
	    else
!VDIR NODEP
	      do i=1,n1g
	        zb8(isz,j,k,i)=za8(isz,i,j,k)
	      enddo
	    endif
	  enddo
	  enddo
	  enddo
	  return
	endif

	allocate(ta(size,n2,n3w,n1partiel,npe))
	ta8_=loc(ta)

	do isz=1,isz2
	i0 = 0
	do k=1,n3g,n3partiel
	 i0 = 1+i0
	 do i=1,n1partiel
	  if(odd) then
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      ta(isz,j,min3,i,i0)=za(isz,i,j,k)
	    enddo
	  else
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      ta8(isz,j,min3,i,i0)=za8(isz,i,j,k)
	    enddo
	  endif
	 enddo
	enddo
	enddo
!
	call MPI_ALLTOALL(ta, size*n1partiel*n2*n3w, MPI_INTEGER, &
                          zb, size*n1partiel*n2*n3w, MPI_INTEGER, pecomm, ierr)
!
	deallocate(ta)
	return
	end
	SUBROUTINE RPN_COMM_Xpose2(n3partiel,npe,pecomm,n1partiel,za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
	use rpn_comm
!
!	backward transpose, zb to za
!	(last dimension of arrays is always in-processor)
!
	implicit none
	integer n3partiel,n1partiel,npe,pecomm,size
	integer min1,max1,n1g,n2,min3,max3,n3g
	integer za(size,min1:max1,n2,n3g)
	integer zb(size,n2,min3:max3,n1g)
!
	integer *8 za8(size/2,min1:max1,n2,n3g)
	integer *8 zb8(size/2,n2,min3:max3,n1g)
	pointer (za8_,za8)
	pointer (zb8_,zb8)
	integer, allocatable :: ta(:,:,:,:,:)
	integer *8, dimension(size/2,n2,min3:max3,n1partiel,npe) :: ta8
	pointer(ta8_,ta8)
	integer i,j,k,iter,i0,ierr,isz,isz2,n3w
	logical odd
!
	za8_ = loc(za)
	zb8_ = loc(zb)
	isz2 = size/2
	odd = isz2*2 .ne. size
	if(odd) isz2=size
	n3w = max3-min3+1
	if(npe.eq.1)then
	  do isz=1,isz2
	  do k=1,n3g
	  do j=1,n2
	    if(odd) then
	      do i=1,n1g
	        za(isz,i,j,k)=zb(isz,j,k,i)
	      enddo
	    else
!VDIR NODEP
	      do i=1,n1g
	        za8(isz,i,j,k)=zb8(isz,j,k,i)
	      enddo
	    endif
	  enddo
	  enddo
	  enddo
	  return
	endif
	allocate(ta(size,n2,min3:max3,n1partiel,npe))
	ta8_=loc(ta)
	call MPI_ALLTOALL(zb, size*n1partiel*n2*n3w, MPI_INTEGER, &
                          ta, size*n1partiel*n2*n3w, MPI_INTEGER, pecomm, ierr)
	do isz=1,isz2
	i0 = 0
	do k=1,n3g,n3partiel
	 i0 = 1+i0
	 do i=1,n1partiel
	  if(odd) then
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      za(isz,i,j,k)=ta(isz,j,min3,i,i0)
	    enddo
	  else
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      za8(isz,i,j,k)=ta8(isz,j,min3,i,i0)
	    enddo
	  endif
	 enddo
	enddo
	enddo
	deallocate(ta)
	return
	end
