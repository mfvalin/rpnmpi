!
!    Test verifying that rpn_comm_transpose is giving correct results.
!
        subroutine get_n_domains(ndomains, offset, err)
        integer n_domains
        common/ndomains/n_domains
        integer ndomains, offset, err
        character (len=128) SYSTEM_COMMAND
	SYSTEM_COMMAND="1"
        call get_environment_variable("TEST_DOMAINS",SYSTEM_COMMAND)
        if(SYSTEM_COMMAND == "" )SYSTEM_COMMAND="1"
        read(SYSTEM_COMMAND,*)ndomains
        n_domains=ndomains
        offset = 0
        err = 0
        return
        end
        subroutine rpn_mpi_test_102
        implicit none
!
!       dummy test program implemented as a subroutne
!       because some compilers do not like POINTER
!       statements with variable dimensions in a main program
!
        integer :: nptsx,nptsy,nptsz,ihalox,ihaloy
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy  ! common containing problem dimensions and halo size
        integer :: nodex, nodey
        common /pernode/ nodex, nodey

        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
        integer, allocatable, dimension(:)  ::iarr,jarr
        integer, allocatable, dimension(:,:,:) ::data,data2,glob
!
        external sss
        integer npex,npey
        integer min3,max3,ierr, istat
        integer mini,maxi,nil,nilmax,ni0,i0
        integer minj,maxj,njl,njlmax,nj0,j0
        integer nzl, nzlmzx, nz0,irep,iter,irep2
        real*8 time1,time2,MPI_wtime,time_min,time_max,time_tot
        external MPI_wtime
!
        integer RPN_COMM_topo
        integer RPN_COMM_init_multi_level, mygrid, mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm
        integer RPN_COMM_colors
        integer RPN_COMM_comm
        integer test_grids
	integer RPN_COMM_option_L
        external RPN_COMM_mydomain
        external get_n_domains
        external test_grids
	external RPN_COMM_option_L
        integer domains, mydomain,irank,isize
        integer n_domains,j
        integer peercomm, npeers
        common/ndomains/n_domains
        character *6 is_async(2)

        integer cache_flush(2000*2000)
        integer, parameter :: NEXCH=256
        real, parameter :: SCAL=1.0/NEXCH
        real time_ew(0:5*NEXCH), time_ns(0:5*NEXCH)
        integer nbytes, nodebytes
!
!       force a "callback" type initialization
!
        npex=0
        npey=0
        call RPN_COMM_mydomain (get_n_domains, mydomain,ierr)
!
        call RPN_COMM_set_petopo(1,1000)   ! force vertically striped distribution
        mygrid = RPN_COMM_init_multi_level(sss,Pelocal,Petotal,npex,npey,n_domains,1)
!	============= TEST for WORLD/MULTIGRID/GRID ====================
!
        if(test_grids(mygrid,Pelocal) .ne. 0) goto 9999
!
!	================= TEST for transpose ===================
!
           call setup_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl)
           istat =  RPN_COMM_topo(nptsz,min3,max3,nzl,nzlmzx,0,nz0,.true.,.true.)
           allocate(data2((max3-min3+1),(maxj-minj+1),(nptsx+10)*2))
           if(.true.) then
             if(Pelocal.eq.0 )then
              print *,' size of data =',                                      &
     &             (maxi-mini+1),(maxj-minj+1),nptsz,                         &
     &             (maxi-mini+1)*(maxj-minj+1)*nptsz*2
              print *,' size of data2 =',                                     &
     &             (maxj-minj+1),(max3-min3+1),nptsx,                         &
     &             (max3-min3+1)*(maxj-minj+1)*(nptsx+10)*2
              print *,'nptsz,min3,max3,nzl,nzlmzx,nz0',                       &
     &             nptsz,min3,max3,nzl,nzlmzx,nz0
             endif
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
           endif
!
!     THE END !!
!
!
 2222   Continue
        call RPN_COMM_BARRIER('WORLD',ierr)
 9999	continue
        call RPN_COMM_FINALIZE(ierr)
        stop
        end
!
        subroutine vfy_xpos(z,jtab,minj,maxj,njl,nz0,nzl,min3,max3,nptsx,j0)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers
!
!
        integer minj,maxj,njl,min3,max3,nptsx,nzl,nz0,j0
        integer z(2,minj:maxj,min3:max3,nptsx)
        integer jtab(minj:maxj)
        integer i,j,k
        integer error
        error=0
        do i=1,nptsx
        do k=1,nzl
        do j=1,njl
            if(z(1,j,k,i)/100000 .ne. i) error=error+1
            if(mod(z(1,j,k,i), 100000)/100 .ne. j0+j-1) error=error+1
            if(z(2,j,k,i).ne.( nz0-1+k )) error = error + 1
        enddo
        enddo
        enddo
        if(error.ne.0) then 
          print *,error,' VERIFICATION ERRORS for Pe ',Pelocal
        endif
      if(Pelocal.eq.0) print *,'vfy_xpos : Verification performed, number of errors was ',error
        return
        end
!
        subroutine vfy_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers ignoring k
!
        integer minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy
        integer z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k,k0
        integer error
        error=0
        do k=1,nk
         k0=mod(k,3)
          do j=1,njl
          do i=1,nil
            if(z(1,i,j,k).ne.jtab(j)*100+itab(i)*100000) error = error + 1
            if(z(2,i,j,k).ne.( k ))  error = error + 1
          enddo
          enddo
        enddo
        if(error.ne.0) then 
          print 1,'vfy_arr2:',error,' VERIFICATION ERRORS for Pe ',Pelocal,' (',(nil)*(njl)*nk,' points)'
1       format(A,I4,A,I4,A,I5,A)
        endif
        if(Pelocal.eq.0) print *,'vfy_arr2 : Verification performed,',error,' errors on PE 0'
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine set_ijarr(itab,jtab,minx,maxx,miny,maxy,i0,j0,nptsx,nptsy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        precompute initial position tables
!
        integer minx,maxx,miny,maxy,nptsx,nptsy,i0,j0
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j
        do i=minx,maxx
          itab(i)=i0+i-1
          if(itab(i).gt.nptsx) itab(i)=itab(i)-nptsx
          if(itab(i).le.0    ) itab(i)=itab(i)+nptsx
        enddo
        do j=miny,maxy
          jtab(j)=j0+j-1
          if(jtab(j).gt.nptsy) jtab(j)=jtab(j)-nptsy
          if(jtab(j).le.0    ) jtab(j)=jtab(j)+nptsy
        enddo
        return
        end
!
        subroutine setup_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl)
!!!!        use rpn_comm
        implicit none 
!
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!        fill array with markers ignoring k subscript
!
        integer minx,maxx,miny,maxy,nk
        integer nil,njl,z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k
        z = -1
        do k=1,nk
          do j=1,njl
          do i=1,nil
              z(1,i,j,k) = jtab(j)*100
              z(1,i,j,k) = z(1,i,j,k) + itab(i)*100000
              z(2,i,j,k) = k
          enddo
          enddo
        enddo
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        SUBROUTINE sss(nx,ny)
        implicit none 
        integer zouf,nx,ny
!
!        "callback routine" used to get initial topology
!        information
!
        integer :: nodex, nodey
        common /pernode/ nodex, nodey
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy
        integer nptsx,nptsy,nptsz,ihalox,ihaloy
        integer deltai,deltaj
        open(5,file='TEST_data_001',form='FORMATTED')
        print *,'PEs =',nx*ny
        read(5,*)nx,ny,nptsx,nptsy,nptsz,ihalox,ihaloy,deltai,deltaj,nodex,nodey
        print *, ' problem size is ',nptsx,' by ',nptsy,' by ',nptsz
        print *, ' halo size is ',ihalox,' by ',ihaloy
        print *, ' topology = ',nx,' by ',ny
        print *, ' PE block topology = ',deltai,' by ',deltaj
        print *, 'Node tiles = ',nodex,' by ',nodey
        call RPN_COMM_set_petopo(deltai,deltaj)
        return
        end
!
	integer function test_grids(mygrid,Pelocal)
	implicit none
	integer mygrid, Pelocal

	integer RPN_COMM_comm, RPN_COMM_colors
	external RPN_COMM_comm, RPN_COMM_colors
	integer mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm,peercomm
	integer irank, isize, ierr, npeers

        test_grids = -1
        mymultigridcomm = RPN_COMM_comm('MULTIGRID')
        myallgridcomm   = RPN_COMM_comm('ALLGRIDS')
        mygridcomm      = RPN_COMM_comm('GRID')
        peercomm        = RPN_COMM_comm('GRIDPEERS')
        myworldcomm     = RPN_COMM_comm('WORLD')
!
        call RPN_COMM_BARRIER('GRID',ierr)
        call RPN_COMM_rank( 'GRID', irank ,ierr )
        call RPN_COMM_size( 'GRID', isize ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRID',mygridcomm,isize

        call RPN_COMM_BARRIER('MULTIGRID',ierr)
        call RPN_COMM_rank( 'MULTIGRID', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on MULTIGRID',mymultigridcomm

        call RPN_COMM_BARRIER('ALLGRIDS',ierr)
        call RPN_COMM_rank( 'ALLGRIDS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on ALLGRIDS',myallgridcomm

        call RPN_COMM_BARRIER('WORLD',ierr)
        call RPN_COMM_rank( 'WORLD', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on WORLD',myworldcomm

        call RPN_COMM_BARRIER('GRIDPEERS',ierr)
        call RPN_COMM_rank( 'GRIDPEERS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRIDPEERS',peercomm

        call RPN_COMM_size( 'GRIDPEERS', npeers ,ierr )
        if(Pelocal.eq.0 )then
          print *,'mygrid=',mygrid,'myworld',myworldcomm,'mymultigridcomm=',mymultigridcomm,'mygridcomm=',mygridcomm
          print *,'ID=',RPN_COMM_colors('WORLD'),'/',RPN_COMM_colors('MULTIGRID'),'/',RPN_COMM_colors('GRID')
          print *,'peer to peer comm=',peercomm,' npeers=',npeers
        endif
        call RPN_COMM_BARRIER('WORLD',ierr)
        test_grids = 0
	return
	end
        subroutine timing_report(n,times,label)
        implicit none
        integer :: n
        character (len=*) :: label
        real, dimension(3,n) :: times
        real tmin(3), tmax(3)
        real *8 ttot(3), ttot2(3)
        integer i,j
        character (len=10) :: labels(3)

        labels(1)=' PRE COPY '
        labels(2)=' MPI CALL '
        labels(3)=' POST COPY'

!            print *,label,' N=',n
!            print 111,'T1 ',nint(times(1,:)*1000000)
!            print 111,'T2 ',nint(times(2,:)*1000000 - times(1,:)*1000000)
!            print 111,'T3 ',nint(times(3,:)*1000000 - times(2,:)*1000000)
!111         format(A,15I6)
!        return
        times = times * 1000000
        ttot = 0.0
        tmin = 1000000000.0
        tmax = 0.0
        do j = 1 , n
         times(3,j) = times(3,j) - times(2,j)
         times(2,j) = times(2,j) - times(1,j)
         do i = 1 , 3
          tmin(i)=min(tmin(i),times(i,j))
          tmax(i)=max(tmax(i),times(i,j))
          ttot(i)=ttot(i)+times(i,j)
         enddo
        enddo
        ttot = ttot / n
        do I = 1 , 3
          print 222,label//labels(i)//' (min,max,avg) microsec',nint(tmin(i)),nint(tmax(i)),nint(ttot(i))
        enddo
222     format(A,3I8)
        return
        end
