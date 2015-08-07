!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                K H                        \\\\\\\\\\
!
!=======================================================================
!
subroutine kh
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use domain
      use scratch
#ifdef MPI_USED
      use mpiyes
#else 
      use mpino
#endif 
      use mpipar

      implicit none
      integer  :: i,j,k, iseed=-145, seedsize
      integer, dimension(:), allocatable :: seed
      real(rl) :: v0, rho_in, rho_out, ampl, rn, p0, line_val
      real(rl), dimension(:), allocatable :: rannum
      real(rl), dimension(:,:,:), allocatable :: pert
      external ran1
      namelist /pgen/ v0, p0, rho_in, rho_out, ampl, line_val

      v0      = 0.5
      p0      = 2.5
      rho_in  = 2.0
      rho_out = 1.0
      ampl    = 0.01
      line_val= 0.25
      if(myid.eq.0)then
         write(*,*)"K-H problem setup."
      endif
!
!     read in namelist and write it to log file
!
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
        buf_in( 1)  = v0
        buf_in( 2)  = p0
        buf_in( 3)  = rho_in
        buf_in( 4)  = rho_out
        buf_in( 5)  = ampl
        buf_in( 6)  = line_val
      endif
      call MPI_BCAST(buf_in, 6, MPI_FLOAT , 0, comm3d, ierr )
      if (myid .ne. 0) then
        v0      = buf_in(1)
        p0      = buf_in(2)
        rho_in  = buf_in(3)
        rho_out = buf_in(4)
        ampl    = buf_in(5)
        line_val= buf_in(6)
#endif /* MPI_USED */
      endif
!
!     
      allocate(pert(in,jn,kn))
      allocate(rannum(in))
      call random_seed ( SIZE=seedsize)
!
! this stuff lets you set and/or adjust the seed
! multiply the seed by the grid location so that
! you will get a different seed on each processor
! if you are doing a multiprocessor run
      do j=1,jn
        do i=1,in
          pert(i,j,:) = ampl*sin(2.*pi/ABS(x1max-x1min)*x1b(i))
        enddo
      enddo
!
      !allocate(seed(seedsize))
      !call random_seed ( GET=seed(1:seedsize))
      !!seed = seed*agrid%x(1)*agrid%y(1)*agrid%z(1)
      !seed = seed*x1b(1)*x2b(1)*x2b(1)
      !write(*,*) 'seed: ', seed
      !call random_seed ( PUT=seed(1:seedsize))
      !call random_number(rannum)
      !do i=1,in
      !  pert(i,:,:) = ampl*(rannum(i)-0.5)*2.0
      !enddo
      !deallocate(rannum)

      d = tiny
      call ran1(-iseed, rn)
!      do k=ks-2,ke+2
!        do j=js-2,je+2
!          do i=is-2,ie+2
      do k=1,kn
        do j=1,jn
          do i=1,in
            if(abs(x2b(j)) .ge. line_val) then
              d(i,j,k)=rho_out
              v1(i,j,k)=rho_out*(-v0+pert(i,j,k))
            else if (abs(x2b(j)) .lt. line_val) then
              d(i,j,k)=rho_in
              v1(i,j,k)=rho_in*(v0+pert(i,j,k))
            endif
            e(i,j,k)=p0/(gamma -1)
          enddo !i
        enddo !j
      enddo !k
      v2= d*pert
      v3 = 0.0
!
      deallocate(pert)

      if(myid_w .eq. 0) then
        write(*,*) "KH SETUP FINISHED"
        write(*,*) maxval(d), minval(d)
      endif
      return
end subroutine kh
