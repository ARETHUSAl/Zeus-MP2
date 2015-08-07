!
! JReppin 04/08/15
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////            T U R B U L E N C E            \\\\\\\\\\
!=======================================================================
!
subroutine turb_simple
!
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use subs, only: backskip_to_time
      implicit none
      integer :: i,j,k
      integer :: varfile, irho, iux, iuy, iuz, ibx, iby, ibz
      integer :: mxgrid, mygrid, mzgrid, nv
      real(rl):: dx, dy, dz, tvar
      real(rl), dimension (:,:,:,:), allocatable :: ga
      real(rl), dimension (:,:),     allocatable :: buffer
      real(rl), dimension (:),       allocatable :: x, y, z
      integer, dimension (:),        allocatable :: input_vals
      integer, parameter :: tag_ga=675, lun_input=111, ip=7
      integer :: pz, pa, iz, io_len, alloc_err
      integer(kind=8) :: rec_len
      character(len=4) :: filename
      logical :: luse_lnrho
      namelist / pgen  / varfile, irho, iux, iuy, iuz, ibx, iby, ibz, luse_lnrho
!
! default values
!
      mxgrid = in+1
      mygrid = jn+1
      mzgrid = kn+1
      nv = 7
!
! set default values
!
      varfile = 0
      iux = 1
      iuy = 2
      iuz = 3
      irho= 4
      ibx = 5
      iby = 6
      ibz = 7
      luse_lnrho = .true.

      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
        ibuf_in(1) = varfile
        ibuf_in(2) = irho
        ibuf_in(3) = iux
        ibuf_in(4) = iuy
        ibuf_in(5) = iuz
        ibuf_in(6) = ibx
        ibuf_in(7) = iby
        ibuf_in(8) = ibz    
      endif
      call MPI_BCAST( ibuf_in, 8, MPI_INTEGER &
                     , 0, comm3d, ierr )
!      call MPI_BCAST(ibuf_in, 13, MPI_FLOAT &
!                     , 0, comm3d, ierr )
!
! write filename and open file, 
! then read in VAR file
!   
      if (varfile .eq. -1) then
        filename = 'var.dat'
      else
        write(filename, "(A3,I1)") 'VAR',varfile
      endif
!
! allocate arrays 
!
      !if (myid .eq. 0) then
        allocate(ga(mxgrid,mygrid,mzgrid,nv), stat=alloc_err)
        allocate(buffer(mxgrid,mygrid), stat=alloc_err)
!
        if (ip <= 8) print *, 'input_snap: open VAR0'
        inquire (IOLENGTH=io_len) tvar
        open (lun_input, FILE=filename, access='direct', &
              recl=mxgrid*mygrid*io_len, status='old')
!
        if (ip <= 8) print *, 'input_snap: read dim=', mxgrid, mygrid, mzgrid, nv
        ! iterate through variables
        do pa = 1, nv
          ! iterate through xy-planes and read each plane separately
          do iz = 1, mzgrid
            read (lun_input, rec=iz+(pa-1)*mzgrid) buffer
            ga(:,:,iz,pa) = buffer
          enddo
        enddo
!
! allocate x,y & z arrays and read the data
!
        allocate(x(mxgrid),y(mygrid),z(mzgrid), stat=alloc_err)
!
        close (lun_input)
        open (lun_input, FILE=filename, & 
              FORM='unformatted', status='old', position='append')
        call backskip_to_time(lun_input)
!
        read (lun_input) tvar, x, y, z, dx, dy, dz
        if (ip <= 8) print *, tvar, dx, x(1)
    
        if (ip <= 8) print *, shape(v1) 
        if (luse_lnrho) then
          d = exp(ga(1:in, 1:jn, 1:kn, irho))
        else
          d = ga(1:in, 1:jn, 1:kn, irho)
        endif
        v1 = ga(1:in, 1:jn, 1:kn, iux)
        v2 = ga(1:in, 1:jn, 1:kn, iuy)
        v3 = ga(1:in, 1:jn, 1:kn, iuz)
        b1 = ga(1:in, 1:jn, 1:kn, ibx)
        b2 = ga(1:in, 1:jn, 1:kn, iby)
        b3 = ga(1:in, 1:jn, 1:kn, ibz)
      !endif

      write(*,*) maxval(v1), maxval(b1)
      write(*,*) maxval(v2), maxval(b2)
      write(*,*) maxval(v3), maxval(b3)
      write(*,*) minval(d), maxval(d)
      if (allocated(ga)) deallocate(ga)
      return
end subroutine turb_simple
