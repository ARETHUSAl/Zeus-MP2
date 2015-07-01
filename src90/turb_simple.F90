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
      implicit none
      integer :: i,j,k
      integer :: varfile, irho, iux, iuy, iuz, ibx, iby, ibz
      integer :: mx, my, mz, mpar
      integer :: st
      real(rl), dimension(:,:,:,:), allocatable :: f
      namelist / pgen  / varfile, irho, iux, iuy, iuz, ibx, iby, ibz
!
!      default values
!      fix all units according to the values of rho, c_s at z = 0 
!
      mx = in+6
      my = jn+6
      mz = kn+6
!
! set default values
!
      varfile = 0
      irho = 1
      iux = 2
      iuy = 3
      iuz = 4
      ibx = 5
      iby = 6
      ibz = 7

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
! allocate memory for VAR file
!
      allocate(f(mx,my,mz,mpar), stat=st)
!
! write filename and open file, 
! then read in VAR file
!   
      if (varfile .eq. 0) then
        filename = 'VAR0'
      else if (varfile .eq. -1) then
        filename = 'var.dat'
      else
        write(filename, "(A3,I1)") varfile
      endif
      open(unit=22, file=varfile, form='unformatted', status='old')
      read(unit=22) f
     
      d  = f(2:in+1, 2:jn+1, 2:kn+1, irho)
      v1 = f(2:in+1, 2:jn+1, 2:kn+1, iux)
      v2 = f(2:in+1, 2:jn+1, 2:kn+1, iuy)
      v3 = f(2:in+1, 2:jn+1, 2:kn+1, iuz)
      b1 = f(2:in+1, 2:jn+1, 2:kn+1, ibx)
      b2 = f(2:in+1, 2:jn+1, 2:kn+1, iby)
      b3 = f(2:in+1, 2:jn+1, 2:kn+1, ibz)

      write(*,*) maxval(v1), maxval(b1)
      write(*,*) maxval(v2), maxval(b2)
      write(*,*) maxval(v3), maxval(b3)
      if (allocated(f)) deallocate(f)
      return
end subroutine turb_simple
