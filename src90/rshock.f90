!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                R S H O C K                \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine rshock
!
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k
!
      real(rl) :: d0, e0, temp, dtemp
!
      namelist /pgen/ d0, e0
!
!     initialize and read in parameters from PGEN namelist
!
      d0 = 1.D0
      e0 = 1.D-4
!
      if (myid_w .eq. 0) then
       read (1,pgen)
       write(2,pgen)
       buf_in(1) = d0
       buf_in(2) = e0
      endif
      call MPI_BCAST( buf_in, 2, MPI_DOUBLE_PRECISION &
                     , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0 = buf_in(1)
       e0 = buf_in(2)
      endif
!
      gamm1 = gamma - 1.0D0
!
!     set up density distribution
!
      do 1 k = 1, kn
      do 1 j = 1, jn
      do 1 i = 1, in
       d(i,j,k) = d0
       p(i,j,k) = -1.0
1     continue
!
!     initialize gas and radiation energy densities, velocities
!
      fois( 2) = e0*d0*boltz/((gamma-1.0)*mmw*mh)
      temp     = e0
      fois(12) = rad_con * temp**4
      dtemp = -75.D0/7.D10
      do 2 k=1,kn
      do 2 j=1,jn
      do 2 i=is,ie
       temp      = e0 + dtemp*(x1b(i)-7.0D10)
       e (i,j,k) = temp*d0*boltz/((gamma-1.0)*mmw*mh)
       if(lrad .eq. 1) then
        er(i,j,k) = rad_con * temp**4
       endif
       v1(i,j,k) = fois(3)
       v2(i,j,k) = 0.D0
       v3(i,j,k) = 0.D0
2     continue
!
!     set up outer 1-boundary values
!
      do k = 1, kn
      do j = 1, jn
       eoib(j,k,1) = fois( 2)
       eoib(j,k,2) = fois( 2)
       if(lrad .eq. 1) then
        eroib(j,k,1) = fois(12)
        eroib(j,k,2) = fois(12)
       endif
      enddo
      enddo
!
      return
      end
      subroutine rshockres
!
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k
!
      real(rl) :: d0, e0, temp, dtemp
!
      namelist /pgen/ d0, e0
!
!     initialize and read in parameters from PGEN namelist
!
      d0 = 1.D0
      e0 = 1.D-4
!
      if (myid_w .eq. 0) then
       read (1,pgen)
       write(2,pgen)
       buf_in(1) = d0
       buf_in(2) = e0
      endif
      call MPI_BCAST( buf_in, 2, MPI_DOUBLE_PRECISION &
                     , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0 = buf_in(1)
       e0 = buf_in(2)
      endif
!
      gamm1 = gamma - 1.0D0
!
      return
      end
