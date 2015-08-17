!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                   N O H                   \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
!
       subroutine noh
!
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i  , j, k, ip1, jp1, &
                  kp1, idirect
!
      real(rl) :: &
                  d0 , p0 , e0, v0
!
      namelist / pgen     / &
                     p0, d0, v0, idirect
!
!-----------------------------------------------------------------------
!
       d0 =  1.0D0
       p0 =  1.0D-6
       v0 = -1.0D0
!
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
         buf_in(1)   = d0
         buf_in(2)   = p0 
         buf_in(3)   = v0  
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 3, MPI_DOUBLE_PRECISION &
                     , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 1, MPI_INTEGER &
                     , 0, comm3d, ierr )
       if (myid .ne. 0) then
         d0      = buf_in(1)
         p0      = buf_in(2)
         v0      = buf_in(3)
         idirect = ibuf_in( 1)
       endif
!
!      Set up initial state.
!
      do k=ks,ke
       do j=js,je
        do i=is,ie
         d(i,j,k) = d0
         e(i,j,k) = p0/gamm1
        enddo
       enddo
      enddo
!
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         if(idirect .eq. 1) then
          v1(i,j,k) = v0
          v2(i,j,k) = 0.0
          v3(i,j,k) = 0.0
         else if(idirect .eq. 2) then
          v2(i,j,k) = v0
          v1(i,j,k) = 0.0
          v3(i,j,k) = 0.0
         else
          v3(i,j,k) = v0
          v1(i,j,k) = 0.0
          v2(i,j,k) = 0.0
         endif
        enddo
       enddo
      enddo
!
      if(idirect .eq. 1) then
       do k = 1, kn
        do j = 1, jn
         v1oib(j,k,1) = v0
         v1oib(j,k,2) = v0
         v2oib(j,k,1) = 0.0
         v2oib(j,k,2) = 0.0
         v3oib(j,k,1) = 0.0
         v3oib(j,k,2) = 0.0
         doib(j,k,1) = d0
         doib(j,k,2) = d0
         eoib(j,k,1) = p0/gamm1
         eoib(j,k,2) = p0/gamm1
        enddo
       enddo
      else if(idirect .eq. 2) then
       do k = 1, kn
        do i = 1, in
         v2ojb(i,k,1) = v0
         v2ojb(i,k,2) = v0
         v1ojb(i,k,1) = 0.0
         v1ojb(i,k,2) = 0.0
         v3ojb(i,k,1) = 0.0
         v3ojb(i,k,2) = 0.0
         dojb(i,k,1) = d0
         dojb(i,k,2) = d0
         eojb(i,k,1) = p0/gamm1
         eojb(i,k,2) = p0/gamm1
        enddo
       enddo
      else
       do j = 1, jn
        do i = 1, jn
         v3okb(i,j,1) = v0
         v3okb(i,j,2) = v0
         v1okb(i,j,1) = 0.0
         v1okb(i,j,2) = 0.0
         v2okb(i,j,1) = 0.0
         v2okb(i,j,2) = 0.0
         dokb(i,j,1) = d0
         dokb(i,j,2) = d0
         eokb(i,j,1) = p0/gamm1
         eokb(i,j,2) = p0/gamm1
        enddo
       enddo
      endif
!
      return
      end
