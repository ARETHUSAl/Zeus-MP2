!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 S E D O V                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
       subroutine sedov
!
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k 
!
      real(rl) :: d0, t0, e0, eblast, rblast, eb0
!
      namelist / pgen     / d0, t0, eblast, rblast
!
!-----------------------------------------------------------------------
!
       d0     = 1.0d-8
       t0     = 50.0
       eblast = 1.0d50
       rblast = 1.0d12
!
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = d0  
         buf_in(2) = t0  
         buf_in(3) = eblast
         buf_in(4) = rblast
#endif
       endif
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 4, MPI_FLOAT &
                      , 0, comm3d, ierr )
        if (myid .ne. 0) then
         d0     = buf_in(1)
         t0     = buf_in(2)
         eblast = buf_in(3)
         rblast = buf_in(4)
        endif ! myid
#endif
!
!     Set up atmosphere.
!
      e0 = d0*boltz*t0/(gamm1*mmw*mh)
      do 30 k=1,kn
        do 20 j=1,jn
          do 10 i=1,in
            d (i,j,k) = d0
            v1(i,j,k) = 0.0D0
            v2(i,j,k) = 0.0D0
            v3(i,j,k) = 0.0D0
            e (i,j,k) = e0
10        continue
20      continue
30    continue
!
!     Set up central region.
!
      eb0 = 3.0d0*eblast/(4.0D0*pi*rblast**3)
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         if(x1b(i) .lt. rblast) e(i,j,k) = e(i,j,k) + eb0
        enddo
       enddo
      enddo
!
      return
      end
!
