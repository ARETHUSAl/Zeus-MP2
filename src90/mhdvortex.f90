!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////            M H D S H K T U B E            \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!     PURPOSE: initializes Orszag-Tang Vortex test problem.
!
!     Written by: John Hayes
!
!=======================================================================
subroutine mhdvortex
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
      integer  :: i, j, k
!
      real(rl) :: b0, Az(in,jn,kn), aism3(in,jn,kn), ajsm3(in,jn,kn), &
                                    aksm3(in,jn,kn)
!
      real(rl) :: foo ! dummy argument to keep PGEN as a placeholder
!
      namelist / pgen     /  foo
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
       endif
!       call MPI_BCAST( buf_in, 12, MPI_DOUBLE_PRECISION
!     &               , 0, comm3d, ierr )
       if (myid .ne. 0) then
       endif
!
!----------------------------------------------------------------------- 
!     Initialize hydro variables
!----------------------------------------------------------------------- 
!
      do k=1,kn
       do j=1,jn
        do i=1,in
         d (i,j,k) = 25.0D0 / (36.0D0 * pi)
         e (i,j,k) =  5.0D0 / (gamm1*(12.0D0 * pi))
         v1(i,j,k) = -sin(2.0D0*pi*x2b(j))
         v2(i,j,k) =  sin(2.0D0*pi*x1b(i))
         !v3(i,j,k) =  cos(2.0D0*pi*x1b(i))
         v3(i,j,k) = 0.0D0
        enddo
       enddo
      enddo
!
!----------------------------------------------------------------------- 
!     Initialize magnetic field
!----------------------------------------------------------------------- 
!
      b0 = 1.0D0 / sqrt(4.0D0*pi)
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         b1(i,j,k) = -b0*sin(2.0*pi*x2b(j))
         b2(i,j,k) =  b0*sin(4.0*pi*x1b(i))
         b3(i,j,k) = 0.0
        enddo
       enddo
      enddo
!
      return
end subroutine mhdvortex
