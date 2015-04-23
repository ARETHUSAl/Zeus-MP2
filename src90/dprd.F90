!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               L S _ D P R D               \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!     PURPOSE:  compute a vector dot product
!
!     Written by: F. Douglas SwestyF. Douglas Swesty
!
!=======================================================================
      real*8 function ls_dprd(isx,iex,isy,iey,isz,iez,v1,v2)
!
      use real_prec
      use config
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez
!
      real(rl) :: v1(neqm,in,jn,kn), &
                  v2(neqm,in,jn,kn)
      real(rl) ::  psum, dotprd
!
      integer  :: ix, iy, iz, j, ierror
!
      psum = 0.0d0
      do iz = isz,iez,1
       do iy = isy,iey,1
        do ix = isx,iex,1
         do j = 1,neqm,1
           psum = psum+v1(j,ix,iy,iz)*v2(j,ix,iy,iz)
         enddo
        enddo
       enddo
      enddo
!
#ifdef MPI_USED
      call mpi_allreduce(psum,dotprd,1,MPI_FLOAT,mpi_sum, &
                         comm3d,ierror) 
#endif
#ifndef MPI_USED
      dotprd = psum
#endif
!
      ls_dprd = dotprd
!
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D    S U B R O U T I N E      //////////
!    //////////                L S _ D P R D             \\\\\\\\\\
!
!=======================================================================
