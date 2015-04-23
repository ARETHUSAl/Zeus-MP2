!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              R A D S O L V E              \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine rad_solve
!
      use real_prec
      use config
      use root
      use field
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
!-----------------------------------------------------------------------
!     radiation solver menu:
!    
!     LRAD = 0:  No radiation transport
!     LRAD = 1:  Grey, flux-limited (non-equilibrium) diffusion
!     LRAD > 1:  User-supplied options
!-----------------------------------------------------------------------
!
      if(lrad .eq. 1) then
       call grey_fld(w3dd, e, w3de, w3dh, er)
       call rad_imp_dt(w3dh,er)
       dtnri2 = dtnri2 / (1.26*1.26)
      endif
      if(lrad .gt. 1) then
       if(myid .eq. 0) &
          write(*,"('RAD_SOLVE: Unsupported value of LRAD')")
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
!
      return
      end
