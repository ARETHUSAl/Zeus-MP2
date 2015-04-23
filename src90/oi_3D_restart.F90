!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////         O I _ 3 D _ R E S T A R T         \\\\\\\\\\
!
!=======================================================================
!
      subroutine oi_3D_restart
!
!
!  written by: Daniel Whalen 06-16-06
!
!  ported to ZEUS-MP 2.1 by DJW 11.29.06
!
!  PURPOSE:  restarts the minihalo external photoevaporation problem.
!            
!
!  EXTERNALS:  none
!
!-----------------------------------------------------------------------
!
      use real_prec
      use param
      use cons
      use config
      use root
      use chem
      use field
      use bndry
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
      integer  :: nhalo, ihalo 
      real(rl) :: ovrdns, yc, zc, m_halo, mu, f_H2, X_e
      namelist / pgen     / nhalo, ovrdns, xc, yc, zc, &
                            mu, ihalo, m_halo, f_H2, X_e
      nhalo  = 11
      ihalo  = 2
      ovrdns = 1.0
      yc     = 0.0
      zc     = 0.0
      mu     = 1.22
      m_halo = 2.0e05
      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)
#ifdef MPI_USED
        buf_in(1)  = ovrdns 
        buf_in(2)  = xc
        buf_in(3)  = yc
        buf_in(4)  = zc
        buf_in(5)  = mu
        buf_in(6)  = m_halo
        buf_in(7)  = X_e
        buf_in(8)  = f_H2
        ibuf_in(1) = ihalo
        ibuf_in(2) = nhalo
      endif
      call MPI_BCAST( buf_in ,8, MPI_FLOAT &
                    , 0, comm3d, ierr )
      call MPI_BCAST( ibuf_in,2, MPI_INTEGER &
                    , 0, comm3d, ierr )
      if (myid .ne. 0) then
        ovrdns    = buf_in(1)
        xc        = buf_in(2)
        yc        = buf_in(3)
        zc        = buf_in(4)
        mu        = buf_in(5)
        m_halo    = buf_in(6)
        X_e       = buf_in(7)
        f_H2      = buf_in(8)
        ihalo     = ibuf_in(1)
        nhalo     = ibuf_in(2)
#endif /* MPI_USED */
      endif
!#ifdef MPI_USED
!      call MPI_BARRIER(comm3d, ierr)
!#endif /* MPI_USED */
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////         O I _ 3 D _ R E S T A R T         \\\\\\\\\\
!
!=======================================================================
!
