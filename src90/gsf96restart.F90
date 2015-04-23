!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          G S F 9 6 R E S T A R T          \\\\\\\\\\
!
!=======================================================================
!
      subroutine gsf96restart
!
!
!  written by: Daniel Whalen 11-30-05
!
!  ported to ZEUS-MP 2.1 by DJW 11.29.06
!
!  PURPOSE:    Restarts the gsf96 problem by rereading in pgen variables
!              
!            
!            
!  LOCAL VARIABLES:
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
      real(rl) :: t_backg, t_cool, r_trans, n_central, r_core
      namelist / pgen     / usrtag, omega, n_central, r_core, &
                            t_backg, r_trans
      usrtag    =  'usr'
      omega     = -2.0
      n_central = 1.0e05
      r_core    = 2.1e16
      t_backg   = 100.0
      r_trans   = 3.084e17
      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)
#ifdef MPI_USED
        buf_in(1) = omega 
        buf_in(2) = n_central
        buf_in(3) = r_core
        buf_in(4) = t_backg
        buf_in(5) = r_trans
      endif
      call MPI_BCAST( buf_in, 5, MPI_FLOAT &
                     , 0, comm3d, ierr )
      call MPI_BCAST( usrtag, 3, MPI_CHARACTER &
                     , 0, comm3d, ierr )
      if (myid .ne. 0) then
        omega     = buf_in(1)
        n_central = buf_in(2)
        r_core    = buf_in(3)
        t_backg   = buf_in(4)
        r_trans   = buf_in(5)
#endif /* MPI_USED */
      endif
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////          G S F 9 6 R E S T A R T          \\\\\\\\\\
!
!=======================================================================
!
