c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           L C D M R E S T A R T           \\\\\\\\\\
c
c=======================================================================
c
      subroutine lcdmrestart
c
c
c  written by: Daniel Whalen 08-26-04
c
c  ported to ZEUS-MP 2.1 by DJW 11.29.06
c
c  PURPOSE:  Restart module for the lcdm problem (reads in pgen vbls)  
c            
c            
c 
c
c  EXTERNALS:  none
c
c-----------------------------------------------------------------------
c
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
c
      implicit NONE

      integer  ::  nhalo

      real(rl) ::  ovrdns, t_halo, r_trans

      namelist / pgen     / usrtag, ovrdns, r_trans, nhalo

      usrtag  = 'usr'
      ovrdns  = 1.0
      nhalo   = 1
      r_trans = 0.

      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)
#ifdef MPI_USED
        buf_in (1) = ovrdns 
        buf_in (2) = r_trans
        ibuf_in(1) = nhalo
      endif
      call MPI_BCAST( buf_in, 2, MPI_FLOAT
     .              , 0, comm3d, ierr )
      call MPI_BCAST( ibuf_in,1, MPI_INTEGER
     .              , 0, comm3d, ierr )
      call MPI_BCAST( usrtag, 3, MPI_CHARACTER
     .              , 0, comm3d, ierr )
      if (myid .ne. 0) then
        ovrdns    = buf_in(1)
        r_trans   = buf_in(2)
        nhalo     = ibuf_in(1)
#endif /* MPI_USED */
      endif

      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////           L C D M R E S T A R T           \\\\\\\\\\
c
c=======================================================================
c


