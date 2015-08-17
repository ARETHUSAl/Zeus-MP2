!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               S P C D U M P               \\\\\\\\\!
!=======================================================================
!
!
!    jms:zeus3d.spcdump<----------------------- contorols spectra  dumps
!                                                            Jan, 2002
!
!    written by: Robi Banerjee
!
!  PURPOSE: dumps spectra 
!
!  Currently, spectra dumps can be made for:
!
!   variable               meaning            file written to
!
!      b1     1-component of magnetic field       zub1NNNXX
!      b2     2-component of magnetic field       zub2NNNXX
!      b3     3-component of magnetic field       zub3NNNXX
!      d      density                             zud_NNNXX
!      v1     1-component of velocity             zuv1NNNXX
!      v2     2-component of velocity             zuv2NNNXX
!      v3     3-component of velocity             zuv3NNNXX
!      h      helicity power spectrum             zuh_NNNXX
!
!  where NNN is a three digit integer which distinguishes the spc files
!  dumped during the run, and XX is a two-character id tag specified in
!  "rescon".  
!
!  LOCAL VARIABLES:
!
!  EXTERNALS:
!    MK1DSPC
!
!-----------------------------------------------------------------------
! 
      subroutine spcdump
      use real_prec
      use param
      use field
      use root
      use scratch
      use fftpars
      use mpiyes
      use mpipar
      implicit NONE
!
      deallocate(fftdata)
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               S P C D U M P               \\\\\\\\\!
!=======================================================================
!
