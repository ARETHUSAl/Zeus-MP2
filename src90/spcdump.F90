!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               S P C D U M P               \\\\\\\\\\
!
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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      implicit NONE
#ifdef SPECTRA
!
      character(len=9) :: filename
      integer          :: i, j, k, ip, jp, kp
      integer          :: N 
      integer          :: nbins
      real(rl)         :: dbin
      allocate(fftdata(in, jn, kn))
!
!
!      External statements
!
!      external      mk1dspc
!
!-----------------------------------------------------------------------
!
       N      = ijkn
      nbins   = N
      dbin    = sqrt(3.e0)*real(N/2)/real(nbins)
!
!------ initialize fft table -------------------------------------------
!
!      call zzfft3d(0, N, N, N, 0.0, fftdata, in, jn
!     &             ,fftdata, in, jn, table, work, 0)
!-----------------------------------------------------------------------
!------ Generate filename for magenetic power spectra ------------------
!-----------------------------------------------------------------------
      write (filename, 2010) 'b ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'
      pwspc = 0.e0
!
!-----magnetic field 1 component --------------------------------------
!
      do 30 k=1,N
        do 20 j=1,N
          do 10 i=1,N
            fftdata(i,j,k) = cmplx( b1(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
10        continue
20      continue
30    continue
      call mk1dspc
!
!-----magnetic field 2 component --------------------------------------
!
      do 130 k=1,N
        do 120 j=1,N
          do 110 i=1,N
            fftdata(i,j,k) = cmplx( b2(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
110       continue
120     continue
130   continue
      call mk1dspc
!
!-----magnetic field 3 component --------------------------------------
!
      do 230 k=1,N
        do 220 j=1,N
          do 210 i=1,N
            fftdata(i,j,k) = cmplx( b3(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
210       continue
220     continue
230   continue
      call mk1dspc
!
!-----write magnetic power spectra to file -----------------------------
!
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 310 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
310   continue
      close ( unit = iousr )
      if (iolog .gt. 0) &
        write (iolog, 2020) filename,time, nhy
!----------------------------------------------------------------------
!----- Generate filename for velocity power spectra -------------------
!----------------------------------------------------------------------
      write (filename, 2010) 'v ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'
      pwspc = 0.e0
!
!-----velocity 1 component --------------------------------------
!
      do 430 k=1,N
        do 420 j=1,N
          do 410 i=1,N
            fftdata(i,j,k) = cmplx( v1(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
410       continue
420     continue
430   continue
      call mk1dspc
!
!-----velocity 2 component --------------------------------------
!
      do 530 k=1,N
        do 520 j=1,N
          do 510 i=1,N
            fftdata(i,j,k) = cmplx( v2(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
510       continue
520     continue
530   continue
      call mk1dspc
!
!-----velocity 3 component --------------------------------------
!
      do 630 k=1,N
        do 620 j=1,N
          do 610 i=1,N
            fftdata(i,j,k) = cmplx( v3(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
610       continue
620     continue
630   continue
      call mk1dspc
!
!-----write velocity power spectra to file ----------------------
!
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 710 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
710   continue
      close ( unit = iousr )
      if (iolog .gt. 0) &
        write (iolog, 2020) filename, time, nhy
!
!-----------------------------------------------------------------------
!------ Generate filename for helicity power spectra -------------------
!-----------------------------------------------------------------------
      write (filename, 2010) 'h ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'
      pwspc = 0.e0
!
!-----helicity --------------------------------------------------------
!
      do 810 k=ksmn,kemx
        kp = k+1
        do 820 j=jsmn,jemx
          jp = j+1
          do 830 i=ismn,iemx
            ip = i+1
            wa3d(i,j,k) = &
                     ( ( a1(i ,j ,k ) + a1(i ,j ,kp) + &
                         a1(i ,jp,k ) + a1(i ,jp,kp) ) * &
                       ( b1(i ,j ,k ) + b1(ip,j ,k ) ) + &
                       ( a2(i ,j ,k ) + a2(i ,j ,kp) + &
                         a2(ip,j ,k ) + a2(ip,j ,kp) ) * &
                       ( b2(i ,j ,k ) + b2(i ,jp,k ) ) + &
                       ( a3(i ,j ,k ) + a3(i ,jp,k ) + &
                         a3(ip,j ,k ) + a3(ip,jp,k ) ) * &
                       ( b3(i ,j ,k ) + b3(i ,j ,kp) ) )
830       continue
820     continue
810   continue
      do 930 k=1,N
        do 920 j=1,N
          do 910 i=1,N
            fftdata(i,j,k) = cmplx( wa3d(i+ismn-1,j+jsmn-1,k+ksmn-1) &
                               , 0.0e0 )
910       continue
920     continue
930   continue
      call mk1dspc
!
!-----write helicity power spectra to file ----------------------
!
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 1010 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
1010  continue
      close ( unit = iousr )
      if (iolog .gt. 0) &
        write (iolog, 2020) filename, time, nhy
!-----------------------------------------------------------------------
!----- Increment USER file counter -------------------------------------
!-----------------------------------------------------------------------
      iusrdmp = iusrdmp + 1
!
!-----------------------------------------------------------------------
!----------------------- Write format statements -----------------------
!-----------------------------------------------------------------------
!
2010  format('zu',a2,i3.3,a2)
2020  format('SPCDUMP : SPC dump ',a9,' stored at time =' &
            ,1pe12.5,', nhy =',i6)
2030  format(1pe12.5, 1pe12.5)
2040  format('# spectra dump at time =', 'nhy =',i6)
#endif /* SPECTRA */
!
      deallocate(fftdata)
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               S P C D U M P               \\\\\\\\\\
!
!=======================================================================
!
