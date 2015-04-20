c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               S P C D U M P               \\\\\\\\\\
c
c=======================================================================
c
c
c    jms:zeus3d.spcdump<----------------------- contorols spectra  dumps
c                                                            Jan, 2002
c
c    written by: Robi Banerjee
c
c  PURPOSE: dumps spectra 
c
c  Currently, spectra dumps can be made for:
c
c   variable               meaning            file written to
c
c      b1     1-component of magnetic field       zub1NNNXX
c      b2     2-component of magnetic field       zub2NNNXX
c      b3     3-component of magnetic field       zub3NNNXX
c      d      density                             zud_NNNXX
c      v1     1-component of velocity             zuv1NNNXX
c      v2     2-component of velocity             zuv2NNNXX
c      v3     3-component of velocity             zuv3NNNXX
c      h      helicity power spectrum             zuh_NNNXX
c
c  where NNN is a three digit integer which distinguishes the spc files
c  dumped during the run, and XX is a two-character id tag specified in
c  "rescon".  
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c    MK1DSPC
c
c-----------------------------------------------------------------------
c 
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
c
      character(len=9) :: filename
      integer          :: i, j, k, ip, jp, kp
      integer          :: N 
      integer          :: nbins
      real(rl)         :: dbin
      allocate(fftdata(in, jn, kn))
c
c
c      External statements
c
c      external      mk1dspc
c
c-----------------------------------------------------------------------
c
       N      = ijkn
      nbins   = N
      dbin    = sqrt(3.e0)*real(N/2)/real(nbins)

c
c------ initialize fft table -------------------------------------------
c
c      call zzfft3d(0, N, N, N, 0.0, fftdata, in, jn
c     &             ,fftdata, in, jn, table, work, 0)
c-----------------------------------------------------------------------
c------ Generate filename for magenetic power spectra ------------------
c-----------------------------------------------------------------------
      write (filename, 2010) 'b ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'

      pwspc = 0.e0

c
c-----magnetic field 1 component --------------------------------------
c
      do 30 k=1,N
        do 20 j=1,N
          do 10 i=1,N
            fftdata(i,j,k) = cmplx( b1(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
10        continue
20      continue
30    continue

      call mk1dspc
c
c-----magnetic field 2 component --------------------------------------
c
      do 130 k=1,N
        do 120 j=1,N
          do 110 i=1,N
            fftdata(i,j,k) = cmplx( b2(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
110       continue
120     continue
130   continue

      call mk1dspc
c
c-----magnetic field 3 component --------------------------------------
c
      do 230 k=1,N
        do 220 j=1,N
          do 210 i=1,N
            fftdata(i,j,k) = cmplx( b3(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
210       continue
220     continue
230   continue

      call mk1dspc
c
c-----write magnetic power spectra to file -----------------------------
c
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 310 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
310   continue
      close ( unit = iousr )

      if (iolog .gt. 0)
     1  write (iolog, 2020) filename,time, nhy

c----------------------------------------------------------------------
c----- Generate filename for velocity power spectra -------------------
c----------------------------------------------------------------------
      write (filename, 2010) 'v ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'

      pwspc = 0.e0

c
c-----velocity 1 component --------------------------------------
c
      do 430 k=1,N
        do 420 j=1,N
          do 410 i=1,N
            fftdata(i,j,k) = cmplx( v1(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
410       continue
420     continue
430   continue

      call mk1dspc
c
c-----velocity 2 component --------------------------------------
c
      do 530 k=1,N
        do 520 j=1,N
          do 510 i=1,N
            fftdata(i,j,k) = cmplx( v2(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
510       continue
520     continue
530   continue

      call mk1dspc
c
c-----velocity 3 component --------------------------------------
c
      do 630 k=1,N
        do 620 j=1,N
          do 610 i=1,N
            fftdata(i,j,k) = cmplx( v3(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
610       continue
620     continue
630   continue

      call mk1dspc
c
c-----write velocity power spectra to file ----------------------
c
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 710 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
710   continue
      close ( unit = iousr )

      if (iolog .gt. 0)
     1  write (iolog, 2020) filename, time, nhy
c
c-----------------------------------------------------------------------
c------ Generate filename for helicity power spectra -------------------
c-----------------------------------------------------------------------
      write (filename, 2010) 'h ', iusrdmp, id
      if (filename(4:4) .eq. ' ') filename(4:4) = '_'

      pwspc = 0.e0

c
c-----helicity --------------------------------------------------------
c
      do 810 k=ksmn,kemx
        kp = k+1
        do 820 j=jsmn,jemx
          jp = j+1
          do 830 i=ismn,iemx
            ip = i+1
            wa3d(i,j,k) = 
     &               ( ( a1(i ,j ,k ) + a1(i ,j ,kp) +
     &                   a1(i ,jp,k ) + a1(i ,jp,kp) ) * 
     &                 ( b1(i ,j ,k ) + b1(ip,j ,k ) ) +
     &                 ( a2(i ,j ,k ) + a2(i ,j ,kp) +
     &                   a2(ip,j ,k ) + a2(ip,j ,kp) ) *
     &                 ( b2(i ,j ,k ) + b2(i ,jp,k ) ) +
     &                 ( a3(i ,j ,k ) + a3(i ,jp,k ) +
     &                   a3(ip,j ,k ) + a3(ip,jp,k ) ) * 
     &                 ( b3(i ,j ,k ) + b3(i ,j ,kp) ) )
830       continue
820     continue
810   continue

      do 930 k=1,N
        do 920 j=1,N
          do 910 i=1,N
            fftdata(i,j,k) = cmplx( wa3d(i+ismn-1,j+jsmn-1,k+ksmn-1)
     &                         , 0.0e0 )
910       continue
920     continue
930   continue

      call mk1dspc
c
c-----write helicity power spectra to file ----------------------
c
      open ( unit=iousr, file=filename, status='unknown' )
      write(iousr, 2040) time, nhy
      do 1010 i=1,N
        write(iousr, 2030) real(i-1)*dbin, pwspc(i)
1010  continue
      close ( unit = iousr )

      if (iolog .gt. 0)
     1  write (iolog, 2020) filename, time, nhy


c-----------------------------------------------------------------------
c----- Increment USER file counter -------------------------------------
c-----------------------------------------------------------------------
      iusrdmp = iusrdmp + 1
c
c-----------------------------------------------------------------------
c----------------------- Write format statements -----------------------
c-----------------------------------------------------------------------
c
2010  format('zu',a2,i3.3,a2)
2020  format('SPCDUMP : SPC dump ',a9,' stored at time ='
     1      ,1pe12.5,', nhy =',i6)
2030  format(1pe12.5, 1pe12.5)
2040  format('# spectra dump at time =', 'nhy =',i6)
#endif /* SPECTRA */
c
      deallocate(fftdata)
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               S P C D U M P               \\\\\\\\\\
c
c=======================================================================
c
