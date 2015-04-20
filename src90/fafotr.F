c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               F A F O T R                 \\\\\\\\\\
c
c=======================================================================
c
c relabeling of Fourier coefficients, preparing for FFT, call FFT, write
c result into A-matrix
C
C    modified 1: write to fft data stream improved   
C                03/01 Robi Banerjee
C
c    modified 2: use of math-pack fft -> opitmized for vector computers
c                12/2001 Robi Banerjee
c
c    modified 3: reduce of 3D variables
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fafotr( ar, ai )
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use fftpars
      implicit none
      include "fftw_f77.i"

      real(rl), dimension(in,jn,kn), intent(inout)  :: ar, ai
      integer  :: ix, iy, iz, jx, jy, jz, plan
c      integer isign
      real(rl) :: re_av, re_rms
      real(rl) :: im_av, im_rms
  
      allocate(fftdata(in, jn, kn))
      fftdata = (0.0d0, 0.0d0)

c -------------------------------------------------------------------------
c -------- vector field 1 direction ---------------------------------------
c -------------------------------------------------------------------------
c
c     relable & write into string

      do 30 iz=1, kn
        if (iz.le.(kn/2+1)) then
           jz = iz + kn/2 - 1 
        else
           jz = iz - kn/2 - 1
        endif         
        do 20 iy=1, jn       
          if (iy.le.(jn/2+1)) then
              jy = iy + jn/2 - 1 
          else
              jy = iy - jn/2 - 1
          endif
 
          do 10 ix=1, in
            if (ix.le.(in/2+1)) then
                 jx = ix + in/2 - 1 
            else
                 jx = ix - in/2 - 1
            endif                      
            fftdata(ix,iy,iz) = cmplx(  ar(jx,jy,jz)
     &                               , ai(jx,jy,jz) )
10        continue ! ix
20      continue ! iy
30    continue ! iz
c      write(*,*) 'dimensions of Arrays:', in, jn, kn, N
c      write(*,*) 'Data sanity check before FFT'
c      write(*,*) fftdata(1,1,1), fftdata(10,10,10)
c      write(*,*) fftdata(N/2,N/2, N/2),  fftdata(N,N,N)
c
c     Create and execute FFTw (2.1.5) plan
c
      call fftw3d_f77_create_plan(plan,in,jn,kn, 
     &                  FFTW_FORWARD, FFTW_ESTIMATE + FFTW_IN_PLACE)
      call fftwnd_f77_one(plan, fftdata, 0)
      call fftwnd_f77_destroy_plan(plan)       
c      write(*,*) 'Data sanity check after FFT'
c      write(*,*) size(fftdata,1), size(fftdata,2),size(fftdata,3),N,ijkn
c      write(*,*) fftdata(1,1,1), fftdata(10,10,10)
c      write(*,*) fftdata(N/2,N/2, N/2),  fftdata(N,N,N)

c     write result into result-matrix
c
      re_av  = 0.0
      re_rms = 0.0
      im_av  = 0.0
      im_rms = 0.0

      do 60 iz=1, kn
        do 50 iy=1, jn
          do 40 ix=1, in       

            ar(ix,iy,iz) = real(fftdata(ix,iy,iz))
            ai(ix,iy,iz) = imag(fftdata(ix,iy,iz))
         
            re_av  = re_av  + ar(ix,iy,iz)
            re_rms = re_rms + ar(ix,iy,iz)**2
            im_av  = im_av  + ai(ix,iy,iz)
            im_rms = im_rms + ai(ix,iy,iz)**2

40        continue
50      continue
60    continue

      write(6, 2020) 'im_av  = ', im_av/re_av
      write(6, 2020) 'im_rms = ', sqrt(im_rms/re_rms)


2010  format('FAFOTR   : ' )
2020  format('FAFOTR   : ', a10, 1pe13.6 )
2030  format('FAFOTR   : ', a10, 1pi7 )

      deallocate(fftdata)
      return
      end
c=======================================================================
c
c    \\\\\\\\\\      E N D   S U B R O U T I N E         //////////
c    //////////             F A F O T R                  \\\\\\\\\\
c
c=======================================================================
