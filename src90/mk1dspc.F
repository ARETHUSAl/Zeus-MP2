c=======================================================================
c
c    \\\\\\\\\\      B E G I nn   S U B R O U T I nn E      //////////
c    //////////               M K 1 D S P C               \\\\\\\\\\
c
c=======================================================================
c
      subroutine mk1dspc
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
      include "fftw_f77.i"      

c    ???:zeus3d.???????? <----------------------------------------------
c                                                        December, 2001
c
c    written by: Robi Banerjee
c    modified 1: ?
c
c  PURPOSE:  makes 1D spectra of the required variable (b-field, 
c            v-field, density)
c            adds power spectra to pwpsc
c            
c
c  INPUT VARIABLES: 
c        nbins    : number of bins of the 1D spectra
c
c  OUTPUT VARIABLES:
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
c INPUT VARIABLES
c
      integer nbins
c
c  LOCAL VARIABLES:
c
c      integer isign
c      REAL    scale
      integer   :: plan
      integer   :: i,j,k, idx, inn(3)
      integer   :: N 
      integer   :: k1,k2,k3
      real(rl)  :: kmax, kmag
c      real(rl), dimension(N+1),intent(inout) :: pwspc
      inn(1) = in
      inn(2) = jn
      inn(3) = kn
      kmax   = sqrt(3.e0)*real(N/2)
      N      = ijkn
      nbins  = N

c      isign = -1
c      scale = 1.0/(real(N))**3
c
c      Create and execute FFTw (2.1.5) plan
c
       call fftwnd_f77_create_plan(plan,3,inn, 
     &                  FFTW_FORWARD, FFTW_ESTIMATE + FFTW_IN_PLACE)
       call fftwnd_f77_one(plan, fftdata, 0)
       call fftwnd_f77_destroy_plan(plan)       

c      Routines for FFTW-3.3.4
c
c      call dfftw_plan_dft_3d(plan, N, N, N, data, data,
c     &                       FFTW_FORWARD, FFTW_ESTIMATE)
c      call dfftw_execute_dft(plan, data, data )
c      call dfftw_destroy_plan(plan)
c      call zzfft3d(isign, N, N, N, scale, data, in, jn
c    &             ,data, in, jn, table, work, 0)

      do 130 k=1,N
        if (k.gt.N/2) then 
          k3 = k - N - 1
        else
          k3 = k - 1
        endif
        do 120 j=1,N
          if (j.gt.N/2) then 
            k2 = j - N - 1
          else
            k2 = j - 1
          endif
          do 110 i=1,N
            if (i.gt.N/2) then 
              k1 = i - N
            else
              k1 = i - 1
            endif
            kmag = sqrt(real(k1**2) + real(k2**2) + real(k3**2))
            idx  = int(kmag/kmax*real(nbins)) + 1
            pwspc(idx) = pwspc(idx)+abs(fftdata(i,j,k))**2
110       continue
120     continue
130   continue

      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               M K 1 D S P C               \\\\\\\\\\
c
c=======================================================================
