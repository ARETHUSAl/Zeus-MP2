!=======================================================================
!
!    \\\\\\\\\\      B E G I nn   S U B R O U T I nn E      //////////
!    //////////               M K 1 D S P C               \\\\\\\\\!
!=======================================================================
!
      subroutine mk1dspc
      use real_prec
      use param
      use field
      use root
      use scratch
      use fftpars
      use mpiyes
      use mpipar
      implicit NONE
      include "fftw_f77.i"      
!    ???:zeus3d.???????? <----------------------------------------------
!                                                        December, 2001
!
!    written by: Robi Banerjee
!    modified 1: ?
!
!  PURPOSE:  makes 1D spectra of the required variable (b-field, 
!            v-field, density)
!            adds power spectra to pwpsc
!            
!
!  INPUT VARIABLES: 
!        nbins    : number of bins of the 1D spectra
!
!  OUTPUT VARIABLES:
!
!  LOCAL VARIABLES:
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
!
! INPUT VARIABLES
!
      integer nbins
!
!  LOCAL VARIABLES:
!
!      integer isign
!      REAL    scale
      integer   :: plan
      integer   :: i,j,k, idx, inn(3)
      integer   :: N 
      integer   :: k1,k2,k3
      real(rl)  :: kmax, kmag
!      real(rl), dimension(N+1),intent(inout) :: pwspc
      inn(1) = in
      inn(2) = jn
      inn(3) = kn
      kmax   = sqrt(3.e0)*real(N/2)
      N      = ijkn
      nbins  = N
!      isign = -1
!      scale = 1.0/(real(N))**3
!
!      Create and execute FFTw (2.1.5) plan
!
       call fftwnd_f77_create_plan(plan,3,inn, &
                        FFTW_FORWARD, FFTW_ESTIMATE + FFTW_IN_PLACE)
       call fftwnd_f77_one(plan, fftdata, 0)
       call fftwnd_f77_destroy_plan(plan)       
!      Routines for FFTW-3.3.4
!
!      call dfftw_plan_dft_3d(plan, N, N, N, data, data,
!     &                       FFTW_FORWARD, FFTW_ESTIMATE)
!      call dfftw_execute_dft(plan, data, data )
!      call dfftw_destroy_plan(plan)
!      call zzfft3d(isign, N, N, N, scale, data, in, jn
!    &             ,data, in, jn, table, work, 0)
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
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               M K 1 D S P C               \\\\\\\\\!
!=======================================================================
