#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           M N U _ D I R _ F L X           \\\\\\\\\\
c
c=======================================================================
c
       subroutine mnu_dir_flx(kslice) 

c  written by: Daniel Whalen 9.06
c
c  ported to ZEUS-MP 2.1 by DJW 11.29.06
c
c  PURPOSE:  Perform multifrequency static radiative transfer to 
c            compute radiative primordial chemistry rate coefficients. 
c            They are evaluated every chemical timestep.  The processes
c            are:
c 
c            k24:   HI    + gamma   -> HII   + e
c            k25:   HeII  + gamma   -> HeIII + e
c            k26:   HeI   + gamma   -> HeII  + e
c            k27:   HM    + gamma   -> HI    + e
c            k28:   H2II  + gamma   -> HI    + HII
c            k29:   H2I   + gamma   -> H2II  + e
c            k30:   H2II  + gamma   -> 2HII  + e
c            k31:   H2I   + gamma   -> 2HI
c 
c            Only k27 and k28 are calculated and summed over the first
c            nnu1 energy bins that lie below 13.6 eV.  All are evaluated
c            for the nnu2 energy bins above the Lyman edge.
c
c            The code solves the static flux equation in photon
c            conserving form for a point source centered in an RTP
c            coordinate mesh (lgeom = 3) or for plane waves along the 
c            1-direction in XYZ or ZRP coordinates (lgeom = 1 or 2).  
c            Plane waves can be geometrically attenuated by setting 
c            iPWA = 1.
c 
c
c  BOUNDARY VALUES USED:
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

      integer  :: i, j, kslice, n, index, niter

      real(rl) :: mhi, nioniz, chii, hnu, kt, dum, xi, q1, tau, 
     .            taumin, tau_HI, tau_HeI, tau_HeII, tau_H2I, 
     .            tau_H2II1, tau_H2II2, nrm_HI, nrm_HeI, nrm_HeII, 
     .            nrm_H2I, nrm_H2II1, nrm_H2II2, norm, k24p, k25p, 
     .            k26p, N_H2, N_HM, tau0, logtau0, ltau_l, hmnm,  
     .            x,Fshield, b5, tau_HM, nrm_HM, logHMnm,sep,
     .            tau_dust, nrm_dust

      real(rl), dimension(40) :: bin


#ifndef UNROLL_I
#define UNROLL_I
#endif
c
c-----------------------------------------------------------------------
c

      if (iRT .eq. 0) goto 2000

      do j = js, je
CDIR$ UNROLL UNROLL_I
      do i = is, ie

#ifdef H
          k24   (i,j) = tiny
          k24mv (i,j) = tiny
          piHI  (i,j) = tiny
#endif /* H */
#ifdef He
          k25   (i,j) = tiny
          k25mv (i,j) = tiny
          k26   (i,j) = tiny
          k26mv (i,j) = tiny
          piHeI (i,j) = tiny
          piHeII(i,j) = tiny
#endif /* He */
#ifdef H2
          k27   (i,j) = tiny
          k28   (i,j) = tiny
          k29   (i,j) = tiny
          k30   (i,j) = tiny
          k31   (i,j) = tiny
#endif /* H2 */

      enddo
      enddo
     

      mhi       = 1.0/mh
      taumin    = 1.0e-03
      kt        = t_star/tevk

      tau_HI    = 0.
      tau_HeI   = 0.
      tau_HeII  = 0.
      tau_HM    = 0.
      tau_H2I   = 0.
      tau_H2II1 = 0.
      tau_H2II2 = 0.
      tau_dust  = 0.
      nrm_dust  = 0.

      if (lgeom .eq. 3) then  ! RTP coordinates

      if (time .ge. t_on .and. time .lt. t_off) then

        do 30 n = 1, nnu1 + nnu2

         if (nu(n) .lt. 13.602) then 

#ifdef H2
          do 180 j = js  , je

                 i = is

           f(i,j,n)  = fcentral(n)

           tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
           tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
           if (iextinct .eq. 1) then
              tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                    abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
           endif

           tau       = tau_HM  + tau_H2II1 + tau_dust                       

           if (tau .le. taumin) then
             nioniz   = f(i,j,n) * tau * x1a(i) * x1a(i) 
     .                / (vol1a(i+1) - vol1a(i))
           else
             nioniz   = f(i,j,n) * (1.0-dexp(-1.0D0*tau))*x1a(i)*x1a(i) 
     .                / (vol1a(i+1) - vol1a(i))
           endif

           if (tau_HM .le. taumin) then
             nrm_HM  = tau_HM
           else
             nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
           endif

           if (tau_H2II1 .le. taumin) then
             nrm_H2II1  = tau_H2II1
           else
             nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
           endif

           if (iextinct .eq. 1) then
           if (tau_dust .le. taumin) then
             nrm_dust  = tau_dust
           else
             nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
           endif
           endif

           norm      = nioniz / ((nrm_HM  + nrm_H2II1 + nrm_dust)*mhi)

           k27(i,j)  = k27(i,j)+ nrm_HM   *norm/(abun(i,j,kslice,7)    )
           k28(i,j)  = k28(i,j)+ nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)

CDIR$ UNROLL UNROLL_I
           do 170 i = is+1, ie

            if (clight*time .lt. x1a(i)) goto 180
 
            f(i,j,n) = f(i-1,j,n) * dexp(-1.0D0*tau) * 
     .                 (x1a(i-1)*x1a(i-1))* (x1ai(i )*x1ai(i ))

            tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
            tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
            if (iextinct .eq. 1) then
               tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                     abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
            endif

            tau       = tau_HM  + tau_H2II1 + tau_dust                         

            if (tau .le. taumin) then
              nioniz   = f(i,j,n) * tau * x1a(i) * x1a(i) 
     .                 / (vol1a(i+1) - vol1a(i))
            else
              nioniz   = f(i,j,n) * (1.0-dexp(-1.0D0*tau))*x1a(i)*x1a(i) 
     .                 / (vol1a(i+1) - vol1a(i))
            endif

            if (tau_HM .le. taumin) then
              nrm_HM  = tau_HM
            else
              nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
            endif

            if (tau_H2II1 .le. taumin) then
              nrm_H2II1  = tau_H2II1
            else
              nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
            endif

            if (iextinct .eq. 1) then
            if (tau_dust .le. taumin) then
              nrm_dust  = tau_dust
            else
              nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
            endif
            endif

            norm     = nioniz / ((nrm_HM  + nrm_H2II1 + nrm_dust)*mhi)

            k27(i,j) = k27(i,j)+ nrm_HM   *norm/(abun(i,j,kslice,7)    )
            k28(i,j) = k28(i,j)+ nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)

170        continue
180       continue
#endif /* H2 */
         
         goto 30

         endif

         hnu = nu(n) * everg

         do 20 j = js  , je

                 i = is

          f(i,j,n) = fcentral(n)

#ifdef H
          tau_HI    = abun(i,j,kslice,1)     *  sigma24(n)*mhi*dx1a(i)
#endif /* H */
#ifdef He
          tau_HeI   = abun(i,j,kslice,4)*qrt *  sigma26(n)*mhi*dx1a(i)
          tau_HeII  = abun(i,j,kslice,5)*qrt *  sigma25(n)*mhi*dx1a(i) 
#endif /* He */
#ifdef H2
          tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
          tau_H2I   = abun(i,j,kslice,8)*haf *  sigma29(n)*mhi*dx1a(i) 
          tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
          tau_H2II2 = abun(i,j,kslice,9)*haf *  sigma30(n)*mhi*dx1a(i)
#endif /* H2 */
          if (iextinct .eq. 1) then
            tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                  abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
          endif

          tau      = 
#ifdef H
     .               tau_HI 
#endif /* H */
#ifdef He
     .             + tau_HeI + tau_HeII 
#endif /* He */
#ifdef H2
     .             + tau_HM  + tau_H2I + tau_H2II1 + tau_H2II2                         
#endif /* H2 */
     .             + tau_dust

          if (tau .le. taumin) then
            nioniz   = f(i,j,n) * tau * x1a(i) * x1a(i) 
     .               / (vol1a(i+1) - vol1a(i))
          else
            nioniz   = f(i,j,n) * (1.0-dexp(-1.0D0*tau))*x1a(i)*x1a(i) 
     .               / (vol1a(i+1) - vol1a(i))
	  endif

#ifdef H
          if (tau_HI .le. taumin) then
            nrm_HI   = tau_HI
          else
            nrm_HI   = 1.0-dexp(-1.0D0*tau_HI)
          endif
#endif /* H */

#ifdef He
          if (tau_HeI .le. taumin) then
            nrm_HeI  = tau_HeI
          else
            nrm_HeI  = 1.0-dexp(-1.0D0*tau_HeI)
          endif

          if (tau_HeII .le. taumin) then
            nrm_HeII = tau_HeII
          else
            nrm_HeII = 1.0-dexp(-1.0D0*tau_HeII)
          endif

#endif /* He */
#ifdef H2
          if (tau_HM .le. taumin) then
            nrm_HM  = tau_HM
          else
            nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
          endif

          if (tau_H2I .le. taumin) then
            nrm_H2I  = tau_H2I
          else
            nrm_H2I  = 1.0-dexp(-1.0D0*tau_H2I)
          endif

          if (tau_H2II1 .le. taumin) then
            nrm_H2II1  = tau_H2II1
          else
            nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
          endif

          if (tau_H2II2 .le. taumin) then
            nrm_H2II2  = tau_H2II2
          else
            nrm_H2II2  = 1.0-dexp(-1.0D0*tau_H2II2)
          endif
#endif /* H2 */

          if (iextinct .eq. 1) then
          if (tau_dust .le. taumin) then
            nrm_dust  = tau_dust
          else
            nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
          endif
          endif

          norm = nioniz / ((
#ifdef H
     .                      nrm_HI
#endif /* H */
#ifdef He
     .                    + nrm_HeI + nrm_HeII  
#endif /* He */
#ifdef H2
     .                    + nrm_HM  + nrm_H2I + nrm_H2II1 + nrm_H2II2
#endif /* H2 */
     .                    + nrm_dust
     .           )*mhi)

#ifdef H
          k24p        = nrm_HI*norm/abun(i,j,kslice,1)
          k24   (i,j) = k24   (i,j) + k24p
          k24mv (i,j) = k24mv (i,j) + k24p * hnu
          dum         = hnu         - e24  * everg
          piHI  (i,j) = piHI  (i,j) + dum  * k24p
#endif /* H */

#ifdef He
          k26p        = nrm_HeI*norm/(abun(i,j,kslice,4)*qrt)
          k26   (i,j) = k26   (i,j) + k26p
          k26mv (i,j) = k26mv (i,j) + k26p * hnu
          dum         = hnu         - e26  * everg
          piHeI (i,j) = piHeI (i,j) + dum  * k26p 
                      
          k25p        = nrm_HeII*norm/(abun(i,j,kslice,5)*qrt)
          k25   (i,j) = k25   (i,j) + k25p
          k25mv (i,j) = k25mv (i,j) + k25p * hnu
          dum         = hnu         - e25  * everg
          piHeII(i,j) = piHeII(i,j) + dum  * k25p                       
#endif /* He */

#ifdef H2
          k27   (i,j) = k27(i,j)+nrm_HM   *norm/(abun(i,j,kslice,7)    )
          k28   (i,j) = k28(i,j)+nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)
          k29   (i,j) = k29(i,j)+nrm_H2I  *norm/(abun(i,j,kslice,8)*haf)
          k30   (i,j) = k30(i,j)+nrm_H2II2*norm/(abun(i,j,kslice,9)*haf)
#endif /* H2 */

CDIR$ UNROLL UNROLL_I
         do 10 i = is+1, ie

          if (clight*time .lt. x1a(i)) goto 20
 
          f(i,j,n) = f(i-1,j,n) *dexp(-1.0D0*tau) * 
     .                       (x1a(i-1)*x1a(i-1))* (x1ai(i )*x1ai(i ))

#ifdef H
          tau_HI    = abun(i,j,kslice,1)     *  sigma24(n)*mhi*dx1a(i)
#endif /* H */
#ifdef He
          tau_HeI   = abun(i,j,kslice,4)*qrt *  sigma26(n)*mhi*dx1a(i)
          tau_HeII  = abun(i,j,kslice,5)*qrt *  sigma25(n)*mhi*dx1a(i) 
#endif /* He */
#ifdef H2
          tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
          tau_H2I   = abun(i,j,kslice,8)*haf *  sigma29(n)*mhi*dx1a(i) 
          tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
          tau_H2II2 = abun(i,j,kslice,9)*haf *  sigma30(n)*mhi*dx1a(i)
#endif /* H2 */
          if (iextinct .eq. 1) then
            tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                  abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
          endif

          tau      = 
#ifdef H
     .               tau_HI 
#endif /* H */
#ifdef He
     .             + tau_HeI + tau_HeII 
#endif /* He */
#ifdef H2
     .             + tau_HM  + tau_H2I + tau_H2II1 + tau_H2II2                         
#endif /* H2 */
     .             + tau_dust

          if (tau .le. taumin) then
            nioniz   = f(i,j,n) * tau * x1a(i) * x1a(i) 
     .               / (vol1a(i+1) - vol1a(i))
          else
            nioniz   = f(i,j,n) * (1.0-dexp(-1.0D0*tau))*x1a(i)*x1a(i) 
     .               / (vol1a(i+1) - vol1a(i))
          endif

#ifdef H
          if (tau_HI .le. taumin) then
            nrm_HI   = tau_HI
          else
            nrm_HI   = 1.0-dexp(-1.0D0*tau_HI)
          endif
#endif /* H */

#ifdef He
          if (tau_HeI .le. taumin) then
            nrm_HeI  = tau_HeI
          else
            nrm_HeI  = 1.0-dexp(-1.0D0*tau_HeI)
          endif

          if (tau_HeII .le. taumin) then
            nrm_HeII = tau_HeII
          else
            nrm_HeII = 1.0-dexp(-1.0D0*tau_HeII)
          endif

#endif /* He */
#ifdef H2
          if (tau_HM .le. taumin) then
            nrm_HM  = tau_HM
          else
            nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
          endif

          if (tau_H2I .le. taumin) then
            nrm_H2I  = tau_H2I
          else
            nrm_H2I  = 1.0-dexp(-1.0D0*tau_H2I)
          endif

          if (tau_H2II1 .le. taumin) then
            nrm_H2II1  = tau_H2II1
          else
            nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
          endif

          if (tau_H2II2 .le. taumin) then
            nrm_H2II2  = tau_H2II2
          else
            nrm_H2II2  = 1.0-dexp(-1.0D0*tau_H2II2)
          endif
#endif /* H2 */

          if (iextinct .eq. 1) then
          if (tau_dust .le. taumin) then
            nrm_dust  = tau_dust
          else
            nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
          endif
          endif

          norm = nioniz / ((
#ifdef H
     .                      nrm_HI
#endif /* H */
#ifdef He
     .                    + nrm_HeI + nrm_HeII  
#endif /* He */
#ifdef H2
     .                    + nrm_HM  + nrm_H2I + nrm_H2II1 + nrm_H2II2
#endif /* H2 */
     .                    + nrm_dust
     .           )*mhi)

#ifdef H
          k24p        = nrm_HI*norm/abun(i,j,kslice,1)
          k24   (i,j) = k24   (i,j) + k24p
          k24mv (i,j) = k24mv (i,j) + k24p * hnu
          dum         = hnu         - e24  * everg
          piHI  (i,j) = piHI  (i,j) + dum  * k24p
#endif /* H */

#ifdef He
          k26p        = nrm_HeI*norm/(abun(i,j,kslice,4)*qrt)
          k26   (i,j) = k26   (i,j) + k26p
          k26mv (i,j) = k26mv (i,j) + k26p * hnu
          dum         = hnu         - e26  * everg
          piHeI (i,j) = piHeI (i,j) + dum  * k26p 
                      
          k25p        = nrm_HeII*norm/(abun(i,j,kslice,5)*qrt)
          k25   (i,j) = k25   (i,j) + k25p
          k25mv (i,j) = k25mv (i,j) + k25p * hnu
          dum         = hnu         - e25  * everg
          piHeII(i,j) = piHeII(i,j) + dum  * k25p                       
#endif /* He */

#ifdef H2
          k27   (i,j) = k27(i,j)+nrm_HM   *norm/(abun(i,j,kslice,7)    )
          k28   (i,j) = k28(i,j)+nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)
          k29   (i,j) = k29(i,j)+nrm_H2I  *norm/(abun(i,j,kslice,8)*haf)
          k30   (i,j) = k30(i,j)+nrm_H2II2*norm/(abun(i,j,kslice,9)*haf)
#endif /* H2 */

10       continue
20       continue
30     continue

#ifdef H2
       if (iLW .eq. 1) then
       b5 = 9.12
       do 140 j = js  , je

         N_H2     = 0.
         tau_dust = 0.

         do 130 i = is  , ie

	    if (ibkgnd .eq. 1) then
		k31(i,j) = 1.13d08 * J_21 * 4. * pi
	    endif
            if (clight*time .lt. x1a(i)) goto 130
            N_H2     = N_H2 + dx1a(i) * abun(i,j,kslice,8) * mhi * haf
            if (iextinct .eq. 1) then
              tau_dust=tau_dust+(abun(i,j,kslice,1)+abun(i,j,kslice,2)+ 
     .                           abun(i,j,kslice,8))*mhi*xxi(38)*dx1a(i)
            endif
            x        = N_H2/5.0d14
            Fshield  = 0.965/(1 + (x/b5))**2 + 0.035/dsqrt(1+x)
     .               * dexp(-8.5d-4 * dsqrt(1+x))
            k31(i,j) = k31(i,j) + 1.13d08 * fnu * dexp(-1.0*tau_dust) *
     .                                      Fshield / (x1b(i)*x1b(i))

130     continue
140    continue
       endif
#endif /* H2 */

      endif

      endif  ! lgeom = 3 (RTP coordinates)

      if (lgeom .eq. 1 .or. lgeom .eq. 2) then  ! XYZ or ZRP coordinates

      if (time .ge. t_on .and. time .lt. t_off) then

        do 60 n = 1, nnu1 + nnu2

         if (nu(n) .lt. 13.602) then 

#ifdef H2
          do 160 j = js  , je

                 i = is

           f(i,j,n)  = fcentral(n)

           tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
           tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
           if (iextinct .eq. 1) then
             tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                   abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
           endif

           tau       = tau_HM  + tau_H2II1 + tau_dust                       

           if (tau .le. taumin) then
             nioniz  = f(i,j,n) * tau / dx1a(i)
           else
             nioniz  = f(i,j,n) * (1.0-dexp(-1.0D0*tau)) / dx1a(i)
           endif

           if (tau_HM .le. taumin) then
             nrm_HM  = tau_HM
           else
             nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
           endif

           if (tau_H2II1 .le. taumin) then
             nrm_H2II1  = tau_H2II1
           else
             nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
           endif

           if (iextinct .eq. 1) then
           if (tau_dust .le. taumin) then
             nrm_dust  = tau_dust
           else
             nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
           endif
           endif

           norm      = nioniz / ((nrm_HM  + nrm_H2II1 + nrm_dust)*mhi)

           k27(i,j)  = k27(i,j)+ nrm_HM   *norm/(abun(i,j,kslice,7)    )
           k28(i,j)  = k28(i,j)+ nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)

CDIR$ UNROLL UNROLL_I
           do 150 i = is+1, ie

            if (clight*time .lt. dabs(x1a(i) - x1a(is))) goto 160
 
            if (iPWA .eq. 0) then
             f(i,j,n) = f(i-1,j,n) * dexp(-1.0D0*tau) 
            else if (iPWA .eq. 1) then
             f(i,j,n) = f(i-1,j,n) * dexp(-1.0D0*tau) * 
     .       ((r_sep - xc + x1a(i-1)) * (r_sep - xc + x1a(i-1))) /
     .       ((r_sep - xc + x1a(i  )) * (r_sep - xc + x1a(i  ))) 
            endif

            tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
            tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
            if (iextinct .eq. 1) then
              tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                    abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
            endif

            tau       = tau_HM  + tau_H2II1 + tau_dust                         

            if (tau .le. taumin) then
              nioniz  = f(i,j,n) * tau / dx1a(i)
            else
              nioniz  = f(i,j,n) * (1.0-dexp(-1.0D0*tau)) / dx1a(i)
            endif

            if (tau_HM .le. taumin) then
              nrm_HM  = tau_HM
            else
              nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
            endif

            if (tau_H2II1 .le. taumin) then
              nrm_H2II1  = tau_H2II1
            else
              nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
            endif

            if (iextinct .eq. 1) then
            if (tau_dust .le. taumin) then
              nrm_dust  = tau_dust
            else
              nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
            endif
            endif

            norm      = nioniz / ((nrm_HM  + nrm_H2II1 + nrm_dust)*mhi)

            k27(i,j)  = k27(i,j)+nrm_HM   *norm/(abun(i,j,kslice,7)    )
            k28(i,j)  = k28(i,j)+nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)

150        continue
160       continue
#endif /* H2 */

         goto 60

         endif

         hnu = nu(n) * everg

         do 50 j = js  , je

                 i = is

          f(i,j,n) = fcentral(n)

#ifdef H
          tau_HI    = abun(i,j,kslice,1)     *  sigma24(n)*mhi*dx1a(i)
#endif /* H */ 
#ifdef He
          tau_HeI   = abun(i,j,kslice,4)*qrt *  sigma26(n)*mhi*dx1a(i)
          tau_HeII  = abun(i,j,kslice,5)*qrt *  sigma25(n)*mhi*dx1a(i) 
#endif /* He */
#ifdef H2
          tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
          tau_H2I   = abun(i,j,kslice,8)*haf *  sigma29(n)*mhi*dx1a(i) 
          tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
          tau_H2II2 = abun(i,j,kslice,9)*haf *  sigma30(n)*mhi*dx1a(i)
#endif /* H2 */
          if (iextinct .eq. 1) then
            tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                  abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
          endif

          tau      = 
#ifdef H
     .               tau_HI 
#endif /* H */
#ifdef He
     .             + tau_HeI + tau_HeII 
#endif /* He */
#ifdef H2
     .             + tau_HM  + tau_H2I + tau_H2II1 + tau_H2II2                         
#endif /* H2 */
     .             + tau_dust

          if (tau .le. taumin) then
            nioniz  = f(i,j,n) * tau / dx1a(i)
          else
            nioniz  = f(i,j,n) * (1.0-dexp(-1.0D0*tau)) / dx1a(i)
          endif

#ifdef H
          if (tau_HI .le. taumin) then
            nrm_HI  = tau_HI
          else
            nrm_HI  = 1.0-dexp(-1.0D0*tau_HI)
          endif
#endif /* H */

#ifdef He
          if (tau_HeI .le. taumin) then
            nrm_HeI  = tau_HeI
          else
            nrm_HeI  = 1.0-dexp(-1.0D0*tau_HeI)
          endif

          if (tau_HeII .le. taumin) then
            nrm_HeII = tau_HeII
          else
            nrm_HeII = 1.0-dexp(-1.0D0*tau_HeII)
          endif
#endif /* He */

#ifdef H2
          if (tau_HM .le. taumin) then
            nrm_HM  = tau_HM
          else
            nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
          endif

          if (tau_H2I .le. taumin) then
            nrm_H2I  = tau_H2I
          else
            nrm_H2I  = 1.0-dexp(-1.0D0*tau_H2I)
          endif

          if (tau_H2II1 .le. taumin) then
            nrm_H2II1  = tau_H2II1
          else
            nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
          endif

          if (tau_H2II2 .le. taumin) then
            nrm_H2II2  = tau_H2II2
          else
            nrm_H2II2  = 1.0-dexp(-1.0D0*tau_H2II2)
          endif

          if (iextinct .eq. 1) then
          if (tau_dust .le. taumin) then
            nrm_dust  = tau_dust
          else
            nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
          endif
          endif

#endif /* H2 */

          norm = nioniz / ((
#ifdef H
     .                      nrm_HI
#endif /* H */
#ifdef He
     .                    + nrm_HeI + nrm_HeII  
#endif /* He */
#ifdef H2
     .                    + nrm_HM  + nrm_H2I + nrm_H2II1 + nrm_H2II2
#endif /* H2 */
     .                    + nrm_dust
     .           )*mhi)

#ifdef H
          k24p        = nrm_HI*norm/abun(i,j,kslice,1)
          k24   (i,j) = k24   (i,j) + k24p
          dum         = hnu         - e24 * everg
          piHI  (i,j) = piHI  (i,j) + dum * k24p
#endif /* H */

#ifdef He
          k26p        = nrm_HeI*norm/(abun(i,j,kslice,4)*qrt)
          k26   (i,j) = k26   (i,j) + k26p
          dum         = hnu         - e26 * everg
          piHeI (i,j) = piHeI (i,j) + dum * k26p 
                      
          k25p        = nrm_HeII*norm/(abun(i,j,kslice,5)*qrt)
          k25   (i,j) = k25   (i,j) + k25p
          dum         = hnu         - e25 * everg
          piHeII(i,j) = piHeII(i,j) + dum * k25p                       
#endif /* He */

#ifdef H2
          k27   (i,j) = k27(i,j)+nrm_HM   *norm/(abun(i,j,kslice,7)    )
          k28   (i,j) = k28(i,j)+nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)
          k29   (i,j) = k29(i,j)+nrm_H2I  *norm/(abun(i,j,kslice,8)*haf)
          k30   (i,j) = k30(i,j)+nrm_H2II2*norm/(abun(i,j,kslice,9)*haf)
#endif /* H2 */

CDIR$ UNROLL UNROLL_I
         do 40 i = is+1, ie

          if (clight*time .lt. dabs(x1a(i) - x1a(is))) goto 50
 
          if (iPWA .eq. 0) then
            f(i,j,n) = f(i-1,j,n) * dexp(-1.0D0*tau) 
          else if (iPWA .eq. 1) then
            f(i,j,n) = f(i-1,j,n) * dexp(-1.0D0*tau) * 
     .      ((r_sep - xc + x1a(i-1)) * (r_sep - xc + x1a(i-1))) /
     .      ((r_sep - xc + x1a(i  )) * (r_sep - xc + x1a(i  ))) 
          endif
#ifdef H
          tau_HI    = abun(i,j,kslice,1)     *  sigma24(n)*mhi*dx1a(i)
#endif /* H */
#ifdef He
          tau_HeI   = abun(i,j,kslice,4)*qrt *  sigma26(n)*mhi*dx1a(i)
          tau_HeII  = abun(i,j,kslice,5)*qrt *  sigma25(n)*mhi*dx1a(i) 
#endif /* He */
#ifdef H2
          tau_HM    = abun(i,j,kslice,7)     *  sigma27(n)*mhi*dx1a(i) 
          tau_H2I   = abun(i,j,kslice,8)*haf *  sigma29(n)*mhi*dx1a(i) 
          tau_H2II1 = abun(i,j,kslice,9)*haf *  sigma28(n)*mhi*dx1a(i)
          tau_H2II2 = abun(i,j,kslice,9)*haf *  sigma30(n)*mhi*dx1a(i)
#endif /* H2 */
          if (iextinct .eq. 1) then
            tau_dust = (abun(i,j,kslice,1) + abun(i,j,kslice,2) + 
     .                  abun(i,j,kslice,8)) * mhi * xxi(n) *  dx1a(i)
          endif

          tau      = 
#ifdef H
     .               tau_HI 
#endif /* H */
#ifdef He
     .             + tau_HeI + tau_HeII 
#endif /* He */
#ifdef H2
     .             + tau_HM  + tau_H2I + tau_H2II1 + tau_H2II2                         
#endif /* H2 */
     .             + tau_dust

          if (tau .le. taumin) then
            nioniz   = f(i,j,n) * tau / dx1a(i)
          else 
            nioniz   = f(i,j,n) * (1.0-dexp(-1.0D0*tau)) / dx1a(i)
          endif

#ifdef H
          if (tau_HI .le. taumin) then
            nrm_HI   = tau_HI
          else
            nrm_HI   = 1.0-dexp(-1.0D0*tau_HI)
          endif
#endif /* H */

#ifdef He
          if (tau_HeI .le. taumin) then
            nrm_HeI  = tau_HeI
          else
            nrm_HeI  = 1.0-dexp(-1.0D0*tau_HeI)
          endif

          if (tau_HeII .le. taumin) then
            nrm_HeII = tau_HeII
          else
            nrm_HeII = 1.0-dexp(-1.0D0*tau_HeII)
          endif
#endif /* He */

#ifdef H2
          if (tau_HM .le. taumin) then
            nrm_HM  = tau_HM
          else
            nrm_HM  = 1.0-dexp(-1.0D0*tau_HM)
          endif

          if (tau_H2I .le. taumin) then
            nrm_H2I  = tau_H2I
          else
            nrm_H2I  = 1.0-dexp(-1.0D0*tau_H2I)
          endif

          if (tau_H2II1 .le. taumin) then
            nrm_H2II1  = tau_H2II1
          else
            nrm_H2II1  = 1.0-dexp(-1.0D0*tau_H2II1)
          endif

          if (tau_H2II2 .le. taumin) then
            nrm_H2II2  = tau_H2II2
          else
            nrm_H2II2  = 1.0-dexp(-1.0D0*tau_H2II2)
          endif
#endif /* H2 */

          if (iextinct .eq. 1) then
          if (tau_dust .le. taumin) then
            nrm_dust  = tau_dust
          else
            nrm_dust  = 1.0-dexp(-1.0D0*tau_dust)
          endif
          endif

          norm = nioniz / ((
#ifdef H
     .                      nrm_HI
#endif /* H */
#ifdef He
     .                    + nrm_HeI + nrm_HeII 
#endif /* He */
#ifdef H2
     .                    + nrm_HM  + nrm_H2I + nrm_H2II1 + nrm_H2II2
#endif /* H2 */
     .                    + nrm_dust
     .           )*mhi)

#ifdef H
          k24p        = nrm_HI*norm/abun(i,j,kslice,1)
          k24   (i,j) = k24   (i,j) + k24p
          dum         = hnu         - e24 * everg
          piHI  (i,j) = piHI  (i,j) + dum * k24p
#endif /* H */

#ifdef He
          k26p        = nrm_HeI*norm/(abun(i,j,kslice,4)*qrt)
          k26   (i,j) = k26   (i,j) + k26p
          dum         = hnu         - e26 * everg
          piHeI (i,j) = piHeI (i,j) + dum * k26p 
                      
          k25p        = nrm_HeII*norm/(abun(i,j,kslice,5)*qrt)
          k25   (i,j) = k25   (i,j) + k25p
          dum         = hnu         - e25 * everg
          piHeII(i,j) = piHeII(i,j) + dum * k25p                       
#endif /* He */

#ifdef H2
          k27   (i,j) = k27(i,j)+nrm_HM   *norm/(abun(i,j,kslice,7)    )
          k28   (i,j) = k28(i,j)+nrm_H2II1*norm/(abun(i,j,kslice,9)*haf)
          k29   (i,j) = k29(i,j)+nrm_H2I  *norm/(abun(i,j,kslice,8)*haf)
          k30   (i,j) = k30(i,j)+nrm_H2II2*norm/(abun(i,j,kslice,9)*haf)

#endif /* H2 */
40       continue
50       continue
60     continue
#ifdef H2
       if (iLW .eq. 1) then
       b5 = 9.12
       do 100 j = js  , je


         N_H2 = 0.0
         tau_dust = 0.

         do 90 i = is  , ie

	    if (ibkgnd .eq. 1) then
		k31(i,j) = 1.13d08 * J_21 * 4. * pi
	    endif
            if (clight*time .lt. x1a(i)) goto 90
            N_H2     = N_H2 + dx1a(i) * abun(i,j,kslice,8) * mhi * haf
            if (iextinct .eq. 1) then
              tau_dust=tau_dust+(abun(i,j,kslice,1)+abun(i,j,kslice,2)+ 
     .                           abun(i,j,kslice,8))*mhi*xxi(38)*dx1a(i)
            endif
            x        = N_H2/5.0d14
            Fshield  = 0.965/(1 + (x/b5))**2 + 0.035/dsqrt(1+x)
     .               * dexp(-8.5d-4 * dsqrt(1+x))
            if (iPWA .eq. 0) then
             sep = r_sep
            else if (iPWA .eq. 1) then
             sep = r_sep - xc + x1a(i)
            endif
            k31(i,j) = k31(i,j) + 1.13d08 * fnu * dexp(-1.0*tau_dust) * 
     .                                      Fshield / (sep * sep)

90       continue
100    continue
       endif
#endif /* H2 */

      endif  !  t_on < t < t_off

      endif  !  lgeom = 1 or 2 (XYZ or ZRP coordinates)
 2000 return
      end


c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////           M N U _ D I R _ F L X           \\\\\\\\\\
c
c=======================================================================
c
c
