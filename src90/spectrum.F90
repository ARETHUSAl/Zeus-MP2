#include "rtchem.def"
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              S P E C T R U M              \\\\\\\\\\
!
!=======================================================================
!
       subroutine spectrum 
!  written by: Daniel Whalen 06.09.07
!
!  PURPOSE:  Generate radiation spectrum for boundary fluxes.
!
!  BOUNDARY VALUES USED:
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
      integer  :: i, j, kslice, n, index, niter, lbin
      real(rl) :: kt, LSnm,delta,u_bound,l_bound,bnu,delta_u,u,fLW,sep
      real(rl) :: ndot,d_l
      real(rl), dimension(40) :: bin
#ifndef UNROLL_I
#define UNROLL_I
#endif
!
!-----------------------------------------------------------------------
!
      if (ispct .eq. 2) then
      kt   = t_star/tevk
!        if (nspec .gt. 6) then
! Compute the number of photons emitted by the blackbody with 0.755 eV
! < hnu < 13.602, nphdot1.  Use this number to normalize the blackbody 
! number density as a function of energy over this range.  We need this 
! spectrum to evaluate the H- photodetachment rate coefficient k27 as  
! well as the H2+ photodissociation rate k28 below the H ionization edge
!          delta_u = 0.05
!          LSnm    = 0.0
!          do n=1,2000
!            u     = n * delta_u
!            LSnm  = LSnm  + u**3/(dexp(u)-1) * delta_u
!          enddo
!          niter   = 1000000
!          l_bound = (e27 * tevk)/t_star  
!          u_bound = (e24 * tevk)/t_star  
!          delta_u = (u_bound - l_bound)/niter
!          delta   = 0.0
!          do n=1,niter+1
!            u     = l_bound + (n-1) * delta_u
!            delta = delta + u**2/(dexp(u)-1) * delta_u
!          enddo
!          nphdot1 = (L_star * delta)/(kt * everg * LSnm)
!          bbnm1   = 0.0
!          do n=1,nnu1
!            u     = (nu(n) * everg)/(hplanck * clight)
!            bnu   = u**2/(dexp(nu(n)/kt)-1)
!            bbnm1 = bbnm1 + bnu 
!          enddo
!          bbnm1 = (nphdot1 * nsrc) / (4.0 * pi * bbnm1)
!        endif  ! nspec > 6
! Normalize the blackbody ionizing photon number spectrum according  
! to the total number of photons above the Lyman limit exiting the 
! star with temperature t_star per second, which is nphdot (read in
! from zmp_inp and taken from Schaerer, et al 2002).
        bbnm2   = 0.0
        do n=nnu1+1,nnu1+nnu2
          u     = (nu(n) * everg)/(hplanck * clight)
          bnu   = u**2/(dexp(nu(n)/kt)-1)
          bbnm2 = bbnm2 + bnu 
        enddo
        bbnm2 = (nphdot * nsrc) / (4.0 * pi * bbnm2)
!        if (nspec .gt. 6) then
! Calculate the mean intensity J_nu across the Lyman-Werner band
! to compute the H2 photodissociation rate k31 in mnu_dir_flx.F
!          if (iLW .eq. 1) then
!            delta_u = 0.1
!            LSnm    = 0.0
!            do n=1,1000
!              u    = n * delta_u
!              LSnm = LSnm + u**3/(dexp(u)-1) * delta_u
!            enddo
!            niter   = 1000000
!            l_bound = (11.18 * tevk)/t_star  
!            u_bound = (13.6  * tevk)/t_star  
!            delta_u = (u_bound - l_bound)/niter
!            fLW  = 0.0
!            do n=1,niter+1
!              u = l_bound + (n-1) * delta_u
!              bnu = u**3/(dexp(u)-1)
!              fLW = fLW + bnu * delta_u
!            enddo
!            fnu = (L_star*hplanck)/((13.602-11.18)*everg*4.*pi)
!     .          * (fLW/LSnm)
!          endif   ! iLW
!        endif     ! nspec > 6
      endif       ! ispct=2
      if (ispct .eq. 3) then
        qsonm = 0.0
        do n=nnu1+1,nnu1+nnu2
           qsonm = qsonm + (hplanck/(nu(n)*everg))**(alpha+1.)
        enddo
        qsonm = (nphdot * nsrc) / (4.0 * pi * qsonm)
      endif       ! ispct=3
      if (lgeom .eq. 3) then  ! RTP coordinates
        sep = x1a(is)
      else if (lgeom .eq. 1 .or. lgeom .eq. 2) then  ! XYZ or ZRP coordinates
        if (iPWA .eq. 0) then
          sep = r_sep
        else if (iPWA .eq. 1) then
          sep = r_sep - xc + x1a(is)
        endif
      endif     ! lgeom
      if (ispct .eq. 1) then
      !write(*,*) 'nphdot', nphdot
      !write(*,*) 'nsrc', nsrc
      !write(*,*) 'pi', pi
      !write(*,*) 'sep', sep
      !write(*,*) 'fcentral', size(fcentral)
       fcentral(1) = (nphdot * nsrc) / (4 * pi * sep * sep)
      else if (ispct .eq. 2) then
       if (nspec .gt. 6) then
        fnu = 0.
        do n = 1, nnu1 
          fcentral(n) = bbnm2*((nu(n) * everg) / (hplanck*clight))**2 &
                      / ((dexp(nu(n)/kt) - 1.) * sep * sep)
          if (nu(n) .ge. 11.18 .and. nu(n) .lt. 13.6 .and. iLW .eq. 1) &
          fnu = fnu + fcentral(n) * sep * sep * 0.5 * (nu(n)+nu(n+1)) &
              * everg
        enddo
        fnu = fnu * hplanck /((13.602-11.18)*everg)
       endif  ! nspec > 6 
       do n = nnu1 + 1, nnu1 + nnu2
         fcentral(n) = bbnm2*((nu(n) * everg)/(hplanck*clight))**2 &
                     / ((dexp(nu(n)/kt) - 1.) * sep*sep)
       enddo 
      else if (ispct .eq. 3) then
       if (nspec .gt. 6) then
        fnu = 0.
        do n = 1, nnu1 
          fcentral(n) = qsonm * (hplanck/nu(n)*everg)**(alpha+1.) &
                      / (sep * sep)
          if (nu(n) .ge. 11.18 .and. nu(n) .lt. 13.6 .and. iLW .eq. 1) &
          fnu = fnu + fcentral(n) * sep * sep * 0.5 * (nu(n)+nu(n+1)) &
              * everg
        enddo
        fnu = fnu * hplanck /((13.602-11.18)*everg)
       endif  ! nspec > 6 
       do n = nnu1 + 1, nnu1 + nnu2
          fcentral(n) = qsonm * (hplanck/nu(n)*everg)**(alpha+1.) &
                      / (sep * sep)
       enddo 
      else if (ispct .eq. 4) then
       d_l    = 5.8d28
       rdshft = 2.328
       if (time .lt. spect(1,1)) then
         t_lum = spect(1,1) - time
         do i=3,nnu1+nnu2+2
           fcentral(i-2) = (spect(1,i) * nsrc) / (4 * pi * sep * sep)
!           print*,"nphdot is: ", nu(i-2), fcentral(i-2)*(4*pi*sep*sep)
         enddo
         goto 20
       endif
       do n = 1,nspctpts-1
          if (time .ge. spect(n  ,1) .and. &
              time .lt. spect(n+1,1)) then
            lbin = n
            goto 10
          endif
       enddo
 10    continue
       t_lum = spect(lbin+1,1) - time
!       print*,"lbin, time, and t_lum are: ",lbin, time, t_lum
       do i=3,nnu1+nnu2+2
         ndot = spect(lbin,i)  + (spect(lbin+1,i) - &
                spect(lbin,i)) / (spect(lbin+1,1) - &
                spect(lbin,1)) * (time           - &
                spect(lbin,1)) 
         fcentral(i-2) = (ndot * nsrc) / (4 * pi * sep * sep)
!         print*,nu(i-2),spect(lbin,i),spect(lbin+1,i),spect(lbin+1,1),
!     .   spect(lbin,1), time-spect(lbin,1)
!         print*,"nphdot is: ", i, nu(i-2), fcentral(i-2)*(4*pi*sep*sep)
       enddo
 20    continue
       if (iLW .eq. 1 .or. nspec .gt. 6) then
         fnu = 0.
         do n=1,nnu1
           if (nu(n) .ge. 11.18 .and. nu(n) .lt. 13.6) &
           fnu = fnu + fcentral(n) * sep * sep * 0.5 * (nu(n)+nu(n+1)) &
               * everg
         enddo
         fnu = fnu * hplanck /((13.602-11.18)*everg)
!         print*,"before fireball, fnu is: ",fnu
         if (time .le. 691.) then
           fnu = fnu + 3.07d-27 * (d_l * d_l * 1.2311d-33) &
               * (1+rdshft)**(-0.75) * ((time+3600.)*1.754d-3)**(-0.8)
         else
           fnu = fnu + 6.46d-29 * (d_l * d_l * 1.2311d-33) &
               * (1+rdshft)**(-0.67) * ((time+3600.)*8.448d-6)**(-0.72)
         endif  ! fireball
!         print*,"after fireball, fnu is: ",
!     .         3.07d-27 * (d_l * d_l * 1.2311d-33) 
!     .         * (1+rdshft)**(-0.75) * ((time+3600.)*1.754d-3)**(-0.8)
!     .         * sep * sep
       endif    ! iLW
      endif     ! ispct
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////              S P E C T R U M              \\\\\\\\\\
!
!=======================================================================
!
!
