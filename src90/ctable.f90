!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                C T A B L E                \\\\\\\\\!
!=======================================================================
!
      subroutine ctable
!
!  adapted from Zeus 3.4 03.01.01 by djw
!
!  updated 12.29.05 to properly implement Dalgarno McCray (1972) ISM
!  cooling
!
!  upgraded 06.30.06 to support multifrequency/multigroup RT 
!
!  ported to ZEUS-MP 2.1 on 11-27-06
!
!  reaction network and RT upgrades begun on 12-01-12
!
!  PURPOSE:          Construct rate coefficient lookup tables for the
!                    cooling and chemistry modules of the multi-species
!                    algorithm.  Both primordial gas cooling (Anninos
!                    et al 1997) and galactic ISM cooling (Dalgarno &
!                    McCray 1972) are available.
! 
! 9-species primordial chemistry rate coefficients
!                   
! ------- k1:    H    + e  -> H+   + 2e
! ------- k2:    H+   + e  -> H    + gamma
! ------- k3:    He   + e  -> He+  + 2e
! ------- k4:    He+  + e  -> He   + gamma
! ------- k5:    He+  + e  -> He++ + 2e
! ------- k6:    He++ + e  -> He+  + gamma
! ------- k7:    H    + e  -> H-   + gamma
! ------- k8:    H-   + H  -> H2*  + e
! ------- k9:    H    + H+ -> H2+  + gamma
! ------- k10:   H2+  + H  -> H2*  + H+
! --------k11:   H2   + H+ -> H2+  + H
! --------k12:   H2   + e  -> 2H   + e
! --------k13:   H2   + H  -> 3H
! --------k14:   H-   + e  -> H    + 2e
! --------k15:   H-   + H  -> 2H   + e
! --------k16:   H-   + H+ -> 2H
! --------k17:   H-   + H+ -> H2+  + e
! --------k18:   H2+  + e  -> 2H
! --------k19:   H2+  + H- -> H    + H2
! --------k22:   2H   + H  -> H2I  + H
! --------k24:   H    + gamma  -> H+   + e
! --------k25:   He+  + gamma  -> He++ + e
! --------k26:   He   + gamma  -> He+  + e
! --------k27:   H-   + gamma  -> H    + e
! --------k28:   H2+  + gamma  -> H    + H+
! --------k29:   H2   + gamma  -> H2+  + e
! --------k30:   H2+  + gamma  -> 2H+  + e
! --------k31:   H2   + gamma  -> 2H
!       ------ Deuterium rates -----
! --------k32:   D    + gamma  -> D+   + e (set equal to k31)
! --------k33:   HD   + gamma  -> H    + D
! --------k50:   H+  +  D  -> H  +  D+
! --------k51:   H   +  D+ -> H+ +  D
! --------k52:   H2  +  D+ -> HD +  H+
! --------k53:   HD  +  H+ -> H2 +  D+
! --------k54:   H2  +  D  -> HD +  H
! --------k55:   HD  +  H  -> H2 +  D
! --------k56:   D   +  H- -> HD +  e-
! -------[k57:   D-  +  H  -> HD +  e-]  included by multiplying 56 by 2
!
! Also, for 3-body production of HD above n = 10^8 cm^-3, we include
!
! ----0.5*k22:   2H  +  D  -> HD +  H
! ------?----:  H + D + H2 -> HD +  H2
!
!  LOCAL VARIABLES:
!
!  BOUNDARY VALUES USED:
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
!
      use cons
      use config
      use param
      use chem
      use mpiyes
      use root
      use field
      use grid
      implicit NONE
      integer  ::  i, j, k, l, niter, ntau, m, nHpts, nMpts, &
                   nDMpts, nHIIpts, nHeIIpts, nHeIIIpts, index
      real(rl) ::  ttt      , logtev  , logttt  ,tev    , &
                   kt       , xx      , dum     ,logt4  , &
                   log10ttt , t_half  , tcev    ,bnu    , &
                   u        , delta_u , dlognu  ,hnumin , &
                   LSnm     , fLW     , l_bound ,u_bound, &
                   tm       , lt      , t3      ,HDLR   , &
                   HDLV     , tau0    , delta   ,hc     , &
                   gtd      , lt3
      real(rl) ::  O      , C       , N      , Si      , Fe     , S &
      ,           Le_CI   , Le_SiI  , Le_FeI , Le_O    , Le_H   , Lem_O &
      ,           Lem_1FeI, Lem_2FeI, Lem_OI , Lem_CI  , Lem_SI , Lem_N &
      ,           Lem_SiI , LH_CI   , LH_O   , LH_C    , LH_SiI , LH_FeI &
      ,           t_mhalf , t_phalf , ifrac  , logtmin , logtmax, dlogt2 &
      ,           t_upper , t_lower , lgtmax2, lgtmin2 , logt_low 
      real(rl) ::  Hcool (18), cool1  (28), cool2  (28), cool3 (28), &
                   cool4(308), HIIBr  (31), HIIBc  (31), HeIIBr(18), &
                   HeIIBc(18), HeIIIBr(31), HeIIIBc(31), lmbda (6 ), &
                   a     ( 6), b      ( 6), c1     ( 6)
      temstart  = 10
      temend    = 1.0d09
      nHpts     = 18
      nMpts     = 28
      nDMpts    = 308
      nHIIpts   = 31
      nHeIIpts  = 18
      nHeIIIpts = 27
      logtmin = dlog(temstart)
      logtmax = dlog(temend)
      dlogtem = (logtmax-logtmin)/dfloat(nratec-1)
      lgtmin2 = 9.37371      ! temperature range of the upper			
      lgtmax2 = 15.8644      ! portion of the DM cooling curve
      dlogt2  = (lgtmax2 - lgtmin2)/dfloat(nDMpts-1)
      open(unit=10, file='LeH.dat'      , status='unknown')
      open(unit=11, file='newLH.dat'    , status='unknown')
      open(unit=12, file='linecool2.dat', status='unknown')
      open(unit=13, file='HII.dat'      , status='unknown')
      open(unit=14, file='HeII.dat'     , status='unknown')
      open(unit=15, file='HeIII.dat'    , status='unknown')
      if (ispct .eq. 4) then
        open(unit=16, file='spectrum.dat' , status='unknown')
        read(16,*) nspctpts
        read(16,*)
        allocate(spect(nspctpts,nnu1+nnu2+3))
        do i=1,nspctpts
           read(16,*) (spect(i,j), j=1,nnu1+nnu2+3)
        enddo
        close(unit=16)
      endif
      do i=1,nHpts
         read(10,*) ttt,Hcool(i)
         Hcool(i) = dlog(Hcool(i))
      enddo
      do i=1,nMpts
         read(11,*) ttt,cool1(i),cool2(i),cool3(i)
         cool1(i) = dlog(cool1(i))
         cool2(i) = dlog(cool2(i))
         cool3(i) = dlog(cool3(i))
      enddo
      do i=1,nDMpts
         read(12,*) logttt,cool4(i)
      enddo
      do i=1,nHIIpts
         read(13,*) HIIBr(i),HIIBc(i)
      enddo
      do i=1,nHeIIpts
         read(14,*) HeIIBr(i),HeIIBc(i)
      enddo
      do i=1,nHeIIIpts
         read(15,*) HeIIIBr(i),HeIIIBc(i)
      enddo
      close(unit=10)
      close(unit=11)
      close(unit=12)
      close(unit=13)
      close(unit=14)
      close(unit=15)
      ifrac   = 0.0001
      O      = 4.40d-4 * z_sol 
      C      = 3.75d-4 * z_sol
      N      = 8.70d-5 * z_sol
      Si     = 3.20d-5 * z_sol
      Fe     = 3.20d-5 * z_sol
      S      = 1.40d-5 * z_sol
!
!  Construct tables for rate coefficients and some cooling functions 
!  versus the logarithm of temperature.
!
      do i = 1, nratec
       if (nspec .gt. 1) then
        k1a     (i) = 0.0
        k2a     (i) = tiny
        ceHIa   (i) = tiny
        ciHIa   (i) = tiny
        reHIIa  (i) = tiny
        brema   (i) = tiny
        compt   (i) = tiny
       endif
       if (nspec .gt. 3) then
        k3a     (i) = tiny
        k4a     (i) = tiny
        k5a     (i) = tiny
        k6a     (i) = tiny
        ceHeIa  (i) = tiny
        ceHeIIa (i) = tiny
        ciHeIa  (i) = tiny
        ciHeISa (i) = tiny
        ciHeIIa (i) = tiny
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
       endif
       if (nspec .gt. 6) then
        k7a     (i) = tiny
        k8a     (i) = tiny
        k9a     (i) = tiny
        k10a    (i) = tiny
        k11a    (i) = tiny
        k12a    (i) = tiny
        k13a    (i) = tiny
        k14a    (i) = tiny
        k15a    (i) = tiny
        k16a    (i) = tiny
        k17a    (i) = tiny
        k18a    (i) = tiny
        k19a    (i) = tiny
        k22a    (i) = tiny
        hyd01ka (i) = tiny
        h2k01a  (i) = tiny
        vibha   (i) = tiny
        rotha   (i) = tiny
        rotla   (i) = tiny
        gpldla  (i) = tiny
        gphdla  (i) = tiny
       endif
       if (nspec .gt. 9) then
        k50a    (i) = tiny
        k51a    (i) = tiny
        k52a    (i) = tiny
        k53a    (i) = tiny
        k54a    (i) = tiny
        k55a    (i) = tiny
        k56a    (i) = tiny
        hdltea  (i) = tiny
        hdlowa  (i) = tiny
       endif
      enddo
! fill in tables over the range temstart to temend
      do i = 1, nratec
        logttt = logtmin + dfloat(i-1)*dlogtem
        ttt = dexp(logttt)
        tev = ttt/tevk
        logtev = dlog(tev)
        logt4  = dlog(tev/4.0)
        t_half = dsqrt(ttt)
        if (tev .gt. 0.8) then
         if (nspec .gt. 1) then
          k1a(i) = dexp(-32.71396786375 &
                 + 13.53655609057*logtev &
                 - 5.739328757388*logtev**2 &
                 + 1.563154982022*logtev**3 &
                 - 0.2877056004391*logtev**4 &
                 + 0.03482559773736999*logtev**5 &
                 - 0.00263197617559*logtev**6 &
                 + 0.0001119543953861*logtev**7 &
                 - 2.039149852002d-6*logtev**8)
!          k1a(i) = 1.08d-10 * (ttt**(0.5)) *dexp(-1.5789e05/ttt)   !FTB90 benchmarks
         endif
         if (nspec .gt. 3) then
          k3a(i) = dexp(-44.09864886561001 &
                 + 23.91596563469*logtev &
                 - 10.75323019821*logtev**2 &
                 + 3.058038757198*logtev**3 &
                 - 0.5685118909884001*logtev**4 &
                 + 0.06795391233790001*logtev**5 &
                 - 0.005009056101857001*logtev**6 &
                 + 0.0002067236157507*logtev**7 &
                 - 3.649161410833d-6*logtev**8)
          if (iOTS .eq. 0) then
            k4a(i) = 1.54d-9*(1.+0.3/dexp(8.099328789667/tev)) &
                 / (dexp(40.49664394833662/tev)*tev**1.5) &
                 + 3.92d-13/tev**0.6353
          endif
          if (iOTS .eq. 1) then
            log10ttt  = dlog10(ttt)
            if (ttt .ge. 10. .and. ttt .le. 2.5d04) then
             index     = min0(nHeIIpts-1,max0(0, &
                         idint((log10ttt-1.0)/0.2)+1))
             logt_low  = 1. + dfloat(index-1) * 0.2
             k4a(i)    = 10.0**(HeIIBr(index  ) + (log10ttt - logt_low)* &
                               (HeIIBr(index+1) -   HeIIBr(index))/0.2)/ &
                                t_half
            endif
          endif
          k5a(i) = dexp(-68.71040990212001 &
                 + 43.93347632635*logtev &
                 - 18.48066993568*logtev**2 &
                 + 4.701626486759002*logtev**3 &
                 - 0.7692466334492*logtev**4 &
                 + 0.08113042097303*logtev**5 &
                 - 0.005324020628287001*logtev**6 &
                 + 0.0001975705312221*logtev**7 &
                 - 3.165581065665d-6*logtev**8)
         endif ! nspec > 3
        else   ! tev   > 0.8
         if (nspec .gt. 1) then
          k1a(i) = tiny
         endif
         if (nspec .gt. 3) then
          k3a(i) = tiny
          if (iOTS .eq. 0) then
           k4a(i) = 3.92d-13/tev**0.6353
          endif
          if (iOTS .eq. 1) then
           log10ttt  = dlog10(ttt)
           if (ttt .ge. 10. .and. ttt .le. 2.5d04) then
             index     = min0(nHeIIpts-1,max0(0, &
                         idint((log10ttt-1.0)/0.2)+1))
             logt_low  = 1. + dfloat(index-1) * 0.2
             k4a(i)    = 10.0**(HeIIBr(index  ) + (log10ttt - logt_low)* &
                               (HeIIBr(index+1) -   HeIIBr(index))/0.2)/ &
                                t_half   
           endif
          endif
          k5a(i) = tiny
         endif ! nspec > 3
       endif   ! tev   > 0.8
      if (nspec .gt. 1) then
! Anninos, et al 1997 H case A recombination coefficients--incompatible
! with the on-the-spot approximation:
        if (iOTS .eq. 0) then
         if ( ttt .gt. 5500.0 ) then
          k2a(i) = dexp(-28.61303380689232 &
                 - 0.7241125657826851*logtev &
                 - 0.02026044731984691*logtev**2 &
                 - 0.002380861877349834*logtev**3 &
                 - 0.0003212605213188796*logtev**4 &
                 - 0.00001421502914054107*logtev**5 &
                 + 4.989108920299513d-6*logtev**6 &
                 + 5.755614137575758d-7*logtev**7 &
                 - 1.856767039775261d-8*logtev**8 &
                 - 3.071135243196595d-9*logtev**9)
         else
          k2a(i) = k4a(i)
         endif
        endif
! H case B recombination coefficients: 
        if (iOTS .eq. 1) then
!          k2a(i) = tiny                     ! used in test 1 of WN05
!          k2a(i) = 2.4773d-13               ! used in some of the Whalen,
                                             ! et al 2005 method paper tests
!          k2a(i) = 2.015d-13                ! static rjw test for z_str
!          k2a(i) = 3.0d-10 * (ttt)**(-0.75) ! Tenorio-Tagle, et al 1986
                                             ! case B recombination coeff  
!          if (usrtag .eq. 'T1_') then       ! Toronto code tests case B
!            k2a(i) = 2.59d-13               ! recombination coefficient 
!          else if (usrtag .eq. 'T2_' .or. 
!     .             usrtag .eq. 'T5_' .or. usrtag .eq. 'T6_') then
!            k2a(i) = 2.59d-10 * (ttt)**(-0.75) 
!          endif
          log10ttt  = dlog10(ttt)
          if (ttt .ge. 10. .and. ttt .le. 1.0d07) then
           index     = min0(nHIIpts-1,max0(0, &
                       idint((log10ttt-1.0)/0.2)+1))
           logt_low  = 1. + dfloat(index-1) * 0.2
           k2a(i)    = 10.0**(HIIBr(index  ) + (log10ttt - logt_low) * &
                             (HIIBr(index+1) -    HIIBr(index))/0.2) / &
                              t_half   
          endif
        endif
      endif   ! nspec > 1
      if (nspec .gt. 3) then
        if (iOTS .eq. 0) then
         k6a(i)= 3.36d-10/dsqrt(ttt)/(ttt/1.e3)**0.2/(1+(ttt/1.e6)**0.7)
        endif
        if (iOTS .eq. 1) then
          log10ttt  = dlog10(ttt)
          if (ttt .ge. 40. .and. ttt .le. 1.0d07) then
           index     = min0(nHeIIIpts-1,max0(0, &
                       idint((log10ttt-1.6)/0.2)+1))
           logt_low  = 1.6 + dfloat(index-1) * 0.2
           k6a(i)    = 10.0**(HeIIIBr(index  ) + (log10ttt - logt_low) * &
                             (HeIIIBr(index+1) -  HeIIIBr(index))/0.2) / &
                              t_half   
          endif
        endif
!          k6a(i) = 2.0*dexp(-28.61303380689232
!     .           - 0.7241125657826851*logt4
!     .           - 0.02026044731984691*logt4**2
!     .           - 0.002380861877349834*logt4**3
!     .           - 0.0003212605213188796*logt4**4
!     .           - 0.00001421502914054107*logt4**5
!     .           + 4.989108920299513d-6*logt4**6
!     .           + 5.755614137575758d-7*logt4**7
!     .           - 1.856767039775261d-8*logt4**8
!     .           - 3.071135243196595d-9*logt4**9)
      endif
      if (nspec .gt. 6) then
        k7a(i) = 6.77d-15*tev**0.8779
!        if (ttt .le. 1.5d4) then                   ! Shapiro & Kang
!          k7a(i) = 1.0d-18 * ttt
!        else if (ttt .gt. 1.5d4) then
!          k7a(i) = 10.**(-14.10+0.1175*dlog10(ttt)
!     .           - 9.813d-3*(dlog10(ttt))**2.)
!        endif
        if (tev .gt. 0.1) then
          k8a(i) = dexp(-20.06913897587003 &
                 + 0.2289800603272916*logtev &
                 + 0.03599837721023835*logtev**2 &
                 - 0.004555120027032095*logtev**3 &
                 - 0.0003105115447124016*logtev**4 &
                 + 0.0001073294010367247*logtev**5 &
                 - 8.36671960467864d-6*logtev**6 &
                 + 2.238306228891639d-7*logtev**7)
        else
          k8a(i) = 1.43d-9
        endif
!        if (ttt .le. 1.0d4) then
!          k8a(i) = 1.3d-9
!        else if (ttt .gt. 1.0d4) then
!          k8a(i) = 10.**(-8.78+0.113*dlog10(ttt)
!     .           - 3.475d-2*(dlog10(ttt))**2.)
!        endif
        k9a(i) = 1.85d-23*ttt**1.8
        if(ttt .gt. 6.7e3) &
           k9a(i) =5.81d-16*(ttt/56200)**(-0.6657*dlog10(ttt/56200))
        k10a(i) = 6.0d-10
        if (tev .gt. 0.3) then
          k11a(i) = dexp(-24.24914687731536 &
                  + 3.400824447095291*logtev &
                  - 3.898003964650152*logtev**2 &
                  + 2.045587822403071*logtev**3 &
                  - 0.5416182856220388*logtev**4 &
                  + 0.0841077503763412*logtev**5 &
                  - 0.007879026154483455*logtev**6 &
                  + 0.0004138398421504563*logtev**7 &
                  - 9.36345888928611d-6*logtev**8)
!          k12a(i) = 4.38d-10*dexp(-102000.0/ttt)*ttt**0.35 ! Shapiro & Kang 87
          k12a(i) = 5.6d-11*exp(-102124/ttt)*ttt**0.5
          k13a(i) = 1.0670825d-10*tev**2.012/ &
                    (dexp(4.463/tev)*(1+0.2472*tev)**3.512)
        else
          k11a(i) = tiny
          k12a(i) = tiny
          k13a(i) = tiny 
        endif
        if (tev .gt. 0.04) then
          k14a(i) = dexp(-18.01849334273 &
                  + 2.360852208681*logtev &
                  - 0.2827443061704*logtev**2 &
                  + 0.01623316639567*logtev**3 &
                  - 0.03365012031362999*logtev**4 &
                  + 0.01178329782711*logtev**5 &
                  - 0.001656194699504*logtev**6 &
                  + 0.0001068275202678*logtev**7 &
                  - 2.631285809207d-6*logtev**8)
        else
          k14a(i) =  tiny
        endif
        if (tev .gt. 0.1) then
          k15a(i) = dexp(-20.37260896533324 &
                  + 1.139449335841631*logtev &
                  - 0.1421013521554148*logtev**2 &
                  + 0.00846445538663*logtev**3 &
                  - 0.0014327641212992*logtev**4 &
                  + 0.0002012250284791*logtev**5 &
                  + 0.0000866396324309*logtev**6 &
                  - 0.00002585009680264*logtev**7 &
                  + 2.4555011970392d-6*logtev**8 &
                  - 8.06838246118d-8*logtev**9)
        else
          k15a(i) = 2.56d-9*tev**1.78186
        endif
!          k15a(i) = 5.3d-20*(ttt**2.17)*dexp(-8750./ttt) ! Shapiro & Kang 1987
        k16a(i) = 6.5d-9/dsqrt(tev)
!        k16a(i) = 4.0d-6*dsqrt(ttt)      ! Shapiro & Kang 1987
        k17a(i) = 1.0d-8*ttt**(-0.4)     ! Shapiro & Kang 1987
        if(ttt.gt.1.0e4) &
           k17a(i)=4.0d-4*ttt**(-1.4)*dexp(-15100.0/ttt)
        k18a(i) = 5.56396d-8/tev**0.6035
!        k18a(i) = 1.d-8 
!        if (ttt .gt. 617.)
!     .    k18a(i) = 1.32d-6 * ttt**(-0.76) 
        k19a(i) = 5.d-7*sqrt(100./ttt) 
!        k19a(i) = 4.64d-8/dsqrt(tev)
!
!       ------ 3-body H2 rate ----
!       The first bit is my fit to A.E. Orel 1987, J.Chem.Phys., 87, 
!       314. I then match it to the T^-1 of Palla etal (1983) 
!       Which is then 4 times smaller than the Palla rate ! Thats
!       molecule uncertainties ! :-)
!
        if (ttt .le. 300.0) then
           k22a(i) = 1.3d-32 * (ttt/300.0)**(-0.38) 
        else
           k22a(i) = 1.3d-32 * (ttt/300.0)**(-1.0)
        endif
      endif  ! nspec > 6
      if (nspec .gt. 9) then
!
        k50a(i) = 1.0d-9 *exp(-4.1d1/ttt)  
        k51a(i) = 1.0d-9                   
        k52a(i) = 2.1d-9                   
        k53a(i) = 1.0d-9 *exp(-4.57d2/ttt) 
        k54a(i) = 7.5d-11*exp(-3.82d3/ttt) 
        k55a(i) = 7.5d-11*exp(-4.24d3/ttt) 
        k56a(i) = 1.5d-9 *(ttt/300.0)**(-0.1) 
!
      endif
      enddo
!
      do i = 1, nratec
        logttt = logtmin + dfloat(i-1)*dlogtem
        ttt    = dexp(logttt)
        t_half = dsqrt(ttt)
! Dalgarno & McCray (1972) ISM cooling curves 
      if ( iDMco .eq. 1 .and. ttt .gt. t_cutoff) then
         t_mhalf = ttt**(-0.5)
         t_phalf = ttt**( 0.5)
         if (ttt .lt. 3000. .or. ttt .gt. 11500.) then
            Le_H   = tiny * 1.0d10
         else if (ttt .ge. 3000. .and. ttt .le. 11500.) then
            index   = min0(nHpts-1,max0(0,idint((ttt-3000.)/500.)+1))
            t_lower = 3000. + dfloat(index-1) * 500.
            t_upper = 3000. + dfloat(index  ) * 500.
            Le_H    = dexp(Hcool(index) + (logttt - dlog(t_lower)) * &
                      (Hcool(index+1)-Hcool(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
         endif
         if (ttt .lt. 10.) then
            LH_CI    = tiny * 1.0d10
            LH_O     = tiny * 1.0d10
            LH_C     = tiny * 1.0d10
         else if (ttt .ge. 10. .and. ttt .lt. 100.) then
            index   = min0(9,max0(0,idint((ttt-10.)/10.)+1))
            t_lower = 10. + dfloat(index-1) * 10.
            t_upper = 10. + dfloat(index  ) * 10.
            LH_CI   = dexp(cool1(index) + (logttt - dlog(t_lower)) * &
                      (cool1(index+1)-cool1(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
            LH_O    = dexp(cool2(index) + (logttt - dlog(t_lower)) * &
                      (cool2(index+1)-cool2(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
            LH_C    = dexp(cool3(index) + (logttt - dlog(t_lower)) * &
                      (cool3(index+1)-cool3(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
         else if (ttt .ge. 100.) then
            index   = min0(27,max0(0,idint((ttt-100.)/50.)+10))
            t_lower = 100. + dfloat(index-10) * 50.
            t_upper = 100. + dfloat(index- 9) * 50.
            LH_CI   = dexp(cool1(index) + (logttt - dlog(t_lower)) * &
                      (cool1(index+1)-cool1(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
            LH_O    = dexp(cool2(index) + (logttt - dlog(t_lower)) * &
                      (cool2(index+1)-cool2(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
            LH_C    = dexp(cool3(index) + (logttt - dlog(t_lower)) * &
                      (cool3(index+1)-cool3(index  ))/(dlog(t_upper) &
                      -dlog(t_lower)))
         endif
         Le_CI    = 7.9d-20 * t_mhalf * dexp(-92. /ttt)
         if (Le_CI .lt. tiny .or. Le_CI .ge. 1.)  Le_CI  = tiny * 1.0d10
         Le_SiI   = 1.9d-18 * t_mhalf * dexp(-413./ttt)
         if (Le_SiI .lt. tiny .or. Le_SiI .ge. 1.)Le_SiI = tiny * 1.0d10
         Le_FeI   = 1.1d-18 * t_mhalf *(dexp(-554./ttt) + &
                                  1.3 * dexp(-961./ttt))
         if (Le_FeI .lt. tiny .or. Le_FeI .ge. 1.)Le_FeI = tiny * 1.0d10
         Le_O     = 1.74d-24* t_phalf * &
                    ((1.-7.6* t_mhalf)* dexp(-228./ttt) + &
             0.38 *  (1.-7.7* t_mhalf)* dexp(-326./ttt))
         if (Le_O .lt. tiny .or. Le_O .ge. 1.) Le_O = tiny * 1.0d10
         Lem_O    = 9.4d-23 * t_phalf * dexp(-22700./ttt)
         if (Lem_O .lt. tiny .or. Lem_O .ge. 1.) Lem_O = tiny * 1.0d10
         Lem_1FeI = 4.8d-18 * t_mhalf * dexp(-2694. /ttt)
         if (Lem_1FeI .lt. tiny .or. Lem_1FeI .ge. 1.) Lem_1FeI = tiny * &
                                                                  1.0d10
         Lem_2FeI = 7.8d-18 * t_mhalf * dexp(-3496. /ttt)
         if (Lem_2FeI .lt. tiny .or. Lem_2FeI .ge. 1.) Lem_2FeI = tiny * &
                                                                  1.0d10
         Lem_OI   = 1.5d-17 * t_mhalf * dexp(-38600./ttt)
         if (Lem_OI .lt. tiny .or. Lem_OI .ge. 1.)Lem_OI = tiny * 1.0d10
         Lem_CI   = 3.0d-17 * t_mhalf * dexp(-61900./ttt)
         if (Lem_CI .lt. tiny .or. Lem_CI .ge. 1.)Lem_CI = tiny * 1.0d10
         Lem_SI   = 8.4d-18 * t_mhalf * dexp(-21400./ttt)
         if (Lem_SI .lt. tiny .or. Lem_SI .ge. 1.)Lem_SI = tiny * 1.0d10
         Lem_SiI  = 3.0d-17 * t_mhalf * dexp(-63600./ttt)
         if (Lem_SiI .lt. tiny .or. Lem_SiI .ge. 1.) Lem_SiI = tiny * &
                                                               1.0d10
         Lem_N    = 8.2d-22 * t_phalf * dexp(-27700./ttt) * &
                                       (1.-2.7d-9*ttt**2)
         if (Lem_N .lt. tiny .or. Lem_N .ge. 1.) Lem_N = tiny * 1.0d10
         LH_SiI   = 7.4d-23 *  dexp(-413./ttt)
         if (LH_SiI .lt. tiny .or. LH_SiI .ge. 1.)LH_SiI = tiny * 1.0d10
         LH_FeI   = 1.1d-22 * (dexp(-554./ttt) + &
                        1.4 *  dexp(-961./ttt))
         if (LH_FeI .lt. tiny .or. LH_FeI .ge. 1.)LH_FeI = tiny * 1.0d10
         if (ttt .le. 11500) then
         DMa(i)  = C * (ifrac * (Le_CI  + Lem_CI )           + LH_CI ) + &
                   O * (ifrac * (Le_O   + Lem_O    + Lem_OI) + LH_O  ) + &
                   N * (ifrac *  Lem_N                               ) + &
                   Si* (ifrac * (Le_SiI + Lem_SiI)           + LH_SiI) + &
                   Fe* (ifrac * (Le_FeI + Lem_1FeI + Lem_2FeI        ) &
                                                             + LH_FeI) + &
                   S * (ifrac *  Lem_SI                              ) + &
                        ifrac *  Le_H
         else if (ttt .gt. 11500. .and. ttt .lt. 11775.) then
            DMa(i)   = dexp(Hcool(18) + (logttt - dlog(1.1500d4)) * &
                           (cool4(1)  - Hcool(18))/(dlog(1.1775d4)- &
                           dlog(1.1500d4)))
         else if (ttt .ge. 11775) then
            index    = min0(nDMpts-1,max0(0,idint((logttt-lgtmin2) &
                                                      /dlogt2)+1))
            logt_low = lgtmin2 + dfloat(index-1) * dlogt2
            DMa(i)   = dexp(cool4(index) + (logttt - logt_low)* &
                           (cool4(index+1)-cool4(index  ))/dlogt2)
         endif
      else
         DMa(i)    = tiny
      endif
!         write(50,*) ttt, DMa(i), k1a(i), k2a(i)
! Collisional excitations (Black 1981)
      if (iceco .eq. 1 .and. ttt .gt. t_cutoff) then
       if (nspec .gt. 1) then
        ceHIa(i)    = 7.5d-19*dexp(-dmin1(dlog(huge),118348/ttt)) &
                    /(1.+dsqrt(ttt/1.e5))
       endif
       if (nspec .gt. 3) then
        ceHeIa(i)   = 9.1d-27*dexp(-dmin1(dlog(huge),13179/ttt)) &
                    *ttt**(-0.1687)/(1.+dsqrt(ttt/1.e5))
        ceHeIIa(i)  = 5.54d-17*dexp(-dmin1(dlog(huge),473638/ttt)) &
                    *ttt**(-0.397)/(1.+dsqrt(ttt/1.e5))
       endif
      else
       if (nspec .gt. 1) then
        ceHIa(i)    = tiny
       endif
       if (nspec .gt. 3) then
        ceHeIa(i)   = tiny
        ceHeIIa(i)  = tiny
       endif
      endif
! Collisional ionizations (Tom's polynomial fits)
      if (icico .eq. 1 .and. ttt .gt. t_cutoff) then
       if (nspec .gt. 1) then       
        ciHIa(i)    = 2.18d-11*k1a(i)
       endif 
       if (nspec .gt. 3) then
        ciHeIa(i)   = 3.94d-11*k3a(i)
        ciHeISa(i)  = 5.01d-27*(ttt)**(-0.1687)/(1.+sqrt(ttt/1.0d5)) &
                    * dexp(-dmin1(dlog(huge),55338/ttt)) 
        ciHeIIa(i)  = 8.72d-11*k5a(i)
       endif 
      else
       if (nspec .gt. 1) then
        ciHIa(i)    = tiny
       endif 
       if (nspec .gt. 3) then
        ciHeIa(i)   = tiny
        ciHeIIa(i)  = tiny
        ciHeISa(i)  = tiny
       endif 
      endif
! Recombinations (Cen 1992)
      if (ireco.eq.1 .and. ttt.gt.t_cutoff) then
      if (nspec .gt. 1) then
      if (iOTS  .eq. 0) then
        reHIIa(i)   = 8.70d-27*dsqrt(ttt)*(ttt/1000.0)**(-0.2) &
                    / (1.0 + (ttt/1.e6)**(0.7))
      endif
      if (iOTS .eq. 1) then
!       reHIIa(i) = 1.5  * boltz * ttt * k2a(i)  ! FTB90 benchmark
       log10ttt  = dlog10(ttt)
       if (ttt .ge. 10. .and. ttt .le. 1.0d07) then
        index     = min0(nHIIpts-1,max0(0,idint((log10ttt-1.0)/0.2)+1))
        logt_low  = 1. + dfloat(index-1) * 0.2
        reHIIa(i) = 10.0**(HIIBc(index  ) + (log10ttt - logt_low) * &
                          (HIIBc(index+1) - HIIBc(index))/0.2)    * &
                          boltz * t_half
       endif
      endif
      endif  ! nspec > 1       
      if (nspec .gt. 3) then
      if (iOTS  .eq. 0) then
        reHeII1a(i) = 1.55d-26*ttt**0.3647
        reHeII2a(i) = 1.24d-13*ttt**(-1.5) &
                    *dexp(-dmin1(dlog(huge),470000/ttt)) &
                    *(1.+0.3*dexp(-dmin1(dlog(huge),94000/ttt)))
        reHeIIIa(i) = 3.48d-26*dsqrt(ttt)*(ttt/1000.0)**(-0.2) &
                    / (1.0 + (ttt/1.e6)**(0.7))
      endif
      if (iOTS .eq. 1) then
       log10ttt    = dlog10(ttt)
       if (ttt .ge. 10. .and. ttt .le. 2.5d04) then
        index       = min0(nHeIIpts-1,max0(0, &
                      idint((log10ttt-1.0)/0.2)+1))
        logt_low    = 1. + dfloat(index-1) * 0.2
        reHeII1a(i) = 10.0**(HeIIBc(index  ) + (log10ttt - logt_low) * &
                            (HeIIBc(index+1) -  HeIIBc(index))/0.2)  * &
                             boltz * t_half
       endif
        reHeII2a(i) = 1.24d-13*ttt**(-1.5) &
                    * dexp(-dmin1(dlog(huge),470000/ttt)) &
                    * (1.+0.3*dexp(-dmin1(dlog(huge),94000/ttt)))
       if (ttt .ge. 40. .and. ttt .le. 1.0d07) then
        index       = min0(nHeIIIpts-1,max0(0, &
                      idint((log10ttt-1.6)/0.2)+1))
        logt_low    = 1.6 + dfloat(index-1) * 0.2
        reHeIIIa(i) = 10.0**(HeIIIBc(index  ) + (log10ttt - logt_low) * &
                            (HeIIIBc(index+1) -  HeIIIBc(index))/0.2) * &
                             boltz * t_half
       endif  
      endif   ! iOTS = 1
      endif   ! nspec > 3
      else    ! ireco/t_cutoff
       if (nspec .gt. 1) then
        reHIIa (i)  = tiny
       endif 
       if (nspec .gt. 3) then
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
       endif 
      endif   ! ireco/t_cutoff
! Compton cooling (Peebles 1971)
      if (nspec .gt. 1) then
       if ( icmpt .eq. 1 .and. ttt .ge. 2.73*(1+rdshft) &
            .and. ttt .gt. t_cutoff) then
         compt(i) = 5.65D-36*(1+rdshft)**4*(ttt-2.73*(1+rdshft))
       endif
! Bremsstrahlung (Black 1981)(Spitzer & Hart 1979)
       if ( ibrco .eq. 1 .and. ttt .ge. 1.0e04) then
         brema(i)    = 1.43d-27*dsqrt(ttt) &
                     *(1.1+0.34*dexp(-(5.5-dlog10(ttt))**2/3.0))
       endif
      endif !nspec > 1
! Bremsstrahlung (Shapiro & Kang 1987)(Spitzer 1978)
!        if(ttt .lt. 1.0*3.2e5) gaunt = 0.79464 + 0.1243*log10(ttt/1.)
!        if(ttt .ge. 1.0*3.2e5) gaunt = 2.13164 - 0.1240*log10(ttt/1.)
!        brem1a(i) = 1.426d-27*sqrt(ttt)*gaunt
!        if(ttt .lt. 4.0*3.2e5) gaunt = 0.79464 + 0.1243*log10(ttt/4.)
!        if(ttt .ge. 4.0*3.2e5) gaunt = 2.13164 - 0.1240*log10(ttt/4.)
!        brem2a(i) = 1.426d-27*sqrt(ttt)*gaunt
      if (nspec .gt. 6) then
      if (ttt .gt. t_cutoff) then
!  Molecular hydrogen cooling
!
!  Lepp & Shull rates
!
       if (ih2co .eq. 1) then
        xx = dlog10(ttt/1.0e4)
        vibha  (i) = 1.1d-18*dexp(-dmin1(dlog(huge),6744.0/ttt))
        if( ttt .gt. 1635) then
          dum = 1.0d-12*dsqrt(ttt)*dexp(-1000.0/ttt)
        else
          dum = 1.4d-13*dexp((ttt/125.) - (ttt/577.)**2)
        endif
        hyd01ka(i) = dum*dexp(-dmin1(dlog(huge),8.152d-13/ &
                             (1.38d-16*ttt)))
        dum = 8.152d-13*(4.2/(1.38d-16*(ttt+1190))+1./(1.38d-16*ttt))
        h2k01a(i) = 1.45d-12*dsqrt(ttt)*dexp(-dmin1(dlog(huge),dum))
        if(ttt .gt. 4031) then
          rotla(i) = 1.38d-22*dexp(-9243.0/ttt)
        else
          rotla(i) = 10.0**(-22.9 - 0.553*xx - 1.148*xx*xx)
        endif
        if(ttt .gt. 1087) then
          rotha(i) = 3.9d-19*dexp(-6118.0/ttt)
        else
          rotha(i) = 10.0**(-19.24 + 0.474*xx - 1.247*xx*xx)
        endif
       endif
!
!  Galli and Palla (1999) rates as fit by Tom Abel
!
       if (ih2co .eq. 2) then
        tm = dmax1(ttt, 13.0d0)   ! no cooling below 13 Kelvin ...
        tm = dmin1(tm, 1.d5)      ! fixes numerics
        lt = log10(tm)
!     low density limit from Galli and Palla
        gpldla(i)  = 10.**(-103.0+97.59*lt-48.05*lt**2+10.80*lt*lt*lt &
             -0.9032*lt*lt*lt*lt) 
!     high density limit from HM79      
        t3 = tm/1000.
        HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)* &
                  exp(-(0.13/t3)**3)+3.e-24*exp(-0.51/t3))
        HDLV = (7.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3))
        gphdla(i)  = HDLR + HDLV 
       endif
!
!  Glover (2008) H2 cooling rates, including al relevant collision 
!  channels
!
       if (ih2co .eq. 3) then
!
!     e part 3) - Low density rates from Glover & Abel (2008)
!
!  Excitation by HI
!
         tm  = dmax1 (ttt, 10.0d0)
         tm  = dmin1 (tm , 1.d4  )
         lt3 = dlog10(tm / 1.d3  )  
         if      (tm .lt. 1e2) then
           gaHIa(i) = 10**(-16.818342d0 &
                    + 37.383713d0 * lt3 &
                    + 58.145166d0 * lt3**2 &
                    + 48.656103d0 * lt3**3 &
                    + 20.159831d0 * lt3**4 &
                    + 3.8479610d0 * lt3**5)
         else if (tm .lt. 1e3) then
           gaHIa(i) = 10**(-24.311209d0 &
                    + 3.5692468d0 * lt3 &
                    - 11.332860d0 * lt3**2 &
                    - 27.850082d0 * lt3**3 &
                    - 21.328264d0 * lt3**4 &
                    - 4.2519023d0 * lt3**5) 
         else
           gaHIa(i) = 10**(-24.311209d0 &
                    + 4.6450521d0 * lt3 &
                    - 3.7209846d0 * lt3**2 &
                    + 5.9369081d0 * lt3**3 &
                    - 5.5108047d0 * lt3**4 &
                    + 1.5538288d0 * lt3**5)
         endif
!
! Excitation by H2
!
         gaH2a(i) = 10**(-23.962112d0 &
                  + 2.09433740d0  * lt3 &
                  - 0.77151436d0  * lt3**2 &
                  + 0.43693353d0  * lt3**3 &
                  - 0.14913216d0  * lt3**4 &
                  - 0.033638326d0 * lt3**5)
!
! Excitation by He
!
         gaHea(i) = 10**(-23.689237d0 &
                  + 2.1892372d0  * lt3 &
                  - 0.81520438d0 * lt3**2 &
                  + 0.29036281d0 * lt3**3 &
                  - 0.16596184d0 * lt3**4 &
                  + 0.19191375d0 * lt3**5)
!
! Excitation by H+
!
         gaHpa(i) = 10**(-21.716699d0 &
                  + 1.3865783d0   * lt3 &
                  - 0.37915285d0  * lt3**2 &
                  + 0.11453688d0  * lt3**3 &
                  - 0.23214154d0  * lt3**4 &
                  + 0.058538864d0 * lt3**5) 
!
! Excitation by electrons
!
         if (tm .lt. 200) then
           gaela(i) = 10**(-34.286155d0 &
                    - 48.537163d0  * lt3 &
                    - 77.121176d0  * lt3**2 &
                    - 51.352459d0  * lt3**3 &
                    - 15.169160d0  * lt3**4 &
                    - 0.98120322d0 * lt3**5) 
         else
           gaela(i) = 10**(-22.190316 &
                    + 1.5728955  * lt3 &
                    - 0.21335100 * lt3**2 &
                    + 0.96149759 * lt3**3 &
                    - 0.91023195 * lt3**4 &
                    + 0.13749749 * lt3**5) 
         endif       
      endif ! iH2co = 3
      endif ! ttt/t_cutoff
!
      endif ! nspec > 6
      if (nspec .gt. 9) then 
      if (iHDco.eq.1 .and. ttt.gt.t_cutoff) then
!
!        HD Cooling Function (ergs cm3 /s)
!
! Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
!
        hdltea(i) = -35.6998d0 + 15.35716d0*dlog10(ttt) - &
                    5.58513d0   * (dlog10(ttt))**2 + &
                    0.8561149d0 * (dlog10(ttt))**3 - &
                    1.75538d-2  * (dlog10(ttt))**4
        hdltea(i) = (10.0**dmin1(hdltea(i), 15.d0)) 
!      hdlte=10**lhdlte/den
!
! Galli and Palla 1998 low density limit (erg cm3 /s)
! uses HD-He collisional data, so reduced by 1.27
!
        hdlowa(i) = ((3.0d0 * (4.4d-12 + 3.6d-13*ttt**0.77) * &
                    dexp(-128.d0/ttt) * 128.d0 + &
                    (5.0d0/3.0d0) * (4.1d-12+2.1e-13*ttt**0.92) * &
                    dexp(-255.d0/ttt) * 255.d0) * boltz/1.27d0 ) 
!      hdcool=hdlte/(1.0d0+hdlte/hdlow)
!
      endif
      endif
      enddo ! nratec
      if (iRT .gt. 0) then
! Generate array of photon energies
      if (ifreq .eq. 1) then
         nu(1) = ephot
      endif
      if (ifreq .eq. 2) then
        hnumin = 13.602
        dlognu = (dlog10(hnumax*everg/hplanck) &
               -  dlog10(hnumin*everg/hplanck)) / dfloat(nnu2)
! Reserve nnu1 frequency bins for energies below the H ionization edge
! and nnu2 bins for hnu > 13.602 eV
        if (nspec .gt. 6) then
! H- photodetachment and H2+ dissociation energies (0.755 eV to 13.6 eV)
          delnu = (e24 - e27)/nnu1
          do l = 1, nnu1
            nu(l) = e27 + (l-1) * delnu
          enddo
        endif ! nspec > 6
! Ionizing photon energies
        do l = 1, nnu2
          nu(l+nnu1) = 10**(15.518 + (l-1) * dlognu)
          nu(l+nnu1) = hplanck * nu(l+nnu1) / everg
        enddo
	if (ispct .eq. 2) then
          kt     = t_star/tevk
	  L_star = lsol * 10.**(L_star)
	endif
      endif ! ifreq = 2
! Compute the photon interaction cross sections
!
      do l = 1, nnu1+nnu2
       if (nspec .gt. 1) then
!     H photoionization cross sections
        if ( nu(l) .gt. e24 ) then
          dum = dsqrt(nu(l)/e24-1)
          sigma24(l) = 6.3d-18*(e24/nu(l))**4*dexp &
                      (4.0-4.0*datan(dum)/dum) / (1-dexp(-2.0*pi/dum))
        else
          sigma24(l) = tiny
        endif
       endif 
       if (nspec .gt. 3) then
!     HeII hydrogenic photoionization cross sections
        if ( nu(l) .gt. e25 ) then
          dum = dsqrt(nu(l)/e25-1)
          sigma25(l) = 1.58d-18*(e25/nu(l))**4*dexp &
                      (4.0-4.0*datan(dum)/dum) / (1-dexp(-2.0*pi/dum))
        else
          sigma25(l) = tiny
        endif
!     HeI photoionization cross sections
        if ( nu(l) .gt. e26 ) then
          sigma26(l) = 7.42d-18*(1.66*(nu(l)/e26)**(-2.05) &
                               - 0.66*(nu(l)/e26)**(-3.05))
        else
          sigma26(l) = tiny
        endif
       endif 
       if (nspec .gt. 6) then
!     H-  photo detachment cross sections
        if ( nu(l) .gt. e27 .and. iHM .eq. 1) then
          sigma27(l) = 2.11d-16*(nu(l)-e27)**1.5/nu(l)**3   
        else
          sigma27(l) = tiny
        endif
!     H_2+ + gamma ---> H + H+  cross sections  (Shapiro & Kang 1987)
        if ( nu(l).gt.e28 .and. nu(l).le.11.27 ) then
          sigma28(l) = 10**(-40.97+6.03*nu(l)-0.504*nu(l)**2 &
                            +1.387d-2*nu(l)**3)
        elseif ( nu(l).gt.11.27 .and. nu(l).lt.21.0 ) then
          sigma28(l) = 10**(-30.26+2.79*nu(l)-0.184*nu(l)**2 &
                            +3.535d-3*nu(l)**3)
        else
          sigma28(l) = tiny
        endif
!     H_2 + gamma ---> H_2+ + e- cross sections (Shapiro & Kang 1987)
        if ( nu(l).gt.e29 .and. nu(l).le.16.5 ) then
          sigma29(l) = 6.2d-18*nu(l) - 9.4d-17
        elseif ( nu(l).gt.16.5 .and. nu(l).le.17.7 ) then  
          sigma29(l) = 1.4d-18*nu(l) - 1.48d-17
        elseif ( nu(l).gt.17.7 ) then
          sigma29(l) = 2.5d-14*nu(l)**(-2.71)
        else
          sigma29(l) = tiny
        endif
!     H_2+ + gamma ---> 2 H+ + e-  cross sections (Shapiro & Kang 1987)
        if ( nu(l).ge.e30 .and. nu(l).lt.70.0 ) then
          sigma30(l) = 10**(-16.926-4.528d-2*nu(l)+2.238d-4*nu(l)**2 &
                            +4.245d-7*nu(l)**3)
        else
          sigma30(l) = tiny
        endif
        sigma31(l) = 0.0
       endif
       if (nspec .gt. 9) then
!     D photoionization cross sections
        if ( nu(l) .gt. e24 ) then
          sigma32(l) = sigma24(l)
        else
          sigma32(l) = tiny
        endif
       endif
      enddo
! dust extinction (Pei, 1992, ApJ, 395, 130)
      if (idust .eq. 1) then
        allocate(xxi(nnu1+nnu2))
        gtd  = 7.8d-22 ! dust-to-gas ratio for 
                       ! SMC taken from Table 2
        hc = hplanck*clight*1.0d4/everg
! fit parameters are from Table 4 for SMC
        lmbda(1) = 0.042
        lmbda(2) = 0.08
        lmbda(3) = 0.22
        lmbda(4) = 9.7
        lmbda(5) = 18.
        lmbda(6) = 25.
        a(1)     = 185.
        a(2)     = 27.
        a(3)     = 0.005
        a(4)     = 0.010
        a(5)     = 0.012
        a(6)     = 0.030
        b(1)     = 90.
        b(2)     = 5.5
        b(3)     = -1.95
        b(4)     = -1.95
        b(5)     = -1.80
        b(6)     = 0.
        c1(1)    = 2.
        c1(2)    = 4.
        c1(3)    = 2.
        c1(4)    = 2.
        c1(5)    = 2.
        c1(6)    = 2.
! fit parameters are from Table 4 for MW
        lmbda(1) = 0.047
        lmbda(2) = 0.08
        lmbda(3) = 0.22
        lmbda(4) = 9.7
        lmbda(5) = 18.
        lmbda(6) = 25.
        a(1)     = 165.
        a(2)     = 14.
        a(3)     = 0.045
        a(4)     = 0.002
        a(5)     = 0.002
        a(6)     = 0.012
        b(1)     = 90.
        b(2)     = 4.0
        b(3)     = -1.95
        b(4)     = -1.95
        b(5)     = -1.80
        b(6)     = 0.
        c1(1)    = 2.
        c1(2)    = 6.5
        c1(3)    = 2.
        c1(4)    = 2.
        c1(5)    = 2.
        c1(6)    = 2.
! construct xi (=tau_ext/tau_B --eqs 8 & 20)
        do l = 1, nnu1+nnu2
          xxi(l) = 0.
          do i = 1,6
            xxi(l) = xxi(l) + a(i) / ((hc/(nu(l)*lmbda(i)))**c1(i)+ &
                                      (hc*lmbda(i) / nu(l))**c1(i)+ &
                                       b(i))
          enddo
          xxi(l) = xxi(l) * gtd
        enddo
      endif ! idust
      endif ! iRT
! 77   format(41(1pe11.3e3,1x))
!      open(unit=10, file='ctable1.dat',status='unknown')
!      do i=1,nratec
!        write(10,77)k1a (i),k2a (i),k3a (i),k4a (i),k5a (i),k6a (i),
!     .              k7a (i),k8a (i),k9a (i),k10a(i),k11a(i),k12a(i),
!     .              k13a(i),k14a(i),k15a(i),k16a(i),k17a(i),k18a(i),
!     .              k19a(i),k22a(i),
!     .
!     .              ceHIa   (i),ciHIa   (i),ceHeIa (i),ciHeIa  (i),
!     .              ciHeISa (i),ciHeIIa (i),reHIIa (i),reHeII1a(i),
!     .              reHeII2a(i),reHeIIIa(i),hyd01ka(i),h2k01a  (i),
!     .              vibha   (i),rotha   (i),rotla  (i),gpldla  (i),
!     .              gphdla  (i),brema   (i),compt  (i)
!      enddo
!      close(unit=10)
!      stop
      return
      end subroutine ctable
