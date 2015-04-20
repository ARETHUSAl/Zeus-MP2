#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E     //////////
c    //////////            C O O L C H E M 3 D           \\\\\\\\\\
c
c=======================================================================
c
       subroutine coolchem3D 
    
c
c  PURPOSE:  Computes cooling by subcycle and then utilizes the
c            backward difference formula (BDF) to solve the chemistry
c            rate equations, which describe the updates to the hydrogen 
c            and helium densities due to photoionization, collision, 
c            and recombination.  Two options are currently implemented 
c            for advancing the cooling and rate equations.  The first 
c            option utilizes the cooling timestep to govern both the
c            cooling and chemistry updates.  The second option subcycles
c            the rate solve and heating/cooling separately according
c            to their respective timescales (e/edot for heating, 
c            de/dedot for the rate equations where de is the electron
c            density n_e * m_H).  Mass and charge density conservation
c            are enforced every timestep by calling mconserv.  Invoking
c            the RATE_COOL macro in zeusmp.def activates option 1 to
c            advance both the BDF rate solve and cooling update by the
c            cooling timestep.
c              
c            Note on HERCULES unit conversion patterns:  the rate solve
c            and heating code in this subroutine was adapted from the
c            HERCULES cosmological code used in the 1997 Anninos et al. 
c            simulations (New Astronomy 2 (1997) 209-224).  HERCULES
c            hydrodynamics was done in cosmological units with megapc
c            and the Hubble time being the characteristic length and 
c            time scales.  Densities were also scaled by cosmological
c            units so e and d in are in cosmological units rather than
c            cgs in that code.  The chemistry and cooling, however, 
c            were done in cgs units and then converted back into cosmo-
c            logical units to update d and e.  Those unit conversions 
c            appear in the code below.  They can be negated by being set 
c            equal to 1 and this is done in Zeus-MP because all its 
c            calculations are carried out in the cgs units of the rate 
c            solve and cooling.
c
c 
c                cubelmpc = cosmological cube length in megaparsecs
c                cmpmpc   = cm per megaparsec = 3.086e+24
c                h0       = Hubble time/sec   = 3.240e-18
c                omega0   = HERCULES user input parameter
c                h        = Hubble constant at z=0 in units of 100 
c                           km/sec-Mpc
c                pii      = 3.14159265359
c                gee      = 6.67e-08
c
c                xbase1   = cubelmpc*cmpmpc/(# of x-gridpoints)
c                tbase1   = 1./h0
c                dbase1   = 3.*(h0*h)**2*omega0/(8.*pii*gee)
c                dom      = dbase1/mh
c                xb3      = xbase1**3
c                       
c            
c            multiplying a cosmol. length by xbase1 converts it to cm
c                                  time      tbase1                sec
c                                  density   dom               gm/cm^3
c
c            dividing a cm length by xbase1 converts it to cosmol. units
c                      sec time      tbase1                  "       "
c                  gm/cm^3 density   dom                     "       "
c
c            Hence, in ctable the rate constants k (cm^3/sec) are  
c            converted to cosm units by dividing by kunit (xbase1**3/
c            tbase1) (and are later converted back in the BDF eqns to  
c            cgs by multiplying by kunit).  Cosmological densities in 
c            both rate and cooling eqns are converted to cgs by multi-
c            plying by dom.  The heating/cooling rates themselves are
c            then converted back to cosm units by dividing by edotbase
c            to update the energies and the densities are switched 
c            back by dividing by dom.  edotbase is dbase1*xbase1^2/
c            tbase1^2 and is therefore an energy/vol conversion term.
c            t and dt in Hubble times are converted to sec by multi-
c            plying them by tbase1 in the rate equations.
c
c            As noted earlier, ZeusMP cgs densities, times and energies 
c            are already compatible with the rate and cooling equations 
c            so we set xbase1, tbase1, and dbase1 to 1.0.  We retain the 
c            conversion constants in the equations only to maintain 
c            consistency with HERCULES.
c    
c  LOCAL VARIABLES:
c
c  BOUNDARY VALUES USED:
c
c  EXTERNALS:  mconsrv
c
c-----------------------------------------------------------------------
c
      use real_prec
      use param
      use cons
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

      integer  :: i,j,k,n,itmax,imin,jmin,kmin
      integer  :: iswres,iswhdf,iswhst,iswusr

      real(rl) :: qq,vibl,t_photo,t_heat,dtchem,t_chem
     .,           ttmin,cnorm,scoef,acoef,logtem0,dtimin
     .,           logtem9,mhi,tphmin,thtmin,tsubcycle
#ifdef H
     .,           HIp,HIIp,dep
#endif /* H */
#ifdef He
     .,           HeIp,HeIIp,HeIIIp
#endif /* He */
#ifdef HD
     .,           DIp,DIIp,HDIp,hdlte1,hdlow1
#endif /* HD */

      real(rl), dimension(in,jn) :: t1, t2, tdef, f_rad

      external mconsrv, dataio, mnu_dir_flx

#ifndef UNROLL_I
#define UNROLL_I
#endif
c
c-----------------------------------------------------------------------
c
c

      mhi      = 1./mh
      logtem0  = dlog(temstart)
      logtem9  = dlog(temend)
c      cnorm    = boltz*tevk/(1.22*mh)
      cnorm    = boltz*tevk*mhi
      iswres   = 0
      iswhdf   = 0
      iswhst   = 0
      iswusr   = 1

      tphmin   = huge
      thtmin   = huge
      itmax    = 10000

c      do i = is, ie

c         tphi   (i) = tiny
c         tlambda(i) = tiny

c      enddo
      

c     First, compute the heating/cooling timestep by which we'll advance
c     the hydro.  This requires a sweep over the entire problem domain to
c     find the global min timestep

     
      do k = ks, ke

#ifdef RT
c     Solve the radiative transfer equation for the photoionzation
c     rates k24 thru k26 and the photoheating rates piH*
 
          call mnu_dir_flx(k)
#endif /* RT */

c     Next, compute the heat/cool and chemistry coefficients
          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie
        
c          rad_press(i,j,k) = 0.   

          tgas(i,j,k) = gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / (
#ifdef H
     .        HII(i,j,k)       + HI  (i,j,k)       + de   (i,j,k)
#endif /* H */
#ifdef He
     .      + HeI(i,j,k) * qrt + HeII(i,j,k) * qrt + HeIII(i,j,k) * qrt
#endif /* He */
#ifdef H2
     .      + HM (i,j,k)       + H2II(i,j,k) * haf + H2I  (i,j,k) * haf
#endif /* H2 */
     .                )
            logtem(i,j) = dlog(tgas(i,j,k))
            if (logtem(i,j) .lt. logtem0) logtem(i,j) = logtem0
            if (logtem(i,j) .ge. logtem9) logtem(i,j) = logtem9
            indixe(i,j) = min0(nratec-1,
     .                   max0(1,idint((logtem(i,j)-logtem0)/dlogtem)+1))
            t1    (i,j) = (logtem0 + (indixe(i,j) - 1)*dlogtem)
            t2    (i,j) = (logtem0 + (indixe(i,j)    )*dlogtem)
            tdef  (i,j) = t2(i,j) - t1(i,j)
#ifdef H
      ceHI   (i,j) = ceHIa   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ceHIa   (indixe(i,j)+1)-ceHIa   (indixe(i,j)))/tdef(i,j)

      ciHI   (i,j) = ciHIa   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHIa   (indixe(i,j)+1)-ciHIa   (indixe(i,j)))/tdef(i,j)

      reHII  (i,j) = reHIIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHIIa  (indixe(i,j)+1)-reHIIa  (indixe(i,j)))/tdef(i,j)

      DM     (i,j) = DMa     (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(DMa     (indixe(i,j)+1)-DMa     (indixe(i,j)))/tdef(i,j)
#endif /* H */

#ifdef He
      ceHeI  (i,j) = ceHeIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ceHeIa  (indixe(i,j)+1)-ceHeIa  (indixe(i,j)))/tdef(i,j)

      ceHeII (i,j) = ceHeIIa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ceHeIIa (indixe(i,j)+1)-ceHeIIa (indixe(i,j)))/tdef(i,j)

      ciHeI  (i,j) = ciHeIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeIa  (indixe(i,j)+1)-ciHeIa  (indixe(i,j)))/tdef(i,j)

      ciHeIS (i,j) = ciHeISa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeISa (indixe(i,j)+1)-ciHeISa (indixe(i,j)))/tdef(i,j)

      ciHeII (i,j) = ciHeIIa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeIIa (indixe(i,j)+1)-ciHeIIa (indixe(i,j)))/tdef(i,j)

      reHeII1(i,j) = reHeII1a(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeII1a(indixe(i,j)+1)-reHeII1a(indixe(i,j)))/tdef(i,j)

      reHeII2(i,j) = reHeII2a(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeII2a(indixe(i,j)+1)-reHeII2a(indixe(i,j)))/tdef(i,j)

      reHeIII(i,j) = reHeIIIa(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeIIIa(indixe(i,j)+1)-reHeIIIa(indixe(i,j)))/tdef(i,j)
#endif /* He */

#ifdef H2
      if (iH2co .eq. 1) then
      hyd01k (i,j) = hyd01ka (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hyd01ka (indixe(i,j)+1)-hyd01ka (indixe(i,j)))/tdef(i,j)

      h2k01  (i,j) = h2k01a  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(h2k01a  (indixe(i,j)+1)-h2k01a  (indixe(i,j)))/tdef(i,j)

      vibh   (i,j) = vibha   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(vibha   (indixe(i,j)+1)-vibha   (indixe(i,j)))/tdef(i,j)

      roth   (i,j) = rotha   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(rotha   (indixe(i,j)+1)-rotha   (indixe(i,j)))/tdef(i,j)

      rotl   (i,j) = rotla   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(rotla   (indixe(i,j)+1)-rotla   (indixe(i,j)))/tdef(i,j)

      else if (iH2co .eq. 2) then
      gpldl  (i,j) = gpldla  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(gpldla  (indixe(i,j)+1)-gpldla  (indixe(i,j)))/tdef(i,j)

      gphdl  (i,j) = gphdla  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(gphdla  (indixe(i,j)+1)-gphdla  (indixe(i,j)))/tdef(i,j)
      endif
#endif /* H2 */

#ifdef HD
      hdlte  (i,j) = hdltea  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hdltea  (indixe(i,j)+1)-hdltea  (indixe(i,j)))/tdef(i,j)

      hdlow  (i,j) = hdlowa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hdlowa  (indixe(i,j)+1)-hdlowa  (indixe(i,j)))/tdef(i,j)
#endif /* HD */

      brem   (i,j) = brema   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(brema   (indixe(i,j)+1)-brema   (indixe(i,j)))/tdef(i,j)

      cmpt   (i,j) = compt   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(compt   (indixe(i,j)+1)-compt   (indixe(i,j)))/tdef(i,j)

#ifdef H
          k1 (i,j) = k1a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k1a (indixe(i,j)+1) -  k1a(indixe(i,j)))/tdef(i,j)
          k2 (i,j) = k2a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k2a (indixe(i,j)+1) -  k2a(indixe(i,j)))/tdef(i,j)
#endif /* H */
#ifdef He
          k3 (i,j) = k3a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k3a (indixe(i,j)+1) -  k3a(indixe(i,j)))/tdef(i,j)
          k4 (i,j) = k4a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k4a (indixe(i,j)+1) -  k4a(indixe(i,j)))/tdef(i,j)
          k5 (i,j) = k5a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k5a (indixe(i,j)+1) -  k5a(indixe(i,j)))/tdef(i,j)
          k6 (i,j) = k6a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k6a (indixe(i,j)+1) -  k6a(indixe(i,j)))/tdef(i,j)
#endif /* He */
#ifdef H2
          k7 (i,j) = k7a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k7a (indixe(i,j)+1) -  k7a(indixe(i,j)))/tdef(i,j)
          k8 (i,j) = k8a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k8a (indixe(i,j)+1) -  k8a(indixe(i,j)))/tdef(i,j)
          k9 (i,j) = k9a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k9a (indixe(i,j)+1) -  k9a(indixe(i,j)))/tdef(i,j)
          k10(i,j) = k10a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k10a(indixe(i,j)+1) - k10a(indixe(i,j)))/tdef(i,j)
          k11(i,j) = k11a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k11a(indixe(i,j)+1) - k11a(indixe(i,j)))/tdef(i,j)
          k12(i,j) = k12a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k12a(indixe(i,j)+1) - k12a(indixe(i,j)))/tdef(i,j)
          k14(i,j) = k14a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k14a(indixe(i,j)+1) - k14a(indixe(i,j)))/tdef(i,j)
          k16(i,j) = k16a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k16a(indixe(i,j)+1) - k16a(indixe(i,j)))/tdef(i,j)
          k17(i,j) = k17a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k17a(indixe(i,j)+1) - k17a(indixe(i,j)))/tdef(i,j)
          k18(i,j) = k18a(indixe(i,j))   + (logtem(i,j) - t1  (i,j))
     .             *(k18a(indixe(i,j)+1) - k18a(indixe(i,j)))/tdef(i,j)
          k19(i,j) = k19a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k19a(indixe(i,j)+1) - k19a(indixe(i,j)))/tdef(i,j)
          k22(i,j) = k22a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k22a(indixe(i,j)+1) - k22a(indixe(i,j)))/tdef(i,j)
#endif /* H2 */
#ifdef HD
          k50(i,j) = k50a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k50a(indixe(i,j)+1) - k50a(indixe(i,j)))/tdef(i,j)
          k51(i,j) = k51a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k51a(indixe(i,j)+1) - k51a(indixe(i,j)))/tdef(i,j)
          k52(i,j) = k52a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k52a(indixe(i,j)+1) - k52a(indixe(i,j)))/tdef(i,j)
          k53(i,j) = k53a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k53a(indixe(i,j)+1) - k53a(indixe(i,j)))/tdef(i,j)
          k54(i,j) = k54a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k54a(indixe(i,j)+1) - k54a(indixe(i,j)))/tdef(i,j)
          k55(i,j) = k55a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k55a(indixe(i,j)+1) - k55a(indixe(i,j)))/tdef(i,j)
          k56(i,j) = k56a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k56a(indixe(i,j)+1) - k56a(indixe(i,j)))/tdef(i,j)
#endif /* HD */
          enddo
          enddo

c     Now compute the heating and cooling rates of each zone   

          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie

c     Cooling Functions:

      edot(i,j,k)=(

#ifdef H
     . - ceHI   (i,j)* HI   (i,j,k)*de(i,j,k)             ! coll excit of HI
     . - ciHI   (i,j)* HI   (i,j,k)*de(i,j,k)             ! coll ioniz of HI
     . - reHII  (i,j)* HII  (i,j,k)*de(i,j,k)             ! recombin   of HII
     . - DM     (i,j)* HI   (i,j,k)**2                    ! DM ISM cooling
#endif /* H */
#ifdef He
     . - ceHeI  (i,j)* HeII (i,j,k)*de(i,j,k)**2*mhi*qrt  ! coll excit of HeI
     . - ceHeII (i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! coll excit of HeII
     . - ciHeI  (i,j)* HeI  (i,j,k)*de(i,j,k)*qrt         ! coll ioniz of HeI
     . - ciHeII (i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! coll ioniz of HeII
     . - ciHeIS (i,j)* HeII (i,j,k)*de(i,j,k)**2*mhi*qrt  ! coll ioniz of HeIS
     . - reHeII1(i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! recombin   of HeII
     . - reHeII2(i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! recombin   of HeII
     . - reHeIII(i,j)* HeIII(i,j,k)*de(i,j,k)*qrt         ! recombin   of HeIII
#endif /* He */
     . - brem   (i,j)*(
#ifdef H
     .                 HII  (i,j,k)                       ! elec brems cooling
#endif /* H */
#ifdef He
     .               + HeII (i,j,k)*qrt+HeIII(i,j,k)*qrt 
#endif /* He */
     .              )* de   (i,j,k)                       
     . - cmpt   (i,j)* de   (i,j,k)*mh                    ! compton cooling

     .            )*mhi**2

c      tlambda(i)= tlambda(i) + dabs(edot(i,j,k))          ! cooling timescale

#ifdef H2
      if (iH2co .eq. 1) then
      qq        = 1.2*(HI(i,j,k)*mhi)**0.77 + (H2I(i,j,k)*mhi*haf)**0.77
      vibl      = (HI(i,j,k)*hyd01k(i,j) + H2I(i,j,k)*haf*h2k01(i,j))
     .          * mhi*8.18e-13

      edot(i,j,k) = edot(i,j,k) -                         ! ro-vibr cool by H2
     .         ( H2I(i,j,k)*(vibh(i,j)/(1.+vibh(i,j)/vibl)
     .         + roth(i,j)/(1.+roth(i,j)/(qq*rotl(i,j))))*mhi*haf)

      else if (iH2co .eq. 2) then

       gphdl1      = gphdl(i,j  ) / (HI(i,j,k) * mhi)
       edot(i,j,k) = edot (i,j,k) - H2I(i,j,k) *
     .               gphdl(i,j  ) / (1.0 + gphdl1 / gpldl(i,j))
     .             * haf * mhi
      endif

#endif /* H2 */
#ifdef HD
      if (iHDco .eq. 1) then

       hdlte1      = hdlte(i,j)/(HDI(i,j,k)*haf*mhi)
       hdlow1      = max(hdlow(i,j), tiny)
       edot(i,j,k) = edot(i,j,k) - HDI(i,j,k) *
     .               hdlte1/(1.0 + hdlte1/hdlow1)
     .             * haf * mhi
      endif
#endif /* HD */

#ifdef RT
c     Heating functions:

      edot(i,j,k) =   edot(i,j,k) + dfloat(ipiht)*(     ! photoioniz heating
#ifdef H
     .              + piHI  (i,j) * HI  (i,j,k)         ! pi of HI
#endif /* H */
#ifdef He
     .              + piHeI (i,j) * HeI (i,j,k)*qrt     ! pi of HeI
     .              + piHeII(i,j) * HeII(i,j,k)*qrt     ! pi of HeII
#endif /* He */
     .                                            ) * mhi

c      tphi(i  ) =   tphi(i)     + dfloat(ipiht)*(           ! photoioniz timescale
#ifdef H
c     .              + piHI  (i,j) * HI  (i,j,k)         ! pi of HI
#endif /* H */
#ifdef He
c     .              + piHeI (i,j) * HeI (i,j,k)*qrt     ! pi of HeI
c     .              + piHeII(i,j) * HeII(i,j,k)*qrt     ! pi of HeII
#endif /* He */
c     .                                           ) * mhi
#endif /* RT */

          enddo
          enddo

          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie
            dedot(i,j) = (
#ifdef H
     .              k1 (i,j) * HI   (i,j,k) * de (i,j,k)
#endif /* H  */
#ifdef He
     .            + k3 (i,j) * HeI  (i,j,k) * de (i,j,k)*qrt
     .            + k5 (i,j) * HeII (i,j,k) * de (i,j,k)*qrt
#endif /* He */
#ifdef H2
     .            + k8 (i,j) * HM   (i,j,k) * HI (i,j,k)
     .            + k15(i,j) * HM   (i,j,k) * HI (i,j,k)
     .            + k17(i,j) * HM   (i,j,k) * HII(i,j,k)
     .            + k14(i,j) * HM   (i,j,k) * de (i,j,k)
#endif /* H2 */
     .            ) * mhi * mhi - (
#ifdef H
     .              k2 (i,j) * HII  (i,j,k) * de (i,j,k)
#endif /* H  */
#ifdef He
     .            + k4 (i,j) * HeII (i,j,k) * de (i,j,k)*qrt
     .            + k6 (i,j) * HeIII(i,j,k) * de (i,j,k)*qrt
#endif /* He */
#ifdef H2
     .            + k7 (i,j) * HI   (i,j,k) * de (i,j,k)
     .            + k18(i,j) * H2II (i,j,k) * de (i,j,k)*haf
#endif /* H2 */
     .            ) * mhi * mhi
#ifdef RT
     .            + (
#ifdef H
     .                k24(i,j) * HI  (i,j,k)
#endif /* H  */
#ifdef He
     .              + k25(i,j) * HeII(i,j,k)*qrt
     .              + k26(i,j) * HeI (i,j,k)*qrt
#endif /* He */
     .              ) * mhi
#endif /* RT */

c     Compute the photoionization timescale for this zone

            if      (icycle .eq. 1) then
c              t_photo = dabs(0.1*de(i,j,k)/dedot(i,j)*mhi)
              t_photo = dabs(0.1*(de(i,j,k) + 1.0d-3*d(i,j,k))
     .                          /dedot(i,j)*mhi)
            else if (icycle .eq. 2) then
              t_photo = dabs(0.1*(de(i,j,k) + 1.0d-3*d(i,j,k))
     .                          /dedot(i,j)*mhi)

c               t_photo = 0.1 * dmin1((de(i,j,k) + 1.0d-3*d(i,j,k))/
c     .                               (dedot(i,j)*mhi),
c     .                               (HI(i,j,k) + 1.0d-3*d(i,j,k))/      
c     .                               (dedot(i,j)*mhi))

c              t_photo = dmax1(dabs(1.0d-06*HI(i,j,k)/dedot(i,j)*mhi),
c     .                        dabs(0.1*de(i,j,k)/dedot(i,j)*mhi))
c              t_photo = dmax1(dmin1(0.1*dx1a(i)/clight,
c     .                        dabs(0.1*HI(i,j,k)/dedot(i,j)*mhi)),
c     .                        dabs(0.1*de(i,j,k)/dedot(i,j)*mhi))
            endif

c     Compute the heating/cooling time for this zone

            if      (icycle .eq. 1) then
              t_heat  = 0.1 * dabs(e(i,j,k)/edot(i,j,k)) 
            else if (icycle .eq. 2) then
              if (edot(i,j,k) .ge. 0.) then
                t_heat  = dabs(e(i,j,k)/edot(i,j,k)) 
              else
                t_heat  = 0.1 * dabs(e(i,j,k)/edot(i,j,k)) 
              endif
            endif
     

c     Now compare this zone's photoionization timescale to the shortest
c     encountered thus far on the grid and keep the smaller of the two.
c     The shortest photo timescale found on the grid is the maximum
c     timestep by which we can safely advance the entire grid's chemistry
c     and energy.  Do the same with the heating timescale but apply it to
c     instead advance the hydro.  We therefore subcycle the chemistry and
c     rad transfer over photo timesteps until we have covered a heating
c     timestep, at which point it is time to update the hydro.

c            if (t_heat .lt. thtmin) then
c               imin = i
c               jmin = j
c               kmin = k
c            endif

            tphmin = dmin1(t_photo,tphmin)
            thtmin = dmin1(t_heat ,thtmin)

            
          enddo
          enddo

      enddo

c     Adopt the smaller of the nudt hydro dt and dtchem as the hydro timestep

c      dtchem = dmin1(tphmin,thtmin)
      dt     = dmin1(dt    ,thtmin)
c      if (nhy.lt.10000) then 
c        dt = dmin1(dt    ,dmin1(thtmin,0.1*dx1a(is)/clight))
c      else
c        dt     = dmin1(dt,thtmin)
c      endif

c      do i = is, ie

c         tphi   (i) = e (i,js,ks) /(secyr * tphi   (i))
c         tlambda(i) = e (i,js,ks) /(secyr * tlambda(i))
c         trec   (i) = de(i,js,ks) * mhi/
c     .               (k2(i,js   ) * HII(i,js,ks) * 
c     .                de(i,js,ks) * secyr        * mhi**2)
c
c      enddo
      
c
c  Compute min process timescales
c

c      tmin_hydro  = huge
c      tmin_pdv    = huge
c      tmin_phi    = huge
c      tmin_rec    = huge
c      tmin_lambda = huge

c      rmin_hydro  = 0.
c      rmin_pdv    = 0.
c      rmin_phi    = 0.
c      rmin_rec    = 0.
c      rmin_lambda = 0.

c      do i=is,ie
c         if (thydro (i) .lt. tmin_hydro ) then
c            tmin_hydro  = thydro (i)
c            rmin_hydro  = x1b    (i)/cmpc
c         endif
c         if (tpdv   (i) .lt. tmin_pdv   ) then
c            tmin_pdv    = tpdv   (i)
c            rmin_pdv    = x1b    (i)/cmpc
c         endif
c         if (tphi   (i) .lt. tmin_phi   ) then
c            tmin_phi    = tphi   (i)
c            rmin_phi    = x1b    (i)/cmpc
c         endif
c         if (trec   (i) .lt. tmin_rec   ) then
c            tmin_rec    = trec   (i)
c            rmin_rec    = x1b    (i)/cmpc
c         endif
c         if (tlambda(i) .lt. tmin_lambda) then
c            tmin_lambda = tlambda(i)
c            rmin_lambda = x1b    (i)/cmpc
c         endif
c      enddo

#ifdef MPI_USED
c
c Now find the smallest dt among all tiles, and send the result
c to all in buf_out.  This preserves solution concurrency among
c all the tiles
c
       buf_in(1) = dt
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                   , MPI_2FLOAT
     &                   , MPI_MINLOC, comm3d, ierr)
       dt  =   buf_out(1)
#endif /* MPI_USED */


c     Now subcycle the zone energies and species densities over the hydro timestep;
c     we advance the entire grid no further than the shortest chem/heating timestep.
c     


      tsubcycle = 0.

      do iter = 1 , itmax

        dtimin = huge

        do k = ks, ke

c#ifdef RT
          call mnu_dir_flx(k)
c#endif /* RT */

          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie
          tgas(i,j,k) = gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / (
#ifdef H
     .        HII(i,j,k)       + HI  (i,j,k)       + de   (i,j,k)
#endif /* H */
#ifdef He
     .      + HeI(i,j,k) * qrt + HeII(i,j,k) * qrt + HeIII(i,j,k) * qrt
#endif /* He */
#ifdef H2
     .      + HM (i,j,k)       + H2II(i,j,k) * haf + H2I  (i,j,k) * haf
#endif /* H2 */
     .                )
            logtem(i,j) = dlog(tgas(i,j,k))
            if (logtem(i,j) .lt. logtem0) logtem(i,j) = logtem0
            if (logtem(i,j) .ge. logtem9) logtem(i,j) = logtem9
            indixe(i,j) = min0(nratec-1,
     .                   max0(1,idint((logtem(i,j)-logtem0)/dlogtem)+1))
            t1    (i,j) = (logtem0 + (indixe(i,j) - 1)*dlogtem)
            t2    (i,j) = (logtem0 + (indixe(i,j)    )*dlogtem)
            tdef  (i,j) = t2(i,j) - t1(i,j)
#ifdef H
      ceHI   (i,j) = ceHIa   (indixe(i,j)) + (logtem(i,j)  -  t1(i,j))
     .      *(ceHIa   (indixe(i,j)+1)-ceHIa   (indixe(i,j)))/tdef(i,j)

      ciHI   (i,j) = ciHIa   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHIa   (indixe(i,j)+1)-ciHIa   (indixe(i,j)))/tdef(i,j)

      reHII  (i,j) = reHIIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHIIa  (indixe(i,j)+1)-reHIIa  (indixe(i,j)))/tdef(i,j)

      DM     (i,j) = DMa     (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(DMa     (indixe(i,j)+1)-DMa     (indixe(i,j)))/tdef(i,j)
#endif /* H */

#ifdef He
      ceHeI  (i,j) = ceHeIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ceHeIa  (indixe(i,j)+1)-ceHeIa  (indixe(i,j)))/tdef(i,j)

      ceHeII (i,j) = ceHeIIa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ceHeIIa (indixe(i,j)+1)-ceHeIIa (indixe(i,j)))/tdef(i,j)

      ciHeI  (i,j) = ciHeIa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeIa  (indixe(i,j)+1)-ciHeIa  (indixe(i,j)))/tdef(i,j)

      ciHeIS (i,j) = ciHeISa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeISa (indixe(i,j)+1)-ciHeISa (indixe(i,j)))/tdef(i,j)

      ciHeII (i,j) = ciHeIIa (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(ciHeIIa (indixe(i,j)+1)-ciHeIIa (indixe(i,j)))/tdef(i,j)

      reHeII1(i,j) = reHeII1a(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeII1a(indixe(i,j)+1)-reHeII1a(indixe(i,j)))/tdef(i,j)

      reHeII2(i,j) = reHeII2a(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeII2a(indixe(i,j)+1)-reHeII2a(indixe(i,j)))/tdef(i,j)

      reHeIII(i,j) = reHeIIIa(indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(reHeIIIa(indixe(i,j)+1)-reHeIIIa(indixe(i,j)))/tdef(i,j)
#endif /* He */

#ifdef H2
      if (iH2co .eq. 1) then
      hyd01k (i,j) = hyd01ka (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hyd01ka (indixe(i,j)+1)-hyd01ka (indixe(i,j)))/tdef(i,j)

      h2k01  (i,j) = h2k01a  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(h2k01a  (indixe(i,j)+1)-h2k01a  (indixe(i,j)))/tdef(i,j)

      vibh   (i,j) = vibha   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(vibha   (indixe(i,j)+1)-vibha   (indixe(i,j)))/tdef(i,j)

      roth   (i,j) = rotha   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(rotha   (indixe(i,j)+1)-rotha   (indixe(i,j)))/tdef(i,j)

      rotl   (i,j) = rotla   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(rotla   (indixe(i,j)+1)-rotla   (indixe(i,j)))/tdef(i,j)

      else if (iH2co .eq. 2) then
      gpldl  (i,j) = gpldla  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(gpldla  (indixe(i,j)+1)-gpldla  (indixe(i,j)))/tdef(i,j)

      gphdl  (i,j) = gphdla  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(gphdla  (indixe(i,j)+1)-gphdla  (indixe(i,j)))/tdef(i,j)
      endif
#endif /* H2 */

#ifdef HD
      hdlte  (i,j) = hdltea  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hdltea  (indixe(i,j)+1)-hdltea  (indixe(i,j)))/tdef(i,j)

      hdlow  (i,j) = hdlowa  (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(hdlowa  (indixe(i,j)+1)-hdlowa  (indixe(i,j)))/tdef(i,j)

#endif /* HD */
      brem   (i,j) = brema   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(brema   (indixe(i,j)+1)-brema   (indixe(i,j)))/tdef(i,j)

      cmpt   (i,j) = compt   (indixe(i,j)) +  (logtem(i,j)  -  t1(i,j))
     .      *(compt   (indixe(i,j)+1)-compt   (indixe(i,j)))/tdef(i,j)

#ifdef H
          k1 (i,j) = k1a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k1a (indixe(i,j)+1) -  k1a(indixe(i,j)))/tdef(i,j)
          k2 (i,j) = k2a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k2a (indixe(i,j)+1) -  k2a(indixe(i,j)))/tdef(i,j)
#endif /* H */
#ifdef He
          k3 (i,j) = k3a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k3a (indixe(i,j)+1) -  k3a(indixe(i,j)))/tdef(i,j)
          k4 (i,j) = k4a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k4a (indixe(i,j)+1) -  k4a(indixe(i,j)))/tdef(i,j)
          k5 (i,j) = k5a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k5a (indixe(i,j)+1) -  k5a(indixe(i,j)))/tdef(i,j)
          k6 (i,j) = k6a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k6a (indixe(i,j)+1) -  k6a(indixe(i,j)))/tdef(i,j)
#endif /* He */
#ifdef H2
          k7 (i,j) = k7a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k7a (indixe(i,j)+1) -  k7a(indixe(i,j)))/tdef(i,j)
          k8 (i,j) = k8a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k8a (indixe(i,j)+1) -  k8a(indixe(i,j)))/tdef(i,j)
          k9 (i,j) = k9a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k9a (indixe(i,j)+1) -  k9a(indixe(i,j)))/tdef(i,j)
          k10(i,j) = k10a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k10a(indixe(i,j)+1) - k10a(indixe(i,j)))/tdef(i,j)
          k11(i,j) = k11a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k11a(indixe(i,j)+1) - k11a(indixe(i,j)))/tdef(i,j)
          k12(i,j) = k12a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k12a(indixe(i,j)+1) - k12a(indixe(i,j)))/tdef(i,j)
          k14(i,j) = k14a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k14a(indixe(i,j)+1) - k14a(indixe(i,j)))/tdef(i,j)
          k16(i,j) = k16a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k16a(indixe(i,j)+1) - k16a(indixe(i,j)))/tdef(i,j)
          k17(i,j) = k17a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k17a(indixe(i,j)+1) - k17a(indixe(i,j)))/tdef(i,j)
          k18(i,j) = k18a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k18a(indixe(i,j)+1) - k18a(indixe(i,j)))/tdef(i,j)
          k19(i,j) = k19a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k19a(indixe(i,j)+1) - k19a(indixe(i,j)))/tdef(i,j)
          k22(i,j) = k22a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k22a(indixe(i,j)+1) - k22a(indixe(i,j)))/tdef(i,j)
#endif /* H2 */
#ifdef HD
          k50(i,j) = k50a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k50a(indixe(i,j)+1) - k50a(indixe(i,j)))/tdef(i,j)
          k51(i,j) = k51a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k51a(indixe(i,j)+1) - k51a(indixe(i,j)))/tdef(i,j)
          k52(i,j) = k52a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k52a(indixe(i,j)+1) - k52a(indixe(i,j)))/tdef(i,j)
          k53(i,j) = k53a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k53a(indixe(i,j)+1) - k53a(indixe(i,j)))/tdef(i,j)
          k54(i,j) = k54a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k54a(indixe(i,j)+1) - k54a(indixe(i,j)))/tdef(i,j)
          k55(i,j) = k55a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k55a(indixe(i,j)+1) - k55a(indixe(i,j)))/tdef(i,j)
          k56(i,j) = k56a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k56a(indixe(i,j)+1) - k56a(indixe(i,j)))/tdef(i,j)
#endif /* HD */
          enddo
          enddo

c     Now compute the heating and cooling rates of each zone   

          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie

c     Cooling Functions:

      edot(i,j,k)=(

#ifdef H
     . - ceHI   (i,j)* HI   (i,j,k)*de(i,j,k)             ! coll excit of HI
     . - ciHI   (i,j)* HI   (i,j,k)*de(i,j,k)             ! coll ioniz of HI
     . - reHII  (i,j)* HII  (i,j,k)*de(i,j,k)             ! recombin   of HII
     . - DM     (i,j)* HI   (i,j,k)**2                    ! DM ISM cooling
#endif /* H */
#ifdef He
     . - ceHeI  (i,j)* HeII (i,j,k)*de(i,j,k)**2*mhi*qrt  ! coll excit of HeI
     . - ceHeII (i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! coll excit of HeII
     . - ciHeI  (i,j)* HeI  (i,j,k)*de(i,j,k)*qrt         ! coll ioniz of HeI
     . - ciHeII (i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! coll ioniz of HeII
     . - ciHeIS (i,j)* HeII (i,j,k)*de(i,j,k)**2*mhi*qrt  ! coll ioniz of HeIS
     . - reHeII1(i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! recombin   of HeII
     . - reHeII2(i,j)* HeII (i,j,k)*de(i,j,k)*qrt         ! recombin   of HeII
     . - reHeIII(i,j)* HeIII(i,j,k)*de(i,j,k)*qrt         ! recombin   of HeIII
#endif /* He */
     . - brem   (i,j)*(
#ifdef H
     .                 HII  (i,j,k)                       ! elec brems cooling
#endif /* H */
#ifdef He
     .               + HeII (i,j,k)*qrt+HeIII(i,j,k)*qrt 
#endif /* He */
     .              )* de   (i,j,k)                       
     . - cmpt   (i,j)* de   (i,j,k)*mh                    ! compton cooling

     .            )*mhi**2

#ifdef H2
      if (iH2co .eq. 1) then
      qq        = 1.2*(HI(i,j,k)*mhi)**0.77 + (H2I(i,j,k)*mhi*haf)**0.77
      vibl      = (HI(i,j,k)*hyd01k(i,j) + H2I(i,j,k)*haf*h2k01(i,j))
     .          * mhi*8.18e-13

      edot(i,j,k) = edot(i,j,k) -                         ! ro-vibr cool by H2
     .         ( H2I(i,j,k)*(vibh(i,j)/(1.+vibh(i,j)/vibl)
     .         + roth(i,j)/(1.+roth(i,j)/(qq*rotl(i,j))))*mhi*haf)

      else if (iH2co .eq. 2) then

       gphdl1      = gphdl(i,j  ) / (HI(i,j,k) * mhi)
       edot(i,j,k) = edot (i,j,k) - H2I(i,j,k) *
     .               gphdl(i,j  ) / (1.0 + gphdl1 / gpldl(i,j))
     .             * haf * mhi
      endif

#endif /* H2 */
#ifdef HD
      if (iHDco .eq. 1) then

       hdlte1      = hdlte(i,j)/(HDI(i,j,k)*haf*mhi)
       hdlow1      = max(hdlow(i,j), tiny)
       edot(i,j,k) = edot(i,j,k) - HDI(i,j,k) *
     .               hdlte1/(1.0 + hdlte1/hdlow1)
     .             * haf * mhi
      endif
#endif /* HD */

#ifdef RT
c     Heating functions:

      edot(i,j,k) =   edot(i,j,k) + dfloat(ipiht)*(     ! photoioniz heating
#ifdef H
     .              + piHI  (i,j) * HI  (i,j,k)         ! pi of HI
#endif /* H */
#ifdef He
     .              + piHeI (i,j) * HeI (i,j,k)*qrt     ! pi of HeI
     .              + piHeII(i,j) * HeII(i,j,k)*qrt     ! pi of HeII
#endif /* He */
     .                                            ) * mhi

          enddo
          enddo
#endif /* RT */

          do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie
            dedot(i,j) = (
#ifdef H
     .              k1 (i,j) * HI   (i,j,k) * de (i,j,k)
#endif /* H  */
#ifdef He
     .            + k3 (i,j) * HeI  (i,j,k) * de (i,j,k)*qrt
     .            + k5 (i,j) * HeII (i,j,k) * de (i,j,k)*qrt
#endif /* He */
#ifdef H2
     .            + k8 (i,j) * HM   (i,j,k) * HI (i,j,k)
     .            + k15(i,j) * HM   (i,j,k) * HI (i,j,k)
     .            + k17(i,j) * HM   (i,j,k) * HII(i,j,k)
     .            + k14(i,j) * HM   (i,j,k) * de (i,j,k)
#endif /* H2 */
     .            ) * mhi * mhi - (
#ifdef H
     .              k2 (i,j) * HII  (i,j,k) * de (i,j,k)
#endif /* H  */
#ifdef He
     .            + k4 (i,j) * HeII (i,j,k) * de (i,j,k)*qrt
     .            + k6 (i,j) * HeIII(i,j,k) * de (i,j,k)*qrt
#endif /* He */
#ifdef H2
     .            + k7 (i,j) * HI   (i,j,k) * de (i,j,k)
     .            + k18(i,j) * H2II (i,j,k) * de (i,j,k)*haf
#endif /* H2 */
     .            ) * mhi * mhi
#ifdef RT
     .            + (
#ifdef H
     .                k24(i,j) * HI  (i,j,k)
#endif /* H  */
#ifdef He
     .              + k25(i,j) * HeII(i,j,k)*qrt
     .              + k26(i,j) * HeI (i,j,k)*qrt
#endif /* He */
     .              ) * mhi
#endif /* RT */


c     Compute the photoionization timescale for this zone

            if      (icycle .eq. 1) then
c              t_photo = dabs(0.1*de(i,j,k)/dedot(i,j)*mhi)
              t_photo = dabs(0.1*(de(i,j,k) + 1.0d-3*d(i,j,k))
     .                          /dedot(i,j)*mhi)
            else if (icycle .eq. 2) then
              t_photo = dabs(0.1*(de(i,j,k) + 1.0d-3*d(i,j,k))
     .                          /dedot(i,j)*mhi)

c               t_photo = 0.1 * dmin1((de(i,j,k) + 1.0d-3*d(i,j,k))/
c     .                               (dedot(i,j)*mhi),
c     .                               (HI(i,j,k) + 1.0d-3*d(i,j,k))/      
c     .                               (dedot(i,j)*mhi))

c              t_photo = dmax1(dabs(1.0d-06*HI(i,j,k)/dedot(i,j)*mhi),
c     .                        dabs(0.1*de(i,j,k)/dedot(i,j)*mhi))
c              t_photo = dmax1(dmin1(0.1*dx1a(i)/clight,
c     .                        dabs(0.1*HI(i,j,k)/dedot(i,j)*mhi)),
c     .                        dabs(0.1*de(i,j,k)/dedot(i,j)*mhi))
            endif

c     Compute the heating/cooling time for this zone

            if      (icycle .eq. 1) then
              t_heat  = 0.1 * dabs(e(i,j,k)/edot(i,j,k)) 
            else if (icycle .eq. 2) then
              if (edot(i,j,k) .ge. 0.) then
                t_heat  = dabs(e(i,j,k)/edot(i,j,k)) 
              else
                t_heat  = 0.1 * dabs(e(i,j,k)/edot(i,j,k)) 
              endif
            endif


            t_chem = dmin1(t_photo,t_heat)
            dtimin = dmin1(t_chem ,dtimin)

            
          enddo
          enddo

          
        enddo


        dtimin    = dmin1(dtimin,dt-tsubcycle)
 66     format(7(1pe11.3e3,1x))

        do k = ks, ke

c#ifdef RT
         call mnu_dir_flx(k)
c#endif /* RT */

         do j = js, je
CDIR$ UNROLL UNROLL_I
          do i = is, ie

            logtem(i,j) = dlog(tgas(i,j,k))
            if (logtem(i,j) .lt. logtem0) logtem(i,j) = logtem0
            if (logtem(i,j) .ge. logtem9) logtem(i,j) = logtem9
            indixe(i,j) = min0(nratec-1,
     .                   max0(1,idint((logtem(i,j)-logtem0)/dlogtem)+1))
            t1    (i,j) = (logtem0 + (indixe(i,j) - 1)*dlogtem)
            t2    (i,j) = (logtem0 + (indixe(i,j)    )*dlogtem)
            tdef  (i,j) = t2(i,j) - t1(i,j)
#ifdef H
          k1 (i,j) = k1a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k1a (indixe(i,j)+1) -  k1a(indixe(i,j)))/tdef(i,j)
          k2 (i,j) = k2a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k2a (indixe(i,j)+1) -  k2a(indixe(i,j)))/tdef(i,j)
#endif /* H */
#ifdef He
          k3 (i,j) = k3a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k3a (indixe(i,j)+1) -  k3a(indixe(i,j)))/tdef(i,j)
          k4 (i,j) = k4a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k4a (indixe(i,j)+1) -  k4a(indixe(i,j)))/tdef(i,j)
          k5 (i,j) = k5a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k5a (indixe(i,j)+1) -  k5a(indixe(i,j)))/tdef(i,j)
          k6 (i,j) = k6a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k6a (indixe(i,j)+1) -  k6a(indixe(i,j)))/tdef(i,j)
#endif /* He */
#ifdef H2
          k7 (i,j) = k7a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k7a (indixe(i,j)+1) -  k7a(indixe(i,j)))/tdef(i,j)
          k8 (i,j) = k8a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k8a (indixe(i,j)+1) -  k8a(indixe(i,j)))/tdef(i,j)
          k9 (i,j) = k9a (indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k9a (indixe(i,j)+1) -  k9a(indixe(i,j)))/tdef(i,j)
          k10(i,j) = k10a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k10a(indixe(i,j)+1) - k10a(indixe(i,j)))/tdef(i,j)
          k11(i,j) = k11a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k11a(indixe(i,j)+1) - k11a(indixe(i,j)))/tdef(i,j)
          k12(i,j) = k12a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k12a(indixe(i,j)+1) - k12a(indixe(i,j)))/tdef(i,j)
          k14(i,j) = k14a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k14a(indixe(i,j)+1) - k14a(indixe(i,j)))/tdef(i,j)
          k16(i,j) = k16a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k16a(indixe(i,j)+1) - k16a(indixe(i,j)))/tdef(i,j)
          k17(i,j) = k17a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k17a(indixe(i,j)+1) - k17a(indixe(i,j)))/tdef(i,j)
          k18(i,j) = k18a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k18a(indixe(i,j)+1) - k18a(indixe(i,j)))/tdef(i,j)
          k19(i,j) = k19a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k19a(indixe(i,j)+1) - k19a(indixe(i,j)))/tdef(i,j)
          k22(i,j) = k22a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k22a(indixe(i,j)+1) - k22a(indixe(i,j)))/tdef(i,j)
#endif /* H2 */
#ifdef HD
          k50(i,j) = k50a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k50a(indixe(i,j)+1) - k50a(indixe(i,j)))/tdef(i,j)
          k51(i,j) = k51a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k51a(indixe(i,j)+1) - k51a(indixe(i,j)))/tdef(i,j)
          k52(i,j) = k52a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k52a(indixe(i,j)+1) - k52a(indixe(i,j)))/tdef(i,j)
          k53(i,j) = k53a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k53a(indixe(i,j)+1) - k53a(indixe(i,j)))/tdef(i,j)
          k54(i,j) = k54a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k54a(indixe(i,j)+1) - k54a(indixe(i,j)))/tdef(i,j)
          k55(i,j) = k55a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k55a(indixe(i,j)+1) - k55a(indixe(i,j)))/tdef(i,j)
          k56(i,j) = k56a(indixe(i,j))   + (logtem    (i,j) - t1  (i,j))
     .             *(k56a(indixe(i,j)+1) - k56a(indixe(i,j)))/tdef(i,j)
#endif /* HD */

            e(i,j,k) = e(i,j,k) + edot(i,j,k) * dtimin

c            rad_press(i,j,k) = rad_press(i,j,k) +
c     .                         f(i,j,1) * sigma24(1)*0.002*evhz/clight*
c     .                         0.5 * (HI(i-1,j,k) + HI(i,j,k)) * mhi *
c     .                         dtimin


#ifdef H
            scoef  =   (k2 (i,j) * HII (i,j,k) * de (i,j,k) 
#ifdef H2
     .             +    k13(i,j) * HI  (i,j,k) * H2I(i,j,k)
     .             +    k11(i,j) * HII (i,j,k) * H2I(i,j,k)*haf
     .             +    k12(i,j) * de  (i,j,k) * H2I(i,j,k)
     .             +    k14(i,j) * HM  (i,j,k) * de (i,j,k)    
     .             +    k15(i,j) * HM  (i,j,k) * HI (i,j,k)    
     .             + 2.*k16(i,j) * HM  (i,j,k) * HII(i,j,k)    
     .             +    k18(i,j) * H2II(i,j,k) * de (i,j,k)
     .             +    k19(i,j) * H2II(i,j,k) * HM (i,j,k)*haf 
#endif /* H2 */
     .                 ) * mhi 
#if defined RT && defined H2
     .             +    k31(i,j) * H2I (i,j,k)             
#endif /* RT && H2 */

            acoef  =   (k1 (i,j) * de  (i,j,k) 
#ifdef H2
     .             +    k7 (i,j) * de  (i,j,k)   
     .             +    k8 (i,j) * HM  (i,j,k) 
     .             +    k9 (i,j) * HII (i,j,k) 
     .             +    k10(i,j) * H2II(i,j,k)    * haf
     .             + 2.*k22(i,j) * HI  (i,j,k)**2 * mhi
#endif /* H2 */
     .                 ) * mhi 
#ifdef RT
     .             + k24(i,j)
#endif /* RT */
            HIp    = ( scoef*dtimin + HI(i,j,k)) 
     .             / ( 1. + acoef*dtimin )

            scoef  =   (k1 (i,j) * HIp         * de (i,j,k)    
#ifdef H2
     .             +    k10(i,j) * H2II(i,j,k) * HIp * haf
#endif /* H2 */
     .                 ) * mhi 
#ifdef RT
     .             + k24(i,j) * HIp
#endif /* RT */
            acoef  =   (k2 (i,j) * de  (i,j,k) 
#ifdef H2
     .             +    k9 (i,j) * HIp         
     .             +    k11(i,j) * H2I (i,j,k) * haf
     .             +    k16(i,j) * HM  (i,j,k) 
     .             +    k17(i,j) * HM  (i,j,k) 
#endif /* H2 */
     .                 ) * mhi 
            HIIp   = ( scoef*dtimin + HII(i,j,k)) 
     .             / ( 1. + acoef*dtimin) 

            HI   (i,j,k) = HIp
            HII  (i,j,k) = HIIp
#endif /* H */

#ifdef He
            scoef  = k4(i,j)*HeII(i,j,k)*de(i,j,k) * mhi
            acoef  = k3(i,j)*de(i,j,k)*mhi
#ifdef RT
     .             + k26(i,j)
#endif /* RT */
            HeIp   = ( scoef*dtimin + HeI(i,j,k))
     .             / ( 1. + acoef*dtimin) 

            scoef  = (k3(i,j) * HeIp         * de(i,j,k)
     .             +  k6(i,j) * HeIII(i,j,k) * de(i,j,k)
     .               ) * mhi
#ifdef RT
     .             + k26(i,j) * HeIp
#endif /* RT */
            acoef  = ( (k4(i,j) + k5(i,j))   * de(i,j,k)) * mhi
#ifdef RT
     .             + k25(i,j)
#endif /* RT */
            HeIIp  = ( scoef*dtimin + HeII(i,j,k)) 
     .             / ( 1. + acoef*dtimin ) 

            scoef   = k5(i,j)*HeIIp*de(i,j,k) * mhi
#ifdef RT
     .             + k25(i,j)*HeIIp
#endif /* RT */
            acoef   = k6(i,j)*de(i,j,k)*mhi
            HeIIIp  = ( scoef*dtimin + HeIII(i,j,k)) 
     .              / ( 1. + acoef*dtimin) 

            HeI  (i,j,k) = HeIp
            HeII (i,j,k) = HeIIp
            HeIII(i,j,k) = HeIIIp
#endif /* He */

#if defined H || defined He
            scoef =  0.
#ifdef H2
     .            + (k8 (i,j) * HM(i,j,k) * HIp  
     .            +  k15(i,j) * HM(i,j,k) * HIp 
     .            +  k17(i,j) * HM(i,j,k) * HIIp 
     .              ) * mhi 
#endif /* H2 */
#ifdef RT
#ifdef H
     .            + k24(i,j) * HIp   
#endif /* H */
#ifdef He
     .            + k25(i,j) * HeIIp      * qrt
     .            + k26(i,j) * HeIp       * qrt
#endif /* He */
#endif /* RT */

            acoef = - ( 
#ifdef H
     .                  k1 (i,j)*HIp       - k2(i,j)*HIIp
#endif /* H */
#ifdef He
     .              +   k3 (i,j)*HeIp *qrt - k6(i,j)*HeIIIp*qrt
     .              +   k5 (i,j)*HeIIp*qrt - k4(i,j)*HeIIp *qrt
#endif /* He */
#ifdef H2
     .              +   k14(i,j)*HM(i,j,k)
     .              -   k7 (i,j)*HIp
     .              -   k18(i,j)*H2II(i,j,k)*haf
#endif /* H2 */
     .                ) * mhi
            dep   = ( scoef*dtimin + de(i,j,k)) 
     .            / ( 1. + acoef*dtimin) 
  
            de   (i,j,k) = dep
#endif /* H or He */

#ifdef H2
            HM(i,j,k) = ( k7(i,j)*HI(i,j,k)*de(i,j,k) )
     .                / ( (k8 (i,j)+k15(i,j))*HI (i,j,k)
     .           + ( k16(i,j)+k17(i,j))*HII(i,j,k)+k14(i,j)*de(i,j,k)
#ifdef RT
     .           + k27(i,j)*mh
#endif /* RT */
     .                  )

            H2II(i,j,k) = 2.*( k9 (i,j)*HI(i,j,k)*HII(i,j,k)
     .                       + k11(i,j)*H2I(i,j,k)*haf*HII(i,j,k)
     .                       + k17(i,j)*HM(i,j,k)*HII(i,j,k)
#ifdef RT
     .           + k29(i,j)*H2I(i,j,k)*mh
#endif /* RT */
     .                       )
     .                     / ( k10(i,j)*HI (i,j,k) + k18(i,j)*de(i,j,k)
     .                       + k19(i,j)*HM(i,j,k)
#ifdef RT
     .           + (k28(i,j)+k30(i,j))*mh
#endif /* RT */
     .                       )

            scoef = 2.*( k8 (i,j)*HM  (i,j,k)*HI(i,j,k)
     .            + k10(i,j)*H2II(i,j,k)*HI(i,j,k)*haf
     .            + k19(i,j)*H2II(i,j,k)*HM(i,j,k)*haf 
     .            + k22(i,j)*HI  (i,j,k)**3 * mhi ) * mhi          
            acoef = ( k13(i,j)*HI(i,j,k) + k11(i,j)*HII(i,j,k)
     .              + k12(i,j)*de(i,j,k) )*mhi
#ifdef RT
     .              + k29(i,j) + k31(i,j)
#endif /* RT */

            H2I(i,j,k) = ( scoef*dtimin + H2I(i,j,k)) 
     .                 / ( 1. + acoef*dtimin ) 
#endif /* H2 */

#ifdef HD
            scoef = (   k2 (i,j) * DII(i,j,k) * de(i,j,k)
     .            +     k51(i,j) * DII(i,j,k) * HI(i,j,k)
     .            + 2.* k55(i,j) * HDI(i,j,k) * HI(i,j,k) * thd
     .            ) * mhi 
            acoef      (k1 (i,j) * de (i,j,k)
     .            +     k50(i,j) * HII(i,j,k)
     .            +     k54(i,j) * H2I(i,j,k) * haf
     .            +     k56(i,j) * HM (i,j,k)
     .            +     k22(i,j) * HI (i,j,k)**2 * mhi * haf
     .            +     k57(i,j) * HI (i,j,k) * H2I(i,j,k) * mhi * haf
     .            ) * mhi 
#ifdef RT
     .            +     k24(i,j)
#endif /* RT */
            DIp   = ( scoef*dtimin + DI(i,j,k) ) 
     .              ( 1. + acoef*dtimin )

            scoef =    (k1 (i,j) * DIp        * de (i,j,k) 
     .            +     k50(i,j) * HII(i,j,k) * DIp 
     .            + 2.* k53(i,j) * HII(i,j,k) * HDI(i,j,k) * thd
     .            ) * mhi
#ifdef RT
     .            +     k24(i,j) * DIp
#endif /* RT */
            acoef =    (k2 (i,j) * de (i,j,k)
     .            +     k51(i,j) * HI (i,j,k)
     .            +     k52(i,j) * H2I(i,j,k) * haf
     .            ) * mhi

            DIIp  = ( scoef*dtimin + DII(i,j,k) )
     .            / ( 1. + acoef*dtimin )

            scoef = 3.*(k52(i,j) * DIIp * H2I(i,j,k) * qrt
     .            +     k54(i,j) * DIp  * H2I(i,j,k) * qrt
     .            +     k56(i,j) * DIp  * HM (i,j,k) 
     .            +     k22(i,j) * DIp*HI(i,j,k)**2  * mhi * qrt
     .            +     k57(i,j) * DIp*HI(i,j,k)*H2I(i,j,k)*mhi*qrt
     .                 ) * mhi
            acoef =    (k53(i,j) * HII(i,j,k)
     .            +     k55(i,j) * HI(i,j,k)
     .                 ) * mhi
c
            HDIp  = ( scoef*dtimin + HDI(i,j,k) )
     .            / ( 1. + acoef*dtimin )

            DI (i,j,k) = DIp
            DII(i,j,k) = DIIp
            HDI(i,j,k) = HDIp
#endif /* HD */


          enddo
          enddo

        enddo
  
        tsubcycle = tsubcycle + dtimin
        if (dt - tsubcycle .lt. 0.001*dt) goto 1000
        
        
        enddo

1000  continue

      call mconsrv

c
c  Mark multispecies density boundaries out of date
c

       do 40 i=1,6
#ifdef H
         bvstat(i,15) = 0  !  HII
#endif /* H */
#ifdef He
         bvstat(i,16) = 0  !  HI
         bvstat(i,17) = 0  !  de
         bvstat(i,18) = 0  !  HeI
         bvstat(i,19) = 0  !  HeII
         bvstat(i,20) = 0  !  HeIII
#endif /* He */
#ifdef H2
         bvstat(i,21) = 0  !  HM
         bvstat(i,22) = 0  !  H2I
         bvstat(i,23) = 0  !  H2II
#endif /* H2 */
40     continue

#ifdef MPI_USED
      call MPI_BARRIER(comm3d, ierr)
#endif /* MPI_USED */

 2000 return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////            3 D C O O L C H E M            \\\\\\\\\\
c
c=======================================================================
c
c
