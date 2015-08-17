!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              T E X T D M P 2              \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine textdmp2
!
!  written by:   Robert Fiedler
!  Last modified: 01-09-07 by Daniel Whalen for inclusion in ZEUS-MP
!                 2.1
!
!  PURPOSE:  Write out ASCII data dumps of hydrodynamical, chemistry,
!	     and radiation variables as the user specifies in a variety
!            of possible formats.
!
!
!-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
      use chem
      use scratch
      use opac_law
      use mpiyes
      use mpipar
!
      implicit NONE
      integer  :: i, j, k
      real(rl) :: mhi,cnorm,e_kinetic
      mhi = 1. / mh
!
! Assign unit 12 to the ASCII data file whose name is contained in the
! string variable 'usrfile'.
!
      open(12,file=usrfile)
 1    format(   1pe11.3e3    )
 2    format( 2(1pe11.3e3,1x))
 3    format( 3(1pe11.3e3,1x))
 4    format( 4(1pe11.3e3,1x))
! 4    format( 4(1pe16.8e3,1x))
 5    format( 5(1pe11.3e3,1x))
 6    format( 6(1pe11.3e3,1x))
 7    format( 7(1pe11.3e3,1x))
 8    format( 8(1pe13.6e2,1x))
 9    format( 9(1pe11.3e3,1x))
10    format(10(1pe11.3e3,1x))
11    format(11(1pe11.3e3,1x))
12    format(12(1pe16.8e3,1x))
!12    format(12(1pe11.3e3,1x))
13    format(13(1pe11.3e3,1x))
14    format(14(1pe11.3e3,1x))
 15   format(15(1pe16.8e3,1x))
! 
! Write out the data in the desired format
!
!--------------------------1D Sod_CMA outputs---------------------------
!
!      write(12,"('    x1b          X1          X2          X7        ', 
!     .'  X8          X9          X4          X5          X6       
!     .  HSum       HeSum       ChgSum')")
!      if      (ldimen .eq. 1) then
!      write(12,12)(((x1a(i),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
!     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
!     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
!     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
!     .   abun(i,j,k,9),
!     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
!     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
!     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
!     . , i=is,ie), j=js,je), k=ks,ke)
!      else if (ldimen .eq. 2) then
!      write(12,12)(((x2a(j),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
!     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
!     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
!     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
!     .   abun(i,j,k,9),
!     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
!     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
!     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
!     . , i=is,is), j=js,je), k=ks,ks)
!      else if (ldimen .eq. 3) then
!      write(12,12)(((x3a(k),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
!     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
!     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
!     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
!     .   abun(i,j,k,9),
!     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
!     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
!     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
!     . , i=is,is), j=js,js), k=ks,ke)
!      endif
!-----------------------Sedov-Taylor blast output-----------------------
      cnorm   = boltz*tevk*mhi
      write(12,1) time/3.1536e07
      write(12,*) nhy,dt    
      write(12,"('     x1b        dens       energ       metal        v1 &
      ')")
      write(12,"('     x1b        dens       energ')")
      write(12,3) ((( x1b(i)/cmpc, &
        d(i,j,k), &
       (abun(i,j,k,1)     + abun(i,j,k,2)     + &
        abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + &
        abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)* &
        d   (i,j,k)  *mhi, &
        e(i,j,k), &
        tgas(i,j,k) , &
        gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / &
      ((abun(i,j,k,1)       + abun(i,j,k,2)       + abun(i,j,k,3) &
      + abun(i,j,k,4) * qrt + abun(i,j,k,5) * qrt + abun(i,j,k,6) * qrt &
      )*d(i,j,k)), &
        d(i,j,k) * abun(i,j,k,10), &
        v1(i,j,k)/cmkm, &
       7.6d-1*(1.-abun(i,j,k,10))-abun(i,j,k,1)-abun(i,j,k,2), &
       2.4d-1*(1.-abun(i,j,k,10))-abun(i,j,k,4)-abun(i,j,k,5) &
                 -abun(i,j,k,6), &
      fh*(1-abun(i,j,k,10))- abun(i,j,k,1) - abun(i,j,k,2), &
      (1-fh)*(1-abun(i,j,k,10)) - abun(i,j,k,4)-abun(i,j,k,5) &
      -abun(i,j,k,6), &
      abun(i,j,k,1),abun(i,j,k,2),abun(i,j,k,3),abun(i,j,k,4), &
      abun(i,j,k,5),abun(i,j,k,6),abun(i,j,k,10) &
      , i= is,ie+2), j= js,js), k= ks,ks)
!-----------------------6-species radshock output-----------------------
!
!
!      e_kinetic = 0.
!
!      write(12,1) time/3.1536e07
!c      write(12,*) nhy,dt    
!c      write(12,"('     x1b       dens         energ       tgas       HI      
!c     .',
!c     .             '        HII        elec       He        He+        
!c     .  He++')")
!      write(12,13) ((( x1b(i)/cmpc, d(i,j,k), 
!     .  e(i,j,k) ,tgas(i,j,k)
!     ., abun(i,j,k,1) ! / (abun(i,j,k,1) + abun(i,j,k,2))  
!     ., abun(i,j,k,2) ! / (abun(i,j,k,1) + abun(i,j,k,2))
!     ., abun(i,j,k,3) ! *  d   (i,j,k  )
!     ., abun(i,j,k,4) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
!     ., abun(i,j,k,5) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
!     ., abun(i,j,k,6) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
!     ., abun(i,j,k,7) 
!     ., abun(i,j,k,8) 
!     ., abun(i,j,k,9)
!c     .,(gp(i,j,k)-gp(i-1,j,k))*dx1bi(i)*dt   
!c     ., gp(i,j,k)
!     ., i= is,ie), j= js,js), k= ks,ks)
!
!      write(12,"(a)")
!
!      do k=ks,ke
!       do j=js,je
!        do i=is,ie
!         e_kinetic = e_kinetic + haf * qrt * (v1(i,j,k)+v1(i+1,j,k))**2
!     .             * 4. * pi * thd * d(i,j,k) * (x1a(i+1)**3-x1a(i)**3) 
!        enddo
!       enddo
!      enddo
!         
!      write(12,9) ceHItot  ,ciHItot  ,ceHeItot ,ciHeItot ,ceHeIItot,
!     .            ciHeIItot,bremtot  ,ICtot    ,e_kinetic
!
!
!------------------------1D advection test output------------------------
!
!      write(12,1) time
!      write(12,*) nhy,dt    
!      write(12,"('     x1b       dens       energ       tgas       HI      
!     .',
!     .             '        HII        elec       He        He+        
!     .  He++')")
!      write(12,10) ((( x1b(i), d(i,j,k), e(i,j,k) ,tgas(i,j,k)
!     ., abun(i,j,k,1) * d(i,j,k), abun(i,j,k,2) * d(i,j,k)
!     ., abun(i,j,k,3) * d(i,j,k), abun(i,j,k,4) * d(i,j,k) 
!     ., abun(i,j,k,5) * d(i,j,k), abun(i,j,k,6) * d(i,j,k) 
!     ., i= is,ie), j= js,js), k= ks,ks)
!----------------------------2D GSF outputs-----------------------------
!
!      write(12,1) time/3.1536e07
!      write(12,*) nhy,dt    
!      write(12,*) iter
!      write(12,"('    x1b         dens       energ        tgas')")
!      write(12,4) (((x1b(i)/cmpc, d(i,j,k)/mh, e(i,j,k), tgas(i,j,k) 
!     ., i=is,ie), j=js,je), k=ks,ke)
!-------------------3-species hydro outputs with flux-------------------
!
!      write(12,1) time/3.1536e07
!      write(12,*) nhy,dt    
!      write(12,*) iter
!      write(12,"('    x1b         dens        enrg        tgas       ', 
!     .'  HI         HII         elec         k24        flux')")
!      write(12,9) (((x1b(i)/cmpc 
!     .,(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, e(i,j,k)
!     .,d(i,j,k) * mhi, e(i,j,k)
!     ., tgas(i,j,k  ) 
!     ., abun(i,j,k,1) * d(i,j,k) * mhi, abun(i,j,k,2) * d(i,j,k) * mhi
!     ., abun(i,j,k,3) * d(i,j,k) * mhi, k24(i,j), f(i,j,1)
!     ., i=is,ie), j=js,je), k=ks,ke)
!---------------------3-species network diagnostics---------------------
!
!      write(12,1) time/3.1536e07
!      write(12,*) nhy,dt    
!      write(12,*) iter
!      write(12,"('    x1b          k1         k2         k24       ', 
!     .' flux ')")
!      write(12,5) (((x1b(i)/cmpc 
!     ., k1(i,j) , k2(i,j) , k24(i,j), f(i,j,1)
!     ., i=is,ie), j=js,je), k=ks,ke)
!-----------------hydro outputs with n_strom (paper I)------------------
!
!      write(12,2) time/3.1536e07
!      write(12,*) nhy,dt    
!      write(12,*) iter
!      write(12,"('     x1b       dens       energ       tgas       HI      
!     .',
!     .             '        HII        elec       k24        n_strom')")
!      write(12,9) (((x1b(i)/cmpc 
!     .,(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, e(i,j,k) 
!     ., tgas(i,j,k  )
!     ., abun(i,j,k,1) * d(i,j,k) * mhi, abun(i,j,k,2) * d(i,j,k) * mhi
!     ., abun(i,j,k,3) * d(i,j,k) * mhi, k24(i,j)
!     ., 5.470D30 * x1b(i)**(-1.5)
!     ., i=is,ie), j=js,js), k=ks,ks)
!
!---------------------------timescale outputs---------------------------
!
!      write(12,1) time/3.1536e07
!      write(12,*) nhy,dt, iter    
!      write(12,*) itmcool, itmrate
!      write(12,"('     x1b        dens        tgas        elec       thy      
!     .dro       tpdv        tphi        trec      tlambda')")
!      write(12,9) (((x1b(i)/cmpc 
!      write(12,9) (((dfloat(i), dfloat(j) , dfloat(k), 
!     ., (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, tgas(i,j,k)
!     .,  abun(i,j,k,3) * d(i,j,k) * mhi
!     .,  thydro(i), tpdv(i), tphi(i), trec(i), tlambda(i)
!     .,i=is,ie), j=js,je), k=ks,ke)
!------------------Toronto T1 - T7 conference outputs-------------------
!
!      if (usrtag .eq. 'T1_') then
!        write(12,2) (((x1b(i)/cmkpc, abun(i,j,k,1),
!     .                 i=is,ie), j=js,js), k=ks,ks)
!        goto 68
!
!      else if (usrtag .eq. 'T2_') then
!        write(12,6) (((x1b(i)/cmkpc, abun(i,j,k,1), tgas(i,j,k),
!     .                 i=is,ie), j=js,js), k=ks,ks)
!        goto 68
!
!      else if (usrtag .eq. 'T5_') then
!        write(12,5) (((x1b(i)/cmkpc, (abun(i,j,k,1) + abun(i,j,k,2))
!     .                 *mhi*d(i,j,k),
!     .                 abun(i,j,k,1), gamm1 * e(i,j,k), tgas(i,j,k),
!     .                 i=is,ie), j=js,js), k=ks,ks)
!        goto 68
!
!      else if (usrtag .eq. 'T6_') then
!       write(12,4) (((x1b(i)/cmkpc, abun(i,j,k,1),
!     .(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, tgas(i,j,k),
!     .                 i=is,ie), j=js,je), k=ks,ke)
!        goto 68
!
!      else if (usrtag .eq. 'T7_') then
!c        write(12,5) (((x1b(i)/cmkpc,
!c     .                (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi,
!c     .                 abun(i,j,k,1), gamm1 * e(i,j,k), tgas(i,j,k),
!c     .                 i=is,ie), j=js,je), k=ks,ke)
!        write(12,1) (((dsqrt(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2)/
!     .                 dsqrt(gamma*gamm1*e(i,j,k)/d(i,j,k))
!     .                , i=is,ie ) , j=js,je ), k=ks,ke )
!        goto 68
!      endif
!
!-----------------------3D I-front hydro outputs------------------------
!      write(12,*) time/3.1536e07, iter
!      write(12,"('     x1b         x2b         x3b        dens        
!     .energ       tgas         HI         HII         k24         
!     .flx')")
!      write(12,10)(((dfloat(i),dfloat(j),dfloat(k), 
!      write(12,10) (((x1b(i)/cmpc, x2b(j), x3b(k), 
!     .               (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi,
!     .                e(i,j,k), tgas(i,j,k), abun(i,j,k,1)*d(i,j,k)*mhi, 
!     .                abun(i,j,k,2)*d(i,j,k)*mhi, k24(i,j), f(i,j,1),
!     .                i=is,ie), j=js,je), k=ks,ke)
!-----------------------2D halo collapse outputs------------------------
!      write(12,"('     x1b         x2b         x3b        dens        
!     .tgas')")
!      write(12,5) (((x1b(i)/cmpc, x2b(j)/cmpc, x3b(k), 
!     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
!     .  d   (i,j,k)  *mhi , tgas(i,j,k),
!     .                i=is,ie), j=js,je), k=ks,ke)
!-----------------------6-species network outputs-----------------------
!      write(12,*) time/3.1536e07, iter
!      write(12,"('    x1b         dens        tgas         HI         
!     .HII         HeI         HeII        HeIII ')")
!      write(12,8) (((x1b(i)/cmpc, 
!     . (abun(i,j,k,1) + abun(i,j,k,2) + (abun(i,j,k,4)  + 
!     .  abun(i,j,k,5) + abun(i,j,k,6))* qrt) * d(i,j,k) * mhi,
!     .  tgas(i,j,k),
!     .  abun(i,j,k,1)*d(i,j,k)*mhi     , abun(i,j,k,2)*d(i,j,k)*mhi  
!     .  abun(i,j,k,4)*d(i,j,k)*mhi*qrt , abun(i,j,k,5)*d(i,j,k)*mhi*qrt, 
!     .  abun(i,j,k,6)*d(i,j,k)*mhi*qrt , 
!     .i=is,ie ) , j=js,js ), k=ks,ks )
!--------------------tgas and dens outputs for rjw99--------------------
!
!      write(12,15) (((tgas(i,j,k), d(i,j,k)/mh 
!     .                 ,i=is,ie ) , j=js,je ), k=ks,ke )
!
!-----------------------2D acc cutoff diagnostics-----------------------
!      mhi      = 1./mh
!      write(12,9) (((x1b(i),x2a(j),x3a(k),
!     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
!     .  d   (i,j,k)  *mhi , abun(i,j,k,1), abun(i,j,k,2), 
!     .  v1(i,j,k)/cmkm, v2(i,j,k)/cmkm, v3(i,j,k)/cmkm
!     .                 ,i=is,ie ) , j=js,je ), k=ks,ke )
!-----------------------9-species network outputs-----------------------
      mhi      = 1./mh
      cnorm    = boltz*tevk*mhi
!c      write(12,*) time/3.1536e07
!c      write(12,"('    x1b         dens        tgas         HI       ', 
!c     .'  HII         HeI         HeII        HeIII        HM          H2I
!c     .         H2II')")
!      write(12,1) time/3.1536e07
!      write(12,15) (((x1b(i),
!     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
!     .  d   (i,j,k)  *mhi ,
!     .  d   (i,j,k),
!     .  e   (i,j,k)/d   (i,j,k) ,
!     .  v1  (i,j,k),
!c     .  gamm1 * dabs(e(i,j,k)),
!     .  tgas (i,j,k),
!c     .  gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / (1.66052992e-22 * (
!c     .  abun(i,j,k,1)       + abun(i,j,k,2)       + abun(i,j,k,3)
!c     .+ abun(i,j,k,4) * qrt + abun(i,j,k,5) * qrt + abun(i,j,k,6) * qrt
!c     .+ abun(i,j,k,7)       + abun(i,j,k,8) * haf + abun(i,j,k,9) * haf
!c     .                )),
!     . abun(i,j,k,1), ! /(abun(i,j,k,1)+abun(i,j,k,2)+2.*abun(i,j,k,8)),
!c     . (abun(i,j,k,2)     + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt)/
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf), 
!c     .  abun(i,j,k,1),!     /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,2),!/
!     .  abun(i,j,k,3),!/
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,4),!*qrt /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,5),!*qrt /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,6),!*qrt /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,7),!     /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!c     .  gamm1 * e(i,j,k),
!     .  abun(i,j,k,8),!*haf /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!     .  abun(i,j,k,9),!*haf /
!c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
!c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
!c     .  gp(i,j,k), 
!     .  i=is,ie), j=js,js ), k=ks,ks )
!-----------------9-species network diagnostic outputs------------------
!      write(12,*) time/3.1536e07, iter
!      write(12,"('    x1b         k9         k11        k17         k29         
!     .k10         k18         k19         k28        k30')")
!      write(12,10) (((x1b(i)/cmpc, 
!     .  k9 (i,j)*abun(i,j,k,1)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi    ,
!     .  k11(i,j)*abun(i,j,k,8)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi*haf,  
!     .  k17(i,j)*abun(i,j,k,7)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi    ,
!     .  k29(i,j)*abun(i,j,k,8)*              d(i,j,k)   *mhi    *haf,  
!     .  k10(i,j)*abun(i,j,k,9)*abun(i,j,k,1)*d(i,j,k)**2*mhi*mhi*haf,
!     .  k18(i,j)*abun(i,j,k,9)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi*haf,  
!     .  k19(i,j)*abun(i,j,k,9)*abun(i,j,k,7)*d(i,j,k)**2*mhi*mhi*haf,
!     .  k28(i,j)*abun(i,j,k,9)*              d(i,j,k)   *mhi    *haf,  
!     .  k30(i,j)*abun(i,j,k,9)*              d(i,j,k)   *mhi    *haf, 
!     .i=is,ie), j=js,js), k=ks,ks)
!----------------------3D halo reader test outputs----------------------
!      write(12,8) (((d(i,j,k), v1(i,j,k), v2(i,j,k), v3(i,j,k),  
!     .               e(i,j,k), x1b(i), x2b(j), x3b(k),
!     .               i=is,ie), j=js,je), k=ks,ke)
!--------------------multispecies cooling rate dumps--------------------
!      write(12,"(a)")
!      write(12,"('     x1b        ceco        cico       HH/DM    ',   
!     .      '    cmpt         H2         brem        recH')")
!      write(12,8) ((( x1b(i)/cmpc, 
!     .  ceHI (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
!     .  ciHI (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
!     .  DM   (i,j)*abun(i,j,k,1)*abun(i,j,k,1)*d(i,j,k)**2*mhi*mhi,
!     .  cmpt (i,j)*abun(i,j,k,3)              *d(i,j,k)       *mhi,
!     .  abun (i,j,k,8)*d(i,j,k) * gphdl(i,j)/(1.0+gphdl1/gpldl(i,j))
!     .               * haf * mhi,
!     .  brem  (i,j) * (
!     .  ) * abun(i,j,k,3)*d(i,j,k)*d(i,j,k)*mhi*mhi,
!     .  reHII(i,j)*abun(i,j,k,2)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
!     .i=is,ie), j=js,js), k=ks,ks)
!      write(12,"('     x1b        k1          k2        k24',  
!     .'       dedot   de/dedot')")
!      write(12,6) ((( x1a(i), 
!     .  k1   (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
!     .- k2   (i,j)*abun(i,j,k,2)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
!     .  k24  (i,j)*abun(i,j,k,1)              *d(i,j,k)       *mhi,
!     .  dedot(i,j),abun(i,j,k,3)*d(i,j,k)*mhi/dedot(i,j),
!     .i=is,ie), j=js,js), k=ks,ks)
!--------------------------------velocities-----------------------------
 68   write(12,"(a)")
!      write(12,"('         x1a             x2b             x3b',
!     1 '              v1')")
      write(12,2) ( ( ( x1a(i)/cmpc, v1(i,j,k)/cmkm, & !dv_ioniz(i,j,k)/
!     .(guniv * 4./3. * d(i,j,k) * 0.5 * (x1a(i+1) + x1a(i)) * dt), &
                         i=is,ie), j=js,js), k=ks,ks)
!      write(12,"(a)")
!      write(12,2) ( ( ( x1a(i)/cmpc,gp(i,j,k),
!     1                   i=is,ie), j=js,je), k=ks,ke)
!      write(12,"(a)")
!      write(12,"('         x1b             x2a             x3b',
!     1 '              v2')")
!      write(12,3) ( ( ( x1b(i)/cmpc,x2a(j)/cmpc,v2(i,j,k)/cmkm,
!     1                   i=is,ie), j=js,je), k=ks,ke)
!      write(12,"(a)")
!      write(12,"('         x1b             x2b             x3a',
!     1 '              v3')")
!      write(12,4) ( ( ( x1b(i)/cmpc, x2b(j), x3a(k), v3(i,j,k)/cmkm,
!     1                   i=is,ie), j=js,je), k=ks,ke)
      close(12)
      continue
      return
      end
