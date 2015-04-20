c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              T E X T D M P 2              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine textdmp2
c
c  written by:   Robert Fiedler
c  Last modified: 01-09-07 by Daniel Whalen for inclusion in ZEUS-MP
c                 2.1
c
c  PURPOSE:  Write out ASCII data dumps of hydrodynamical, chemistry,
c	     and radiation variables as the user specifies in a variety
c            of possible formats.
c
c
c-----------------------------------------------------------------------

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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE

      integer  :: i, j, k
      real(rl) :: mhi,cnorm,e_kinetic

      mhi = 1. / mh

c
c Assign unit 12 to the ASCII data file whose name is contained in the
c string variable 'usrfile'.
c

      open(12,file=usrfile)

 1    format(   1pe11.3e3    )
 2    format( 2(1pe11.3e3,1x))
 3    format( 3(1pe11.3e3,1x))
 4    format( 4(1pe11.3e3,1x))
c 4    format( 4(1pe16.8e3,1x))
 5    format( 5(1pe11.3e3,1x))
 6    format( 6(1pe11.3e3,1x))
 7    format( 7(1pe11.3e3,1x))
 8    format( 8(1pe13.6e2,1x))
 9    format( 9(1pe11.3e3,1x))
10    format(10(1pe11.3e3,1x))
11    format(11(1pe11.3e3,1x))
12    format(12(1pe16.8e3,1x))
c12    format(12(1pe11.3e3,1x))
13    format(13(1pe11.3e3,1x))
14    format(14(1pe11.3e3,1x))
 15   format(15(1pe16.8e3,1x))

c 
c Write out the data in the desired format
c
c--------------------------1D Sod_CMA outputs---------------------------
c
c      write(12,"('    x1b          X1          X2          X7        ', 
c     .'  X8          X9          X4          X5          X6       
c     .  HSum       HeSum       ChgSum')")
c      if      (ldimen .eq. 1) then

c      write(12,12)(((x1a(i),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
c     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
c     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
c     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
c     .   abun(i,j,k,9),
c     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
c     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
c     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
c     . , i=is,ie), j=js,je), k=ks,ke)

c      else if (ldimen .eq. 2) then

c      write(12,12)(((x2a(j),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
c     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
c     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
c     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
c     .   abun(i,j,k,9),
c     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
c     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
c     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
c     . , i=is,is), j=js,je), k=ks,ks)

c      else if (ldimen .eq. 3) then

c      write(12,12)(((x3a(k),abun(i,j,k,1), abun(i,j,k,2), abun(i,j,k,7),
c     .                      abun(i,j,k,8), abun(i,j,k,9), abun(i,j,k,4),
c     .                      abun(i,j,k,9), abun(i,j,k,4), fh - 
c     .   abun(i,j,k,1)-abun(i,j,k,2)-abun(i,j,k,7)-abun(i,j,k,8) - 
c     .   abun(i,j,k,9),
c     .   1 - fh       -abun(i,j,k,4)-abun(i,j,k,5)-abun(i,j,k,6) ,
c     .   abun(i,j,k,3)+abun(i,j,k,7)-abun(i,j,k,2)-abun(i,j,k,5) *
c     .   0.25         -abun(i,j,k,6)     * 0.5    -abun(i,j,k,9) * 0.5
c     . , i=is,is), j=js,js), k=ks,ke)

c      endif

c-----------------------Sedov-Taylor blast output-----------------------

      cnorm   = boltz*tevk*mhi
c
c      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,"('     x1b        dens       energ       metal        v1
c     .')")
c      write(12,"('     x1b        dens       energ')")
c      write(12,3) ((( x1b(i)/cmpc, 
c     .  d(i,j,k),
c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
c     .  d   (i,j,k)  *mhi,
c     .  e(i,j,k)
c      .  tgas(i,j,k) 
c     .  gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / 
c     .((abun(i,j,k,1)       + abun(i,j,k,2)       + abun(i,j,k,3)
c     .+ abun(i,j,k,4) * qrt + abun(i,j,k,5) * qrt + abun(i,j,k,6) * qrt
c     .)*d(i,j,k)), 
c     .  d(i,j,k) * abun(i,j,k,10),
c     .  v1(i,j,k)/cmkm
c     . 7.6d-1*(1.-abun(i,j,k,10))-abun(i,j,k,1)-abun(i,j,k,2),
c     . 2.4d-1*(1.-abun(i,j,k,10))-abun(i,j,k,4)-abun(i,j,k,5)
c     .           -abun(i,j,k,6),
c     .fh*(1-abun(i,j,k,10))- abun(i,j,k,1) - abun(i,j,k,2),
c     .(1-fh)*(1-abun(i,j,k,10)) - abun(i,j,k,4)-abun(i,j,k,5)
c     .-abun(i,j,k,6),
c     .abun(i,j,k,1),abun(i,j,k,2),abun(i,j,k,3),abun(i,j,k,4),
c     .abun(i,j,k,5),abun(i,j,k,6),abun(i,j,k,10)
c     ., i= is,ie+2), j= js,js), k= ks,ks)

c-----------------------6-species radshock output-----------------------
c

      e_kinetic = 0.

      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,"('     x1b       dens         energ       tgas       HI      
c     .',
c     .             '        HII        elec       He        He+        
c     .  He++')")
      write(12,13) ((( x1b(i)/cmpc, d(i,j,k), 
     .  e(i,j,k) ,tgas(i,j,k)
     ., abun(i,j,k,1) ! / (abun(i,j,k,1) + abun(i,j,k,2))  
     ., abun(i,j,k,2) ! / (abun(i,j,k,1) + abun(i,j,k,2))
     ., abun(i,j,k,3) ! *  d   (i,j,k  )
     ., abun(i,j,k,4) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
     ., abun(i,j,k,5) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
     ., abun(i,j,k,6) ! / (abun(i,j,k,4) + abun(i,j,k,5) + abun(i,j,k,6))
     ., abun(i,j,k,7) 
     ., abun(i,j,k,8) 
     ., abun(i,j,k,9)
c     .,(gp(i,j,k)-gp(i-1,j,k))*dx1bi(i)*dt   
c     ., gp(i,j,k)
     ., i= is,ie), j= js,js), k= ks,ks)

      write(12,"(a)")

      do k=ks,ke
       do j=js,je
        do i=is,ie
         e_kinetic = e_kinetic + haf * qrt * (v1(i,j,k)+v1(i+1,j,k))**2
     .             * 4. * pi * thd * d(i,j,k) * (x1a(i+1)**3-x1a(i)**3) 
        enddo
       enddo
      enddo
         
      write(12,9) ceHItot  ,ciHItot  ,ceHeItot ,ciHeItot ,ceHeIItot,
     .            ciHeIItot,bremtot  ,ICtot    ,e_kinetic


c------------------------1D advection test output------------------------
c
c      write(12,1) time
c      write(12,*) nhy,dt    
c      write(12,"('     x1b       dens       energ       tgas       HI      
c     .',
c     .             '        HII        elec       He        He+        
c     .  He++')")
c      write(12,10) ((( x1b(i), d(i,j,k), e(i,j,k) ,tgas(i,j,k)
c     ., abun(i,j,k,1) * d(i,j,k), abun(i,j,k,2) * d(i,j,k)
c     ., abun(i,j,k,3) * d(i,j,k), abun(i,j,k,4) * d(i,j,k) 
c     ., abun(i,j,k,5) * d(i,j,k), abun(i,j,k,6) * d(i,j,k) 
c     ., i= is,ie), j= js,js), k= ks,ks)

c----------------------------2D GSF outputs-----------------------------
c
c      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,*) iter
c      write(12,"('    x1b         dens       energ        tgas')")
c      write(12,4) (((x1b(i)/cmpc, d(i,j,k)/mh, e(i,j,k), tgas(i,j,k) 
c     ., i=is,ie), j=js,je), k=ks,ke)

c-------------------3-species hydro outputs with flux-------------------
c
c      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,*) iter
c      write(12,"('    x1b         dens        enrg        tgas       ', 
c     .'  HI         HII         elec         k24        flux')")
c      write(12,9) (((x1b(i)/cmpc 
c     .,(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, e(i,j,k)
c     .,d(i,j,k) * mhi, e(i,j,k)
c     ., tgas(i,j,k  ) 
c     ., abun(i,j,k,1) * d(i,j,k) * mhi, abun(i,j,k,2) * d(i,j,k) * mhi
c     ., abun(i,j,k,3) * d(i,j,k) * mhi, k24(i,j), f(i,j,1)
c     ., i=is,ie), j=js,je), k=ks,ke)

c---------------------3-species network diagnostics---------------------
c
c      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,*) iter
c      write(12,"('    x1b          k1         k2         k24       ', 
c     .' flux ')")
c      write(12,5) (((x1b(i)/cmpc 
c     ., k1(i,j) , k2(i,j) , k24(i,j), f(i,j,1)
c     ., i=is,ie), j=js,je), k=ks,ke)

c-----------------hydro outputs with n_strom (paper I)------------------
c
c      write(12,2) time/3.1536e07
c      write(12,*) nhy,dt    
c      write(12,*) iter
c      write(12,"('     x1b       dens       energ       tgas       HI      
c     .',
c     .             '        HII        elec       k24        n_strom')")
c      write(12,9) (((x1b(i)/cmpc 
c     .,(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, e(i,j,k) 
c     ., tgas(i,j,k  )
c     ., abun(i,j,k,1) * d(i,j,k) * mhi, abun(i,j,k,2) * d(i,j,k) * mhi
c     ., abun(i,j,k,3) * d(i,j,k) * mhi, k24(i,j)
c     ., 5.470D30 * x1b(i)**(-1.5)
c     ., i=is,ie), j=js,js), k=ks,ks)
c
c---------------------------timescale outputs---------------------------
c
c      write(12,1) time/3.1536e07
c      write(12,*) nhy,dt, iter    
c      write(12,*) itmcool, itmrate
c      write(12,"('     x1b        dens        tgas        elec       thy      
c     .dro       tpdv        tphi        trec      tlambda')")
c      write(12,9) (((x1b(i)/cmpc 
c      write(12,9) (((dfloat(i), dfloat(j) , dfloat(k), 
c     ., (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, tgas(i,j,k)
c     .,  abun(i,j,k,3) * d(i,j,k) * mhi
c     .,  thydro(i), tpdv(i), tphi(i), trec(i), tlambda(i)
c     .,i=is,ie), j=js,je), k=ks,ke)

c------------------Toronto T1 - T7 conference outputs-------------------

      if (usrtag .eq. 'T1_') then
        write(12,2) (((x1b(i)/cmkpc, abun(i,j,k,1),
     .                 i=is,ie), j=js,js), k=ks,ks)
        goto 68

      else if (usrtag .eq. 'T2_') then
        write(12,6) (((x1b(i)/cmkpc, abun(i,j,k,1), tgas(i,j,k),
     .                 i=is,ie), j=js,js), k=ks,ks)
        goto 68

      else if (usrtag .eq. 'T5_') then
        write(12,5) (((x1b(i)/cmkpc, (abun(i,j,k,1) + abun(i,j,k,2))
     .                 *mhi*d(i,j,k),
     .                 abun(i,j,k,1), gamm1 * e(i,j,k), tgas(i,j,k),
     .                 i=is,ie), j=js,js), k=ks,ks)
        goto 68

      else if (usrtag .eq. 'T6_') then
       write(12,4) (((x1b(i)/cmkpc, abun(i,j,k,1),
     .(abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi, tgas(i,j,k),
     .                 i=is,ie), j=js,je), k=ks,ke)
        goto 68

      else if (usrtag .eq. 'T7_') then
c        write(12,5) (((x1b(i)/cmkpc,
c     .                (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi,
c     .                 abun(i,j,k,1), gamm1 * e(i,j,k), tgas(i,j,k),
c     .                 i=is,ie), j=js,je), k=ks,ke)
        write(12,1) (((dsqrt(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2)/
     .                 dsqrt(gamma*gamm1*e(i,j,k)/d(i,j,k))
     .                , i=is,ie ) , j=js,je ), k=ks,ke )
        goto 68
      endif

c-----------------------3D I-front hydro outputs------------------------

c      write(12,*) time/3.1536e07, iter
c      write(12,"('     x1b         x2b         x3b        dens        
c     .energ       tgas         HI         HII         k24         
c     .flx')")
c      write(12,10)(((dfloat(i),dfloat(j),dfloat(k), 
c      write(12,10) (((x1b(i)/cmpc, x2b(j), x3b(k), 
c     .               (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k) * mhi,
c     .                e(i,j,k), tgas(i,j,k), abun(i,j,k,1)*d(i,j,k)*mhi, 
c     .                abun(i,j,k,2)*d(i,j,k)*mhi, k24(i,j), f(i,j,1),
c     .                i=is,ie), j=js,je), k=ks,ke)

c-----------------------2D halo collapse outputs------------------------

c      write(12,"('     x1b         x2b         x3b        dens        
c     .tgas')")
c      write(12,5) (((x1b(i)/cmpc, x2b(j)/cmpc, x3b(k), 
c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
c     .  d   (i,j,k)  *mhi , tgas(i,j,k),
c     .                i=is,ie), j=js,je), k=ks,ke)

c-----------------------6-species network outputs-----------------------

c      write(12,*) time/3.1536e07, iter
c      write(12,"('    x1b         dens        tgas         HI         
c     .HII         HeI         HeII        HeIII ')")
c      write(12,8) (((x1b(i)/cmpc, 
c     . (abun(i,j,k,1) + abun(i,j,k,2) + (abun(i,j,k,4)  + 
c     .  abun(i,j,k,5) + abun(i,j,k,6))* qrt) * d(i,j,k) * mhi,
c     .  tgas(i,j,k),
c     .  abun(i,j,k,1)*d(i,j,k)*mhi     , abun(i,j,k,2)*d(i,j,k)*mhi  
c     .  abun(i,j,k,4)*d(i,j,k)*mhi*qrt , abun(i,j,k,5)*d(i,j,k)*mhi*qrt, 
c     .  abun(i,j,k,6)*d(i,j,k)*mhi*qrt , 
c     .i=is,ie ) , j=js,js ), k=ks,ks )

c--------------------tgas and dens outputs for rjw99--------------------
c
c      write(12,15) (((tgas(i,j,k), d(i,j,k)/mh 
c     .                 ,i=is,ie ) , j=js,je ), k=ks,ke )
c
c-----------------------2D acc cutoff diagnostics-----------------------
c      mhi      = 1./mh
c      write(12,9) (((x1b(i),x2a(j),x3a(k),
c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
c     .  d   (i,j,k)  *mhi , abun(i,j,k,1), abun(i,j,k,2), 
c     .  v1(i,j,k)/cmkm, v2(i,j,k)/cmkm, v3(i,j,k)/cmkm
c     .                 ,i=is,ie ) , j=js,je ), k=ks,ke )
c-----------------------9-species network outputs-----------------------

      mhi      = 1./mh
      cnorm    = boltz*tevk*mhi
cc      write(12,*) time/3.1536e07
cc      write(12,"('    x1b         dens        tgas         HI       ', 
cc     .'  HII         HeI         HeII        HeIII        HM          H2I
cc     .         H2II')")

c      write(12,1) time/3.1536e07

c      write(12,15) (((x1b(i),
c     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
c     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
c     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
c     .  d   (i,j,k)  *mhi ,
c     .  d   (i,j,k),
c     .  e   (i,j,k)/d   (i,j,k) ,
c     .  v1  (i,j,k),
cc     .  gamm1 * dabs(e(i,j,k)),
c     .  tgas (i,j,k),
cc     .  gamm1 * dabs(e(i,j,k)) * (tevk/cnorm) / (1.66052992e-22 * (
cc     .  abun(i,j,k,1)       + abun(i,j,k,2)       + abun(i,j,k,3)
cc     .+ abun(i,j,k,4) * qrt + abun(i,j,k,5) * qrt + abun(i,j,k,6) * qrt
cc     .+ abun(i,j,k,7)       + abun(i,j,k,8) * haf + abun(i,j,k,9) * haf
cc     .                )),
c     . abun(i,j,k,1), ! /(abun(i,j,k,1)+abun(i,j,k,2)+2.*abun(i,j,k,8)),
cc     . (abun(i,j,k,2)     + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt)/
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf), 
cc     .  abun(i,j,k,1),!     /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,2),!/
c     .  abun(i,j,k,3),!/
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,4),!*qrt /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,5),!*qrt /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,6),!*qrt /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,7),!     /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
cc     .  gamm1 * e(i,j,k),
c     .  abun(i,j,k,8),!*haf /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
c     .  abun(i,j,k,9),!*haf /
cc     . (abun(i,j,k,1)     + abun(i,j,k,2)     + 
cc     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
cc     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf),
cc     .  gp(i,j,k), 
c     .  i=is,ie), j=js,js ), k=ks,ks )

c-----------------9-species network diagnostic outputs------------------

c      write(12,*) time/3.1536e07, iter
c      write(12,"('    x1b         k9         k11        k17         k29         
c     .k10         k18         k19         k28        k30')")
c      write(12,10) (((x1b(i)/cmpc, 
c     .  k9 (i,j)*abun(i,j,k,1)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi    ,
c     .  k11(i,j)*abun(i,j,k,8)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi*haf,  
c     .  k17(i,j)*abun(i,j,k,7)*abun(i,j,k,2)*d(i,j,k)**2*mhi*mhi    ,
c     .  k29(i,j)*abun(i,j,k,8)*              d(i,j,k)   *mhi    *haf,  
c     .  k10(i,j)*abun(i,j,k,9)*abun(i,j,k,1)*d(i,j,k)**2*mhi*mhi*haf,
c     .  k18(i,j)*abun(i,j,k,9)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi*haf,  
c     .  k19(i,j)*abun(i,j,k,9)*abun(i,j,k,7)*d(i,j,k)**2*mhi*mhi*haf,
c     .  k28(i,j)*abun(i,j,k,9)*              d(i,j,k)   *mhi    *haf,  
c     .  k30(i,j)*abun(i,j,k,9)*              d(i,j,k)   *mhi    *haf, 
c     .i=is,ie), j=js,js), k=ks,ks)

c----------------------3D halo reader test outputs----------------------

c      write(12,8) (((d(i,j,k), v1(i,j,k), v2(i,j,k), v3(i,j,k),  
c     .               e(i,j,k), x1b(i), x2b(j), x3b(k),
c     .               i=is,ie), j=js,je), k=ks,ke)

c--------------------multispecies cooling rate dumps--------------------

c      write(12,"(a)")
c      write(12,"('     x1b        ceco        cico       HH/DM    ',   
c     .      '    cmpt         H2         brem        recH')")
c      write(12,8) ((( x1b(i)/cmpc, 
c     .  ceHI (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
c     .  ciHI (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
c     .  DM   (i,j)*abun(i,j,k,1)*abun(i,j,k,1)*d(i,j,k)**2*mhi*mhi,
c     .  cmpt (i,j)*abun(i,j,k,3)              *d(i,j,k)       *mhi,
c     .  abun (i,j,k,8)*d(i,j,k) * gphdl(i,j)/(1.0+gphdl1/gpldl(i,j))
c     .               * haf * mhi,

c     .  brem  (i,j) * (
#ifdef H
c     .      abun(i,j,k,2)                       
#endif /* H */
#ifdef He
c     .    + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt 
#endif /* He */
c     .  ) * abun(i,j,k,3)*d(i,j,k)*d(i,j,k)*mhi*mhi,
                       
c     .  reHII(i,j)*abun(i,j,k,2)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
c     .i=is,ie), j=js,js), k=ks,ks)

c      write(12,"('     x1b        k1          k2        k24',  
c     .'       dedot   de/dedot')")
c      write(12,6) ((( x1a(i), 
c     .  k1   (i,j)*abun(i,j,k,1)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
c     .- k2   (i,j)*abun(i,j,k,2)*abun(i,j,k,3)*d(i,j,k)**2*mhi*mhi,
c     .  k24  (i,j)*abun(i,j,k,1)              *d(i,j,k)       *mhi,
c     .  dedot(i,j),abun(i,j,k,3)*d(i,j,k)*mhi/dedot(i,j),
c     .i=is,ie), j=js,js), k=ks,ks)

c--------------------------------velocities-----------------------------

 68   write(12,"(a)")
c      write(12,"('         x1a             x2b             x3b',
c     1 '              v1')")
      write(12,2) ( ( ( x1a(i)/cmpc, v1(i,j,k)/cmkm, !dv_ioniz(i,j,k)/
c     .(guniv * 4./3. * d(i,j,k) * 0.5 * (x1a(i+1) + x1a(i)) * dt),
     1                   i=is,ie), j=js,js), k=ks,ks)
c      write(12,"(a)")
c      write(12,2) ( ( ( x1a(i)/cmpc,gp(i,j,k),
c     1                   i=is,ie), j=js,je), k=ks,ke)
c      write(12,"(a)")
c      write(12,"('         x1b             x2a             x3b',
c     1 '              v2')")
c      write(12,3) ( ( ( x1b(i)/cmpc,x2a(j)/cmpc,v2(i,j,k)/cmkm,
c     1                   i=is,ie), j=js,je), k=ks,ke)
c      write(12,"(a)")
c      write(12,"('         x1b             x2b             x3a',
c     1 '              v3')")
c      write(12,4) ( ( ( x1b(i)/cmpc, x2b(j), x3a(k), v3(i,j,k)/cmkm,
c     1                   i=is,ie), j=js,je), k=ks,ke)

      close(12)

      continue

      return
      end
