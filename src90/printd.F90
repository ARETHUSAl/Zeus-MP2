!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
subroutine printd
!
!  FORMATTED WRITE OF SELECTED VARIABLES
!
!     written by: Jim Stone
!     date:       February,1993
!     modified1:  DAC from ZEUS-3D
!     modified2:  RAF 1/12/96 for ZEUS-MP
!     modified3:  RAF 3/12/97, for ablation test problem
!     modified4:  JCH 4/17/97, cons.h included even when RAD_EXP
!                 not defined
!     modified5:  M-MML 5/18/98, expanded to include the full ZEUS-3D 
!         suite, including max and min value, as well as things like 
!         rms velocity, and to run under MPI
!
!  PURPOSE:  Dumps scalar "history" variables in a formatted write
!  for analysis.  Currently implemented variables are:
!
!   scal( 1) = time
!   scal( 2) = dt
!   scal( 3) = mass
!   scal( 4) = total energy
!   scal( 5) = kinetic energy
!   scal( 6) = magnetic energy
!   scal( 7) = internal energy
!   scal( 8) = radiation energy
!   scal( 9) = gravitational potential energy
!   scal(10) = 0.5*d*v1**2
!   scal(11) = 0.5*d*v2**2
!   scal(12) = 0.5*d*v3**2
!   scal(13) = b1**2
!   scal(14) = b2**2
!   scal(15) = b3**2
!   scal(16) = b1 flux
!   scal(17) = b2 flux
!   scal(18) = b3 flux
!   scal(19) = 1 angular momentum
!   scal(20) = 2 angular momentum
!   scal(21) = 3 angular momentum
!   scal(22) = 1 center of mass
!   scal(23) = 2 center of mass
!   scal(24) = 3 center of mass
!
!
!   scal(25) =    dmin   - minimum density
!   scal(26) =    dmax   - maximum density
!   scal(27) =    emin   - minimum specific energy
!   scal(28) =    emax   - maximum specific energy
!   scal(29) =    pmin   - minimum pressure
!   scal(30) =    pmax   - maximum pressure
!   scal(31) =    v1min  - minimum velocity in 1-direction
!   scal(32) =    v1max  - maximum velocity in 1-direction
!   scal(33) =    v2min  - minimum velocity in 2-direction
!   scal(34) =    v2max  - maximum velocity in 2-direction
!   scal(35) =    v3min  - minimum velocity in 3-direction
!   scal(36) =    v3max  - maximum velocity in 3-direction
!   scal(37) =    vmin   - minimum speed (magnitude of velocity vector)
!   scal(38) =    vmax   - maximum speed (magnitude of velocity vector)
!   scal(39) =    b1min  - minimum mag. field in 1-direction
!   scal(40) =    b1max  - maximum mag. field in 1-direction
!   scal(41) =    b2min  - minimum mag. field in 2-direction
!   scal(42) =    b2max  - maximum mag. field in 2-direction
!   scal(43) =    b3min  - minimum mag. field in 3-direction
!   scal(44) =    b3max  - maximum mag. field in 3-direction
!   scal(45) =    bmin   - minimum mag. field strength (magnitude of B-vector)
!   scal(46) =    bmax   - maximum mag. field strength (magnitude of B-vector)
!   scal(47) =    dvbmin - minimum div(B) (normalised)
!   scal(48) =    dvbmax - maximum div(B) (normalised)
!   scal(49) =    dtcs   - sound speed time stop
!      dt**   - timestep limit caused by ** (eg., dtcs, dtqq, etc.)
!  LOCALS:
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use cons
      use radiation
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else 
      use mpino
#endif
      use mpipar
      implicit NONE
      integer  :: i,j,k
      integer  :: id_ctr
      integer, parameter :: nscal=24
      real(rl) ::  scal(nscal)
      real(rl) ::  mass,x1cm,x2cm,x3cm,eint,ek1,ek2,ek3,epot
      real(rl) ::  emag1,emag2,emag3,b1flx,b2flx,b3flx
      real(rl) ::  l1,l2,l3,erad,dvb
      real(rl) ::  dmin    , dmax    , semin   , semax   , pmin 
      real(rl) ::  pmax    , v1min   , v1max   , v2min   , v2max 
      real(rl) ::  v3min   , v3max   , vsqrmin , vsqrmax , vmin 
      real(rl) ::  vmax    , b1min   , b1max   , b2min   , b2max
      real(rl) ::  b3min   , b3max   , bsqrmin , bsqrmax , bmin 
      real(rl) ::  bmax    , dvbmin  , dvbmax
      integer  ::  imin(12), jmin(12), kmin(12)
      integer  ::  imax(12), jmax(12), kmax(12)
      real(rl) ::  rv,dvola,dm,g3b,q1,q2,q3,dar1,dar2,dar3 
      real(rl) ::  v1av,v2av,v3av,b1av,b2av,b3av
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================
!
!
!  Integrate quantities
!
!      Zero totals.
!
       id_ctr = 9999 
       mass  = tiny
       x1cm  = 0.0
       x2cm  = 0.0
       x3cm  = 0.0
       eint  = 0.0
       ek1   = 0.0
       ek2   = 0.0
       ek3   = 0.0
       epot  = 0.0
       emag1 = 0.0
       emag2 = 0.0
       emag3 = 0.0
       b1flx = 0.0
       b2flx = 0.0
       b3flx = 0.0
       l1    = 0.0
       l2    = 0.0
       l3    = 0.0
       erad  = 0.0
!
!      Compute totals over k-sweeps
!
       do 40 k=ks,ke
         do 30 j=js,je
           do 20 i=is,ie
             dvola     = dvl1a(i) * dvl2a(j) * dvl3a(k)
             dm        = d(i,j,k) * dvola
             g3b       = g31b(i) * g32b(j)
             q1        = x1b(i)
             q2        = x2b(j) * g2b(i)
             q3        = x3b(k) * g3b
             dar1      = g2a(i) * dx2b(j) * g3b    * dx3b(k)
             dar2      = g3b    * dx3b(k)          * dx1b(i)
             dar3      =          dx1b(i) * g2a(i) * dx2b(j)
             v1av      = v1(i,j,k) + v1(i+1,j  ,k  )
             v2av      = v2(i,j,k) + v2(i  ,j+1,k  )
             v3av      = v3(i,j,k) + v3(i  ,j  ,k+1)
             if (xmhd) then
               b1av      = b1(i,j,k) + b1(i+1,j  ,k  )
               b2av      = b2(i,j,k) + b2(i  ,j+1,k  )
               b3av      = b3(i,j,k) + b3(i  ,j  ,k+1)
             else
               b1av      = 0.0
               b2av      = 0.0
               b3av      = 0.0
             endif
             mass      = mass      + dm
             x1cm      = x1cm      + dm * q1
             x2cm      = x2cm      + dm * q2
             x3cm      = x3cm      + dm * q3
             ek1       = ek1       + dm * v1av**2
             ek2       = ek2       + dm * v2av**2
             ek3       = ek3       + dm * v3av**2
             if (.not. xiso) then
!JH
               if (.not. xtotnrg) then
                 eint      = eint      + e(i,j,k) * dvola
               else 
                 eint      = eint      + dvola * &
                         ( e(i,j,k) - 0.125 * d(i,j,k) * &
                                      (v1av**2 + v2av**2 + v3av**2) )
               endif 
             endif  
             if (lrad .ne. 0)      erad  = erad  + er(i,j,k) * dvola
             if (xgrav .or. xgrvfft) epot  = epot  - dm * gp(i,j,k)
             if (lgeom .eq. 1) then
               rv        = sqrt ( ( x1b(i) - x1ptm )**2 &
                                + ( x2b(j) - x2ptm )**2 &
                                + ( x3b(k) - x3ptm )**2 )
             else if (lgeom .eq. 2) then
              rv        = sqrt ( ( x1b(i) - x1ptm )**2 &
                               + ( x2b(j) - x2ptm )**2 &
                               + 2.0 * x2b(j) * x2ptm &
                               * ( 1.0 - cos(x3b(k) - x3ptm) ) )
             else 
               rv        = ( x1b(i)*cos(x2b(j)) - x1ptm*cos(x2ptm) )**2 &
                         + ( x1b(i)*sin(x2b(j)) - x1ptm*sin(x2ptm) )**2 &
                         + 2.0 * x1b(i)*sin(x2b(j)) * x1ptm*sin(x2ptm) &
                         * ( 1.0 - cos(x3b(k) - x3ptm) )
               rv        = sqrt ( rv )
             endif 
             if (xgrav .or. xgrvfft) epot = epot - dm * guniv * ptmass / rv
             if (xmhd) then
               emag1     = emag1     + b1av**2 * dvola
               emag2     = emag2     + b2av**2 * dvola
               emag3     = emag3     + b3av**2 * dvola
               b1flx     = b1flx     + b1(i,j,k) * dar1
               b2flx     = b2flx     + b2(i,j,k) * dar2
               b3flx     = b3flx     + b3(i,j,k) * dar3
             endif
             l1        = l1        + dm * ( q2 * v3av - q3 * v2av )
             l2        = l2        + dm * ( q3 * v1av - q1 * v3av )
             l3        = l3        + dm * ( q1 * v2av - q2 * v1av )
             if (lrad .ne. 0) then         
               if (  abs(x1a(i)-zro).lt.0.25*dx1a(i) &
               .and. abs(x2a(j)-zro).lt.0.25*dx2a(j) &
               .and. abs(x3a(k)-zro).lt.0.25*dx3a(k)) then
                 scal(19) = d (i,j,k)
                 scal(20) = e (i,j,k)
                 scal(21) = er(i,j,k)
                 id_ctr = myid
               endif
             endif ! lrad
!
20         continue ! i
30       continue ! j
40     continue ! k
!
       scal(1)  = time
       scal(2)  = dt
       scal(3)  = mass
!
       scal(7)  = eint
       scal(8)  = erad
       scal(9)  = epot
       scal(10) = ek1   / 8.0
       scal(11) = ek2   / 8.0
       scal(12) = ek3   / 8.0
       scal(13) = emag1 / 8.0
       scal(14) = emag2 / 8.0
       scal(15) = emag3 / 8.0
       scal(16) = b1flx
       scal(17) = b2flx
       scal(18) = b3flx
!
       if (lrad .ne. 0) then      
!......................................................................
!
! Special output for ablation problem:
!
        if(myid .eq. id_ctr) then
         scal(22) = gamm1 * scal(20) * 1.0e-12   !  pressure in Megabars
         scal(23) = scal(20) * (gamm1*mmw*mh) / (scal(19) * everg) ! T_gas
         scal(24) = sqrt(sqrt(scal(21)/rad_con)) * boltz/everg     ! T_rad
        endif
!
!......................................................................
!
       else                
!
         scal(19) = l1    / 2.0
         scal(20) = l2    / 2.0
         scal(21) = l3    / 2.0
!
         scal(22) = x1cm  / mass
         scal(23) = x2cm  / mass
         scal(24) = x3cm  / mass
!
       endif                  
!
! Total up various forms of energy.
!
       scal(5) = scal(10) + scal(11) + scal(12)             ! kinetic
       scal(6) = scal(13) + scal(14) + scal(15)             ! magnetic
       scal(4) = scal( 5) + scal( 6) + scal( 7) + scal( 8) + scal( 9)
!
!  Write out variables to file connected to unit 3 opened in MAIN
!  program unit (zeusmp.src)
!
      if (lrad .ne. 0 ) then           
        if (myid.eq.id_ctr) write(3,2001) (scal(i),i=1,nscal)
      else                    
        write(3,2001) (scal(i),i=1,nscal)
      endif                
2001  format(24(1p,6e13.5))
!2001  format(1p6e13.5/1p6e13.5/1p6e13.5/1p6e13.5/)
!
!
      return
end subroutine printd
