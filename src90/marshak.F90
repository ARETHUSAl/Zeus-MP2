!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              M A R S H A K                \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!     PURPOSE: initializes Marshak wave test problem.
!
!     Written by: Robert Fiedler and John Hayes
!
!=======================================================================
#define ABSORPTION        0.57735
      subroutine marshak
!
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k
      real(rl) :: d1, t0
      real(rl) :: third
!
      real(rl) :: d0     , e0     , er0    , soeps
      real(rl) :: ros_mfp, flx_lim, dmc_max, dx_min
!
      real(rl) :: so_eps , tau    , eriibn, timarg
!
      common /soui/ so_eps, tau, eriibn
!
      REAL souis
!
      namelist /pgen/ d0, e0, er0, soeps
!
!=======================================================================
!
!
!     initialize and read in parameters from PGEN namelist
!
       d0    =  1.0      !  mass density in gm/cc
       e0    =  0.0      !  gas energy
       er0   =  0.0      !  rad energy
       soeps =  0.1      !  epsilon defined in Su and Olson eq. (13)
!
      if (myid_w .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
       buf_in(1) = d0
       buf_in(2) = e0
       buf_in(3) = er0
       buf_in(4) = soeps
      endif
       call MPI_BCAST( buf_in, 4, MPI_FLOAT &
                     , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0    = buf_in(1)
       e0    = buf_in(2)
       er0   = buf_in(3)
       soeps = buf_in(4)
#endif /* MPI_USED */
      endif
      so_eps = soeps
!
! Compute the time-dependent BC (eriib) at time = 0.0.
!
       timarg = 0.0D0
       eriibn = souis (timarg)
!
! Copy the BC into the boundary value array.
!
       do k=ks-1,ke+1
         do j=js-1,je+1
           eriib(j,k,1) = eriibn
         enddo ! j
       enddo ! k
!
!     initialize field arrays
!
      do 1 k = 1, kn
      do 1 j = 1, jn
      do 1 i = 1, in
       d (i,j,k) = d0
       v1(i,j,k) = 0.0
       v2(i,j,k) = 0.0
       v3(i,j,k) = 0.0
       e (i,j,k) = e0
       er(i,j,k) = er0
1     continue
!
      return
      end
!
      subroutine source
!
      use real_prec
      use config
      use param
      use root
      use bndry
      use grid
!
      implicit NONE
!
      integer  :: i,j,k
!
      real(rl) :: souis, eriibn
!
      eriibn = souis(time)
      do k=ks-1,ke+1
        do j=js-1,je+1
          eriib(j,k,1) = eriibn
          eriib(j,k,2) = eriibn
        enddo ! j
      enddo ! k
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
!    //////////                 S O U I S                 \\\\\\\\\\
!
!=======================================================================
!
       real*8 function souis (time)
!
! PURPOSE
!   For non-equilibrium Marshak wave test problem (Su and Olson).
!   Returns the value of the radiation energy density (called u(x,tau))
!   at the surface x=0 at the current problem time for a given value
!   of the parameter so_eps (passed through common block soui).
!
!   Assumes that 4 * F_inc / c has been set to unity so that u(0,tau)
!   equals E(0,t).
!
!   This output is used to specify the time-dependent BC.
!
! AUTHOR
!   Robert A. Fiedler, 1/16/97
!
! LAST MODIFIED
!   1/16/97
!
! USES
!   qromo (Numerical Recipes), so_u_i1, so_u_i2 -- the integrands.
!
! TEST PROGRAM
!       program sou
!c
!c Test main program for routine souis.
!c
!       implicit NONE
!#include "param.h"
!#include "root.h"
!#include "cons.h"
!       integer i
!       REAL kap, so_eps, tau, eriibn
!       common /soui/ so_eps, tau, eriibn
!       REAL souis
!       external souis
!       REAL taus(11)
!       data taus /0.001, 0.003, 0.01, 0.03, 0.1, 0.3
!     1          ,1.0  , 3.0  , 10.0, 30.0, 100./
!       write(*,"('Enter epsilon')")
!       read (*,*) so_eps
!       kap = ABSORPTION
!       do i=1,11
!         time = taus(i) / epsilon / clight / kap
!         eriibn = souis(time)
!         write(*,"(' u(',1pe10.3,') = ',e15.5)") taus(i), eriibn
!       enddo
!       end
!.......................................................................
!
      use real_prec
      use config
      use param
!      use root
      use grid
      use cons
!
! Do not include bndry.h when using with test program sou.
!
      use bndry
!
       implicit NONE
!
       real(rl) :: so_eps, tau, eriibn
       common /soui/ so_eps, tau, eriibn
       real(rl) :: time
       real(rl) :: kap, qa, qb
       integer  :: j,k
!
       real(rl) :: so_u_i1, so_u_i2
!
       external so_u_i1, so_u_i2
!
       kap = ABSORPTION
       tau = so_eps * clight * kap * time
!
! Numerically evaluate the integrals of equation (36) in Su and Olson.
! They are functions of so_eps and tau(time).
!
       call qromo ( so_u_i1, zro, one, qa )
       call qromo ( so_u_i2, zro, one, qb )
!
! Find u(0,tau).
!
       souis = one - (sqrt(3.0)/pi) * ( two * qa + exp(-tau) * qb )
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                S O U I S                  \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
!    //////////               S O U I 1                   \\\\\\\\\\
!
!=======================================================================
!
       real*8 function so_u_i1 ( eta )
!
! First integrand in Su and Olson, eq. (36).
!
! Written by Robert Fiedler, 1/16/97.
!
      use real_prec
      use param
!
      implicit NONE
!
      real(rl) :: eta, g1, th
!
      real(rl) :: so_eps, tau, eriibn
      common /soui/ so_eps, tau, eriibn
!
      g1 = eta * sqrt ( so_eps + one / ( max ( one - eta**2, tiny ) ) )
      th = acos ( sqrt ( 3.0 / ( 3.0 + 4.0 * g1**2 ) ) )
      so_u_i1 = exp ( -tau * eta**2 ) * sin ( th ) &
               / ( eta * sqrt ( 3.0 + 4.0 * g1**2 ) )
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   F U N C T I O N            //////////
!    //////////                S O U I 1                  \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
!    //////////               S O U I 2                   \\\\\\\\\\
!
!=======================================================================
!
       real*8 function so_u_i2 ( eta )
!
! Second integrand in Su and Olson, eq. (36).
!
! Written by Robert Fiedler, 1/16/97.
!
      use real_prec
      use param
      implicit NONE
!
      real(rl) :: eta, g2, th
!
      real(rl) :: so_eps, tau, eriibn
      common /soui/ so_eps, tau, eriibn
!
       g2 = sqrt ( max ( one - eta, tiny ) * ( so_eps + one &
          / max ( eta, tiny ) ) )
       th = acos ( sqrt ( 3.0 / ( 3.0 + 4.0 * g2**2 ) ) )
       so_u_i2 = exp ( -tau / ( so_eps * max ( eta, tiny ) ) ) &
               * sin ( th ) &
               / ( max ( eta, tiny ) * ( one + so_eps * eta ) &
               * sqrt ( 3.0 + 4.0 * g2**2 ) )
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   F U N C T I O N            //////////
!    //////////                S O U I 2                  \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 Q R O M O                 \\\\\\\\\\
!
!=======================================================================
!
       subroutine qromo ( func, x1, x2, ss )
!
!    dac:zeus3d.qromb <------------------- integrates func from x1 to x2
!    from whp:numerical recipes                           december, 1992
!
!    written by: David Clarke
!    modified 1: Robert Fiedler 1/7/97, for improper integrals; 
!                see Numerical Recipes in Fortran, 2nd ed., p. 137.
!    modified 2: Robert Fiedler 1/16/97, quit for small errors dss;
!                tuned for Su and Olson non-equilibrium Marshak wave.
!
!  PURPOSE:  Returns the definite integral of the function "func"
!  between specified limits "x1" and "x2" using Romberg's method of
!  order 2k, where k=2 is Simpson's rule (see Numerical Recipes, 1st
!  edition for FORTRAN, page 114).
!
!  INPUT VARIABLES:
!    func     name of external function describing the univariate
!             function to be integrated.
!    x1, x2   integration limits
!
!  OUTPUT VARIABLES:
!    ss       value of definite integral
!
!  EXTERNALS:
!    FUNC
!    MIDPNT
!    POLINT
!
!-----------------------------------------------------------------------
      use real_prec
      use param
!
      implicit NONE
!
      integer, parameter :: isig=8
      integer, parameter :: jmx=20
      integer, parameter :: k=10
      integer, parameter :: km=k-1
!
      real(rl),parameter :: eps=0.1**isig
!
       integer  :: j
       real(rl) :: x1      , x2      , ss      , dss
!
       real(rl) :: s (jmx+1), h (jmx+1)
!
!      External statements
!
       real(rl) :: func
       external      func    , midpnt  , polint
!
!-----------------------------------------------------------------------
!
       h(1) = 1.0
       do 10 j=1,jmx
         call midpnt ( func, x1, x2, s(j), j )
         if (j .ge. k) then
           call polint ( h(j-km), s(j-km), k, zro, ss, dss )
!           if (abs(dss).lt.eps*abs(ss)) return
!
! Change from Numerical Recipes -- give up if the error is very small.
!
           if (abs(dss).lt.eps*abs(ss) .or. abs(dss).lt.tiny) return
         endif
         s(j+1) = s(j)
         h(j+1) = h(j) / 9.0  !  Step tripling and even error series.
10     continue
       write (6, 2000) jmx, dss, ss
       return
!
!-----------------------------------------------------------------------
!----------------------- Write format statements -----------------------
!-----------------------------------------------------------------------
!
2000   format('QROMO   : *** WARNING *** Romberg Integration failed ' &
             ,'to converge in ',i2,' steps.',/ &
             ,'QROMO   : error =',1pg12.5,', definite integral =',g12.5)
!
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 Q R O M O                 \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                M I D P N T                \\\\\\\\\\
!
!=======================================================================
!
       subroutine midpnt ( func, x1, x2, s, n )
!
!    dac:zeus3d.trapzd <-------------- n'th refinement of trapezoid rule
!    from whp:numerical recipes                           december, 1992
!
!    written by: David Clarke
!    modified 1: Robert Fiedler, 1/8/97, changed to midpnt for
!                improper integrals; see Numerical Recipes in Fortran,
!                2nd. ed., p. 136.
!
!  PURPOSE:  This routine returns in "s" the n'th stage of refinement of
!  an extended trapezoidal rule.  "func" is input as the name of the
!  function to be integrated between limits "x1" and "x2", also input.
!  When called with n=1, the crudest estimate of the definite integral
!  is returned.  Subsequent calls with n=2,3,... (in that order) will
!  improve the accuracy of "s" by adding 2**(n-2) additional interior
!  points.  The value of "s" should not be modified between successive
!  calls (see Numerical Recipes, 1st edition for FORTRAN, page 111).
!
!  INPUT VARIABLES:
!    func     name of external function describing the univariate
!             function to be integrated.
!    x1, x2   integration limits
!    n        stage of refinement
!
!  OUTPUT VARIABLES:
!    s        value of definite integral after n'th stage of refinement.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use param
!
      implicit NONE
!
      integer  :: n       , it      , j
      real(rl) :: x1      , x2      , s      , tnm     , del &
                   , x       , sum   , ddel
!
!      External statements
!
       real(rl) ::   func
       external      func
!
!-----------------------------------------------------------------------
!
       if (n .eq. 1) then
         s  = ( x2 - x1 ) * func ( haf * ( x1 + x2 ) )
       else
         it   = 3**( n - 2 )
         tnm  = real ( it )
         del  = ( x2 - x1 ) / ( 3.0 * tnm )
         ddel = del + del  !  Added pts alternate in spacing: del, ddel.
         x    = x1 + 0.5 * del
         sum  = 0.0
         do 10 j=1,it
           sum = sum + func ( x )
           x   = x + ddel
           sum = sum + func ( x )
           x   = x +  del
10       continue
         s  = ( s + ( x2 - x1 ) * sum / tnm ) / 3.0
       endif
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                M I D P N T                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                P O L I N T                \\\\\\\\\\
!
!=======================================================================
!
       subroutine polint ( xa, ya, n, x, y, dy )
!
!    dac:zeus3d.polint <---------------------- interpolates gridded data
!    from whp:numerical recipes                           december, 1992
!
!    written by: David Clarke
!    modified 1: Robert Fiedler, 1/7/97, for ZEUS-MP.
!
!  PURPOSE:  Given arrays "xa" and "ya", each of length "n", and a given
!  value "x", this routine returns a value "y" and an error estimate
!  "dy" (see Numerical Recipes, 1st edition for FORTRAN, page 82).
!
!  INPUT VARIABLES:
!    xa       input independent variable array
!    ya       input dependent   variable array
!    n        length of "xa" and "ya"
!    x        value of independent variable at which dependent variable
!             is to be interpolated
!
!  OUTPUT VARIABLES:
!    y        interpolated value of dependent variable
!    dy       estimate of error of interpolation.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use param
!
      implicit NONE
!
      integer, parameter :: nmax=20
!
      integer  :: n       , ns      , i       , m
      real(rl) :: x       , y       , dy      , dif     , dift &
                   , ho      , hp      , w       , den
!
      real(rl) :: xa      (   n), ya      (   n) &
                   , c       (nmax), d       (nmax)
!
!-----------------------------------------------------------------------
!
       ns  = 1
       dif = abs ( x - xa(1) )
       do 10 i=1,n
         dift = abs ( x - xa(i) )
         if (dift .lt. dif) then
           ns  = i
           dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
10     continue
       y  = ya(ns)
       ns = ns - 1
       do 30 m=1,n-1
         do 20 i=1,n-m
           ho   = xa(i  ) - x
           hp   = xa(i+m) - x
           w    = c (i+1) - d(i)
           den  = w / ( ( ho - hp ) + tiny )
           c(i) = ho * den
           d(i) = hp * den
20       continue
         if (2*ns .lt. n-m) then
           dy = c(ns+1)
         else
           dy = d(ns  )
           ns = ns - 1
         endif
         y = y + dy
30     continue
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                P O L I N T                \\\\\\\\\\
!
!=======================================================================
!
