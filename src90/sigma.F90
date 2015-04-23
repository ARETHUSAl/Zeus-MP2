!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
!    //////////                 S I G M A                 \\\\\\\\\\
!
!=======================================================================
!
!
!
!    written by: 
!    modified 1: ?
!
!  PURPOSE: returns the variance sigma(k)
!
!  INPUT VARIABLES:   k    magnitude of the wave vector
!                     idx  spectral index
!
!  OUTPUT VARIABLES:
!
!  LOCAL VARIABLES:
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
       function sigma(k, idx)
       use real_prec
       implicit none
       real(rl) :: sigma
       real(rl) :: k, idx
       real(rl) :: km
       km = k + 1d-99
       if ( idx.eq.0.0 ) then 
         sigma = 1.0
       else
         sigma = km**(idx)
       endif
       return
       end
