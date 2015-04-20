c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////                 S I G M A                 \\\\\\\\\\
c
c=======================================================================
c
c
c
c    written by: 
c    modified 1: ?
c
c  PURPOSE: returns the variance sigma(k)
c
c  INPUT VARIABLES:   k    magnitude of the wave vector
c                     idx  spectral index
c
c  OUTPUT VARIABLES:
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
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
