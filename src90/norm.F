c=======================================================================
c
c    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
c    //////////                  N O R M                  \\\\\\\\\\
c
c=======================================================================
c
       function norm(idum)
       use real_prec
c
c    dac:zeus3d.norm <------------- returns normally distributed deviate
c    from whp:numerical recipes                           december, 1992
c
c    written by: David Clarke
c    modified 1:
c
c  PURPOSE:  Returns a normally distrbuted deviate with zero mean and
c  unit variance (see Numerical Recipes, 1st edition for FORTRAN, page
c  203, function GASDEV).
c
c  INPUT VARIABLES:
c    idum     .ge. 0 => get next random number from sequence
c             .lt. 0 => initialise sequence depending on absolute value
c                       of "idum", and gets next random number
c
c  OUTPUT VARIABLES:
c    norm     deviate with normal distribution
c
c  EXTERNALS: RAN1
c
c-----------------------------------------------------------------------
c
c
       integer     ::   idum, inext
       real(rl)    ::   v1, v2, r, fact, next, norm          
c
c      Data statements
c
       data          inext    / 0 /
       save          inext, next
c
c      External statements
c
       real(rl)    ::    ran1
       external          ran1
c
c-----------------------------------------------------------------------
c
       if (inext .eq. 0) then
10       continue
         v1    = 2.0 * ran1(idum) - 1.0
         v2    = 2.0 * ran1(idum) - 1.0
         r     = v1**2 + v2**2
         if (r .ge. 1.0) go to 10
         fact  = sqrt ( -2.0 *  log(r) / r )
         next  = v1 * fact
         norm  = v2 * fact
         inext = 1
       else
         norm  = next
         inext = 0
       endif
c
       return
       end
c
c
c=======================================================================
c
c    \\\\\\\\\\          E N D   F U N C T I O N          //////////
c    //////////                  N O R M                  \\\\\\\\\\
c
c=======================================================================
c
