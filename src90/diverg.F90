!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                D I V E R G                \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine diverg ( c1, c2, c3, inorm, isum, div, sumd )
!
!    dac:zeus3d.diverg <------------ computes divergence of vector field
!    from jms:zeus2d.divb                                september, 1990
!
!    written by: David Clarke
!    modified 1: 11-17-99 by PSLi, for ZeusMP.
!
!  PURPOSE:  Computes divergence of vector field (c1, c2, c3), where
!  each of c1, c2, and c3 are face-centred quantities.  Therefore, the
!  divergence will be a zone-centred quantity, and is computed in the
!  range: ismn to iemx, jsmn to jemx, and ksmn to kemx.
!
!  INPUT VARIABLES:
!    c1       1-component of vector field
!    c2       2-component of vector field
!    c3       3-component of vector field
!    inorm    =1 => normalise divergence
!             =0 => do not normalise divergence
!    isum     =1 => perform reduction sum on "div"
!             =0 => no reduction
!
!  OUTPUT VARIABLES:
!    div      divergence of vector field
!    sumd     sum of "div"
!
!  LOCAL VARIABLES:
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
!
      use real_prec
      use config
      use param
      use grid
!
      implicit NONE
!
      integer  :: i, j, k, l, ip1, jp1, kp1, inorm, isum
      integer  :: kone, km1  !asif
!
      real(rl) :: sumd, sumn, nmax, nmaxi
!
      real(rl) :: norm(ijkn), sumdk(ijkn), sumnk(ijkn), nmaxk(ijkn)
!
      real(rl) :: c1 (in,jn,kn), c2(in,jn,kn), c3(in,jn,kn), &
                  div(in,jn,kn)
!
!      Careful!  "wa1d" - "wi1d" are used by BININT.
!
!       equivalence   ( norm    , ww1d     ), ( sumdk   , wx1d     )
!     1             , ( sumnk   , wy1d     ), ( nmaxk   , wz1d     )
!
!-----------------------------------------------------------------------
!	asif
        if (ldimen .eq. 3 )then
        kone=1
        else
        kone=0
        endif
!
!      Compute divergence of vector field ("div").
!
       do 30 k=ks-2*kone,ke+2*kone
         kp1 = k + kone
         do 20 j=js-2,je+2
           jp1 = j + 1
           do 10 i=is-2,ie+2
             ip1        = i + 1
             div(i,j,k) = ( g2a(ip1) * g31a(ip1) * c1(ip1,j,k) &
                          - g2a(i  ) * g31a(i  ) * c1(i  ,j,k) ) &
                        *                         dvl1ai(i) &
                        + ( g32a(jp1) * c2(i,jp1,k) &
                          - g32a(j  ) * c2(i,j  ,k) ) &
                        *   g2bi(i)             * dvl2ai(j) &
                        + ( c3(i,j,kp1) - c3(i,j,k) ) &
                        *   g31bi(i) * g32bi(j) * dvl3ai(k)
10         continue
20       continue
30     continue
!
!      Perform reduction on "div" if desired.  The reduction is done
!  first over each k-sweep, then for all k.  This is done to aid the
!  EDITOR autotasking process.
!
       sumd = 0.0
       if (isum .eq. 1) then
         do 60 k=ks,ke
           do 50 j=js,je
             do 40 i=is,ie
               sumd = sumd + div(i,j,k)
40           continue
50         continue
60       continue
       endif
!
!      Normalise divergence field if desired.
!
       if (inorm .eq. 1) then
!
!      Evaluate two normalising constants:
!
!  1.  sumn = sum over all grid zones the ratio: (average absolute
!      magnetic field) / (sum of grid zone dimensions), and
!
!  2.  nmax = maximum over the grid of the above ratios.
!
         sumn = 0.0
         nmax = 0.0
         do 90 k=ks,ke
           kp1      = k + kone
           do 80 j=js,je
             jp1 = j + 1
             do 70 i=is,ie
               ip1      = i + 1
               norm (1) = 0.5 * ( abs ( c1(ip1,j  ,k  ) + c1(i,j,k) ) &
                                + abs ( c2(i  ,jp1,k  ) + c2(i,j,k) ) &
                                + abs ( c3(i  ,j  ,kp1) + c3(i,j,k) ) ) &
                        / (dx1a(i) + g2b(i) * dx2a(j)   + &
                           g31b(i) * g32b(j) * dx3a(k) )
               sumn     = sumn     + norm(1)
               nmax     =   max ( nmax    , norm(1) )
70           continue
80         continue
90       continue
!
!      Apply normalising constant "sumn" to the scalar "sumd" and
!  "nmax" to the array "div".
!
         if (sumn .eq. 0.0) sumd  = 0.0
         if (sumn .ne. 0.0) sumd  = sumd / sumn
         if (nmax .eq. 0.0) nmaxi = 0.0
         if (nmax .ne. 0.0) nmaxi = 1.0 / nmax
         do 120 k=ks,ke
           do 110 j=js,je
             do 100 i=is,ie
               div(i,j,k) = div(i,j,k) * nmaxi
100          continue
110        continue
120      continue
!
       endif
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                D I V E R G                \\\\\\\\\\
!
!=======================================================================
!
