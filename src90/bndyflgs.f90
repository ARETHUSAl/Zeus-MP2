!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              B N D Y F L G S              \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine bndyflgs
!
!    dac:zeus3d.bndyflgs <---------------- sets secondary boundary flags
!                                                          october, 1990
!
!    written by: David Clarke
!    modified 1: by RAF 2/13/96 for ZEUS-MP
!
!  PURPOSE:  This subroutine sets the secondary integer boundary flags
!  ("niib2", "niib3", "niib23", etc.) which are used in addition to the
!  primary flags ("niib", etc.) to account for the staggered grid.
!  Thus, the following chart indicates the integer flag which sets the
!  various variables at the various boundaries:
!
!      boundary            variable(s)         integer flag
!
!      inner (outer) i     d, e, v1, b1        niib   (noib  )
!                          v2, b2              niib2  (noib2 )
!                          v3, b3              niib3  (noib3 )
!                          emf1, emf2, emf3    niib23 (noib23)
!      inner (outer) j     d, e, v2, b2        nijb   (nojb  )
!                          v3, b3              nijb3  (nojb3 )
!                          v1, b1              nijb1  (nojb1 )
!                          emf1, emf2, emf3    nijb31 (nojb31)
!      inner (outer) k     d, e, v3, b3        nikb   (nokb  )
!                          v1, b1              nikb1  (nokb1 )
!                          v2, b2              nikb2  (nokb2 )
!                          emf1, emf2, emf3    nikb12 (nokb12)
!
!  Note that there is a "pecking order" for the boundary types.  If two
!  adjacent 1s have different values for the primary integer flag
!  (see discussion in BVALD), the secondary integer flag is set
!  according to this order, which is currently:
!
!      3, 5, 1,-1, 4, 2
!
!  Thus, if niib(j,k) = 1 while niib(j-1,k) = 4, then niib2(j,k) is set
!  to 1.  But if niib(j,k) = 1 while niib(j-1,k) = 5, then niib2(j,k) is
!  set to 5.  The pecking order is determined by the order in which the
!  "if tests" are made in the loops below.
!
!  LOCAL VARIABLES:
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use config
      use param
      use bndry
      use root
      use grid
!
      implicit none
!
      integer :: i, j, k, im1, jm1, km1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
!      if(.NOT. xmhd) return
!
       do 20 k=ks-2,ke+3
         km1 =  max ( ks-2, k - 1 )
         do 10 j=js-2,je+3
           jm1 =  max ( js-2, j - 1 )
!
!      Inner i boundary
!
           if ( (niib(j,k) .eq. 4) .or. (niib(jm1,k) .eq. 4) ) &
             niib2(j,k) = 4
           if ( (niib(j,k) .eq. 1) .or. (niib(jm1,k) .eq. 1) ) &
             niib2(j,k) = 1
           if ( (niib(j,k) .eq.-1) .or. (niib(jm1,k) .eq.-1) ) &
             niib2(j,k) =-1
           if ( (niib(j,k) .eq. 5) .or. (niib(jm1,k) .eq. 5) ) &
             niib2(j,k) = 5
           if ( (niib(j,k) .eq. 3) .or. (niib(jm1,k) .eq. 3) ) &
             niib2(j,k) = 3
!
           if ( (niib(j,k) .eq. 4) .or. (niib(j,km1) .eq. 4) ) &
             niib3(j,k) = 4
           if ( (niib(j,k) .eq. 1) .or. (niib(j,km1) .eq. 1) ) &
             niib3(j,k) = 1
           if ( (niib(j,k) .eq.-1) .or. (niib(j,km1) .eq.-1) ) &
             niib3(j,k) =-1
           if ( (niib(j,k) .eq. 5) .or. (niib(j,km1) .eq. 5) ) &
             niib3(j,k) = 5
           if ( (niib(j,k) .eq. 3) .or. (niib(j,km1) .eq. 3) ) &
             niib3(j,k) = 3
!
           if ( (niib2(j,k) .eq. 4) .or. (niib2(j,km1) .eq. 4) ) &
             niib23(j,k) = 4
           if ( (niib2(j,k) .eq. 1) .or. (niib2(j,km1) .eq. 1) ) &
             niib23(j,k) = 1
           if ( (niib2(j,k) .eq.-1) .or. (niib2(j,km1) .eq.-1) ) &
             niib23(j,k) =-1
           if ( (niib2(j,k) .eq. 5) .or. (niib2(j,km1) .eq. 5) ) &
             niib23(j,k) = 5
           if ( (niib2(j,k) .eq. 3) .or. (niib2(j,km1) .eq. 3) ) &
             niib23(j,k) = 3
!
!      Outer i boundary
!
           if ( (noib(j,k) .eq. 4) .or. (noib(jm1,k) .eq. 4) ) &
             noib2(j,k) = 4
           if ( (noib(j,k) .eq. 1) .or. (noib(jm1,k) .eq. 1) ) &
             noib2(j,k) = 1
           if ( (noib(j,k) .eq.-1) .or. (noib(jm1,k) .eq.-1) ) &
             noib2(j,k) =-1
           if ( (noib(j,k) .eq. 5) .or. (noib(jm1,k) .eq. 5) ) &
             noib2(j,k) = 5
           if ( (noib(j,k) .eq. 3) .or. (noib(jm1,k) .eq. 3) ) &
             noib2(j,k) = 3
!
           if ( (noib(j,k) .eq. 4) .or. (noib(j,km1) .eq. 4) ) &
             noib3(j,k) = 4
           if ( (noib(j,k) .eq. 1) .or. (noib(j,km1) .eq. 1) ) &
             noib3(j,k) = 1
           if ( (noib(j,k) .eq.-1) .or. (noib(j,km1) .eq.-1) ) &
             noib3(j,k) =-1
           if ( (noib(j,k) .eq. 5) .or. (noib(j,km1) .eq. 5) ) &
             noib3(j,k) = 5
           if ( (noib(j,k) .eq. 3) .or. (noib(j,km1) .eq. 3) ) &
             noib3(j,k) = 3
!
           if ( (noib2(j,k) .eq. 4) .or. (noib2(j,km1) .eq. 4) ) &
             noib23(j,k) = 4
           if ( (noib2(j,k) .eq. 1) .or. (noib2(j,km1) .eq. 1) ) &
             noib23(j,k) = 1
           if ( (noib2(j,k) .eq.-1) .or. (noib2(j,km1) .eq.-1) ) &
             noib23(j,k) =-1
           if ( (noib2(j,k) .eq. 5) .or. (noib2(j,km1) .eq. 5) ) &
             noib23(j,k) = 5
           if ( (noib2(j,k) .eq. 3) .or. (noib2(j,km1) .eq. 3) ) &
             noib23(j,k) = 3
!
10       continue
20     continue
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       do 40 k=ks-2,ke+3
         km1 =  max ( ks-2, k - 1 )
         do 30 i=is-2,ie+3
           im1 =  max ( is-2, i - 1 )
!
!      Inner j boundary
!
           if ( (nijb(i,k) .eq. 4) .or. (nijb(i,km1) .eq. 4) ) &
             nijb3(i,k) = 4
           if ( (nijb(i,k) .eq. 1) .or. (nijb(i,km1) .eq. 1) ) &
             nijb3(i,k) = 1
           if ( (nijb(i,k) .eq.-1) .or. (nijb(i,km1) .eq.-1) ) &
             nijb3(i,k) =-1
           if ( (nijb(i,k) .eq. 5) .or. (nijb(i,km1) .eq. 5) ) &
             nijb3(i,k) = 5
           if ( (nijb(i,k) .eq. 3) .or. (nijb(i,km1) .eq. 3) ) &
             nijb3(i,k) = 3
!
           if ( (nijb(i,k) .eq. 4) .or. (nijb(im1,k) .eq. 4) ) &
             nijb1(i,k) = 4
           if ( (nijb(i,k) .eq. 1) .or. (nijb(im1,k) .eq. 1) ) &
             nijb1(i,k) = 1
           if ( (nijb(i,k) .eq.-1) .or. (nijb(im1,k) .eq.-1) ) &
             nijb1(i,k) =-1
           if ( (nijb(i,k) .eq. 5) .or. (nijb(im1,k) .eq. 5) ) &
             nijb1(i,k) = 5
           if ( (nijb(i,k) .eq. 3) .or. (nijb(im1,k) .eq. 3) ) &
             nijb1(i,k) = 3
!
           if ( (nijb3(i,k) .eq. 4) .or. (nijb3(im1,k) .eq. 4) ) &
             nijb31(i,k) = 4
           if ( (nijb3(i,k) .eq. 1) .or. (nijb3(im1,k) .eq. 1) ) &
             nijb31(i,k) = 1
           if ( (nijb3(i,k) .eq.-1) .or. (nijb3(im1,k) .eq.-1) ) &
             nijb31(i,k) =-1
           if ( (nijb3(i,k) .eq. 5) .or. (nijb3(im1,k) .eq. 5) ) &
             nijb31(i,k) = 5
           if ( (nijb3(i,k) .eq. 3) .or. (nijb3(im1,k) .eq. 3) ) &
             nijb31(i,k) = 3
!
!      Outer j boundary
!
           if ( (nojb(i,k) .eq. 4) .or. (nojb(i,km1) .eq. 4) ) &
             nojb3(i,k) = 4
           if ( (nojb(i,k) .eq. 1) .or. (nojb(i,km1) .eq. 1) ) &
             nojb3(i,k) = 1
           if ( (nojb(i,k) .eq.-1) .or. (nojb(i,km1) .eq.-1) ) &
             nojb3(i,k) =-1
           if ( (nojb(i,k) .eq. 5) .or. (nojb(i,km1) .eq. 5) ) &
             nojb3(i,k) = 5
           if ( (nojb(i,k) .eq. 3) .or. (nojb(i,km1) .eq. 3) ) &
             nojb3(i,k) = 3
!
           if ( (nojb(i,k) .eq. 4) .or. (nojb(im1,k) .eq. 4) ) &
             nojb1(i,k) = 4
           if ( (nojb(i,k) .eq. 1) .or. (nojb(im1,k) .eq. 1) ) &
             nojb1(i,k) = 1
           if ( (nojb(i,k) .eq.-1) .or. (nojb(im1,k) .eq.-1) ) &
             nojb1(i,k) =-1
           if ( (nojb(i,k) .eq. 5) .or. (nojb(im1,k) .eq. 5) ) &
             nojb1(i,k) = 5
           if ( (nojb(i,k) .eq. 3) .or. (nojb(im1,k) .eq. 3) ) &
             nojb1(i,k) = 3
!
           if ( (nojb3(i,k) .eq. 4) .or. (nojb3(im1,k) .eq. 4) ) &
             nojb31(i,k) = 4
           if ( (nojb3(i,k) .eq. 1) .or. (nojb3(im1,k) .eq. 1) ) &
             nojb31(i,k) = 1
           if ( (nojb3(i,k) .eq.-1) .or. (nojb3(im1,k) .eq.-1) ) &
             nojb31(i,k) =-1
           if ( (nojb3(i,k) .eq. 5) .or. (nojb3(im1,k) .eq. 5) ) &
             nojb31(i,k) = 5
           if ( (nojb3(i,k) .eq. 3) .or. (nojb3(im1,k) .eq. 3) ) &
             nojb31(i,k) = 3
!
30       continue
40     continue
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       do 60 j=js-2,je+3
         jm1 =  max ( js-2, j - 1 )
         do 50 i=is-2,ie+3
           im1 =  max ( is-2, i - 1 )
!
!      Inner k boundary
!
           if ( (nikb(i,j) .eq. 4) .or. (nikb(im1,j) .eq. 4) ) &
             nikb1(i,j) = 4
           if ( (nikb(i,j) .eq. 1) .or. (nikb(im1,j) .eq. 1) ) &
             nikb1(i,j) = 1
           if ( (nikb(i,j) .eq.-1) .or. (nikb(im1,j) .eq.-1) ) &
             nikb1(i,j) =-1
           if ( (nikb(i,j) .eq. 5) .or. (nikb(im1,j) .eq. 5) ) &
             nikb1(i,j) = 5
           if ( (nikb(i,j) .eq. 3) .or. (nikb(im1,j) .eq. 3) ) &
             nikb1(i,j) = 3
!
           if ( (nikb(i,j) .eq. 4) .or. (nikb(i,jm1) .eq. 4) ) &
             nikb2(i,j) = 4
           if ( (nikb(i,j) .eq. 1) .or. (nikb(i,jm1) .eq. 1) ) &
             nikb2(i,j) = 1
           if ( (nikb(i,j) .eq.-1) .or. (nikb(i,jm1) .eq.-1) ) &
             nikb2(i,j) =-1
           if ( (nikb(i,j) .eq. 5) .or. (nikb(i,jm1) .eq. 5) ) &
             nikb2(i,j) = 5
           if ( (nikb(i,j) .eq. 3) .or. (nikb(i,jm1) .eq. 3) ) &
             nikb2(i,j) = 3
!
           if ( (nikb1(i,j) .eq. 4) .or. (nikb1(i,jm1) .eq. 4) ) &
             nikb12(i,j) = 4
           if ( (nikb1(i,j) .eq. 1) .or. (nikb1(i,jm1) .eq. 1) ) &
             nikb12(i,j) = 1
           if ( (nikb1(i,j) .eq.-1) .or. (nikb1(i,jm1) .eq.-1) ) &
             nikb12(i,j) =-1
           if ( (nikb1(i,j) .eq. 5) .or. (nikb1(i,jm1) .eq. 5) ) &
             nikb12(i,j) = 5
           if ( (nikb1(i,j) .eq. 3) .or. (nikb1(i,jm1) .eq. 3) ) &
             nikb12(i,j) = 3
!
!      Outer k boundary
!
           if ( (nokb(i,j) .eq. 4) .or. (nokb(im1,j) .eq. 4) ) &
             nokb1(i,j) = 4
           if ( (nokb(i,j) .eq. 1) .or. (nokb(im1,j) .eq. 1) ) &
             nokb1(i,j) = 1
           if ( (nokb(i,j) .eq.-1) .or. (nokb(im1,j) .eq.-1) ) &
             nokb1(i,j) =-1
           if ( (nokb(i,j) .eq. 5) .or. (nokb(im1,j) .eq. 5) ) &
             nokb1(i,j) = 5
           if ( (nokb(i,j) .eq. 3) .or. (nokb(im1,j) .eq. 3) ) &
             nokb1(i,j) = 3
!
           if ( (nokb(i,j) .eq. 4) .or. (nokb(i,jm1) .eq. 4) ) &
             nokb2(i,j) = 4
           if ( (nokb(i,j) .eq. 1) .or. (nokb(i,jm1) .eq. 1) ) &
             nokb2(i,j) = 1
           if ( (nokb(i,j) .eq.-1) .or. (nokb(i,jm1) .eq.-1) ) &
             nokb2(i,j) =-1
           if ( (nokb(i,j) .eq. 5) .or. (nokb(i,jm1) .eq. 5) ) &
             nokb2(i,j) = 5
           if ( (nokb(i,j) .eq. 3) .or. (nokb(i,jm1) .eq. 3) ) &
             nokb2(i,j) = 3
!
           if ( (nokb1(i,j) .eq. 4) .or. (nokb1(i,jm1) .eq. 4) ) &
             nokb12(i,j) = 4
           if ( (nokb1(i,j) .eq. 1) .or. (nokb1(i,jm1) .eq. 1) ) &
             nokb12(i,j) = 1
           if ( (nokb1(i,j) .eq.-1) .or. (nokb1(i,jm1) .eq.-1) ) &
             nokb12(i,j) =-1
           if ( (nokb1(i,j) .eq. 5) .or. (nokb1(i,jm1) .eq. 5) ) &
             nokb12(i,j) = 5
           if ( (nokb1(i,j) .eq. 3) .or. (nokb1(i,jm1) .eq. 3) ) &
             nokb12(i,j) = 3
!
50       continue
60     continue
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////              B N D Y F L G S              \\\\\\\\\!
!=======================================================================
!
!
