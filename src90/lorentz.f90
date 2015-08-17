!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               L O R E N T Z               \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine lorentz (ibeg,iend,jbeg,jend,kbeg,kend &
                         ,u1,u2,u3,w1,w2,w3)
!
!    dac:zeus3d.lorentz <----- MoC estimate of transverse lorentz forces
!                                                          october, 1992
!
!    written by: David Clarke
!    modified 1: Mordecai-Mark Mac Low, December 1997 - March 1998
!                rewritten for ZEUS-MP 
!    modified 2: John Hayes, years ago; Rewritten for f90
!    modified 3: John Hayes, 09/01/2006; corrected metric factor
!                in second loop for "st3"
!
!  PURPOSE:  Uses the Method of Characteristics (MoC) to compute the
!  transverse components of the Lorentz force and then accelerate the
!  velocities accordingly.  After a suggestion by John Hawley, these
!  source terms are computed separately from the MoC estimate of the
!  emfs (MOCEMFS) and before the magnetic field update.  This improves
!  stability of multi-dimensional Alfven waves.
!
!  MOCEMFS solves the induction equation without the use of operator
!  splitting.  Thus, the relevant characteristic velocities are v+/-va
!  and includes the flow velocity v (from the advection term v grad(B))
!  and the alfven velocity va (from the induction term va grad(B)).
!  Conversely, Euler's equation is operator split, and the advection
!  term (v grad(v)) is handled elsewhere (MOMX*).  Thus, the effective
!  characteristic velocities when estimating the transverse Lorentz
!  terms is simply +/-va (from the induction term va grad(v)).
!
!  See comments in MOCEMFS for further ideas regarding the Method of
!  Characteristics.
!
!  INPUT VARIABLES:
!
!    u[1,2,3] = old velocities
!
!  OUTPUT VARIABLES:
!
!    w[1,2,3] = velocities from ibeg to iend, jbeg to jend, kbeg to kend
!
!  LOCAL VARIABLES:
!    srd[n]    sqrt of spatially averaged density at [n]-face n=1,2,3
!    bave      spatially averaged density at edge
!    srdpi     1.0/srd[n] along the plus  alfven characteristic (A+)
!    srdmi     1.0/srd[n] along the minus alfven characteristic (A-)
!    valp      characteristic velocity along A+ (-va)
!    valm      characteristic velocity along A- (+va)
!    vpal      velocity interpolated to the time-centred footpoint of A+
!    vmal      velocity interpolated to the time-centred footpoint of A-
!    bpal      B-field  interpolated to the time-centred footpoint of A+
!    bmal      B-field  interpolated to the time-centred footpoint of A-
!    bstar     MoC estimate of b[n] used to evaluate st[n], n=1,2,3
!    st[n]     Transverse Lorentz force on the [n]-velocity, n=1,2,3
!
!  EXTERNALS:
!    
!   -----> these are now inlined (M-MML): X1ZC1D  , X2ZC1D  , X3ZC1D
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use field
      use root
      use grid
      use scratch
      use lor_scr
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i       , j       , k      
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend
      integer  :: kone,km1,kp1  !asif
!
      real(rl) :: u1(in,jn,kn), u2(in,jn,kn), u3(in,jn,kn), &
                  w1(in,jn,kn), w2(in,jn,kn), w3(in,jn,kn)
!
      real(rl) :: absb ,q1, qv1, qv2, qb1, qb2, q3, dqm, dqp, &
                  xi, fact, &
                  dv(ijkn), db(ijkn)
!
      real(rl) :: bave(ijkn), srdpi(ijkn), srdmi(ijkn), &
                  valp(ijkn), valm (ijkn), vtmp (ijkn), &
                  btmp(ijkn), vpal (ijkn), vmal (ijkn), &
                  bpal(ijkn), bmal (ijkn), bstar(ijkn)
!
      real(rl) :: st1(in,jn,kn), st2(in,jn,kn), &
                  st3(in,jn,kn)
!
!       equivalence   
!     .               ( bave    , w1da     ),
!     .               ( srdpi   , w1db     ),
!     1               ( srdmi   , w1dc     ), ( vtmp    , w1dd     ),
!     2               ( vpal    , w1de     ), ( vmal    , w1df     )
!     2               ( valp    , bpal    , w1dg     ),
!     1               ( valm    , bmal    , w1dh     ),
!     8               ( btmp    , bstar   , w1di     ),
!     9               ( db      , w1dj     ),  ( dv     , w1dk     )
!
!      Routine MOCEMFS is expecting to find "srd1", "srd2", and "srd3" 
!  intact in worker arrays "wd3d", "we3d", and "wf3d" respectively.
!  [but MOCEMFS is in ZEUS-3D. does HSMOC also need these? M-MML 9.12.97 
!       so far, yes 4.3.98]
!
!PS: prevent overwritten of momenta.
!       equivalence   ( st1     , w3da     ), ( st2     , w3db     )
!     1             , ( st3     , w3dc     )
!PS -- end
!
!       equivalence   ( srd1    , w3di     ), 
!     .               ( srd2    , w3dj     ),
!     1               ( srd3    , w3df     )
!
!      External statements
!
!
!
!-----------------------------------------------------------------------
!
!      Compute face-centred averages of density, initialise source terms
!  to zero.
!
       if(ldimen.ge.3) then
       kone=1
       else
       kone=0
       endif   !asif: reduce 3D to 2D
       do 3 k=kbeg-kone,kend+kone
       km1=k-kone
         do 2 j=jbeg-1,jend+1
           do 1 i=ibeg-1,iend+1
             srd1(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i-1,j,k) ) )
             srd2(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i,j-1,k) ) )
             srd3(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i,j,km1) ) )
             st1 (i,j,k) = 0.0
             st2 (i,j,k) = 0.0
             st3 (i,j,k) = 0.0
1          continue
2        continue
3      continue
!
!-----------------------------------------------------------------------
!---- 1-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
!      By following the Alfven velocity in the 2-direction, evaluate
!  "bstar" from "b1" to estimate the 2-transverse Lorentz force.
!
       do 50 k=kbeg,kend
        km1=k-kone
         do 45 i=ibeg,iend+1
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 2-direction.
!
           do 10 j=jbeg,jend+1
             bave (j) = 0.5 * ( b2(i,j,k) + b2(i-1,j,k) )
             absb     = abs ( bave(j) )
             srdpi(j) = 1.0 / srd1(i,j  ,k)
             srdmi(j) = 1.0 / srd1(i,j-1,k)
             valp (j) =-absb * srdpi(j)
             valm (j) = absb * srdmi(j)
10         continue
!
!      Interpolate 1-D vectors of "u1" and "b1" in the 2-direction to
!  the footpoints of both Alfven characteristics.
!
           do 20 j=jbeg-2,jend+2
             vtmp(j) = u1(i,j,k) - vg1(i)
             btmp(j) = b1(i,j,k)
20         continue
!          call x2zc1d ( vtmp, valp, valm, iords1, istps1, k, i
!    1                 , g2a, g2ai, vpal, vmal )
!          call x2zc1d ( btmp, valp, valm, iordb1, istpb1, k, i
!    1                 , g2a, g2ai, bpal, bmal )
!
!  Inline interpolation (following tranx*)
!
           do 22 j=jbeg-1,jend+1
             dqm      = ( vtmp(j  ) - vtmp(j-1) ) * dx2bi(j  )
             dqp      = ( vtmp(j+1) - vtmp(j  ) ) * dx2bi(j+1)
             dv (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm     = ( btmp(j  ) - btmp(j-1) ) * dx2bi(j  )
             dqp     = ( btmp(j+1) - btmp(j  ) ) * dx2bi(j+1)
             db (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
22         continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           fact = dt * g2ai(i)
           do 24 j=jbeg,jend+1
             qv1    = vtmp(j-1) + dx2a(j-1) * dv (j-1)
             qv2    = vtmp(j  ) - dx2a(j  ) * dv (j  )
             qb1    = btmp(j-1) + dx2a(j-1) * db (j-1)
             qb2    = btmp(j  ) - dx2a(j  ) * db (j  )
!
             xi     = valp(j) * fact
             q3     = sign ( haf, xi )
             vpal(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bpal(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
!
             xi     = valm(j) * fact
             q3     = sign ( haf, xi )
             vmal(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bmal(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
 24       continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 30 j=jbeg,jend+1
             q1       = sign ( one, bave(j) )
             bstar(j) = ( bpal(j) * srdpi(j) + bmal(j) * srdmi(j) &
                        + q1 * ( vpal(j) - vmal(j) ) ) &
                      / ( srdpi(j) + srdmi(j) )
30         continue
!
!      Evaluate transverse Lorentz force.
!
           do 40 j=jbeg,jend
!JH
             if(lgeom .eq. 3) then
              if(x1a(i) .eq. 0.0) then
               st1(i,j,k) = 0.0D0
               go to 40
              endif
             endif
             st1(i,j,k) = ( bave (j+1) + bave (j) ) * g2ai(i) &
                        * ( bstar(j+1) - bstar(j) ) * dx2ai(j)
40         continue
45       continue
50     continue
!
!-----------------------------------------------------------------------
!	asif: for 3D only
        if (ldimen.eq.3) then
!
!      By following the Alfven velocity in the 3-direction, evaluate
!  "bstar" from "b1" to estimate the 3-transverse Lorentz force.
!
       do 100 j=jbeg,jend
         do 95 i=ibeg,iend+1
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 3-direction.
!
           do 60 k=kbeg,kend+kone
             km1=k-kone   !asif
             bave (k) = 0.5 * ( b3(i,j,k) + b3(i-1,j,k) )
             absb     = abs ( bave(k) )
             srdpi(k) = 1.0 / srd1(i,j,k  )
             srdmi(k) = 1.0 / srd1(i,j,km1)
             valp (k) =-absb * srdpi(k)
             valm (k) = absb * srdmi(k)
60         continue
!
!      Interpolate 1-D vectors of "u1" and "b1" in the 3-direction to
!  the footpoints of both Alfven characteristics.
!
           do 70 k=kbeg-2*kone,kend+2*kone    !asif added kone
             vtmp(k) = u1(i,j,k) - vg1(i)
             btmp(k) = b1(i,j,k)
70         continue
!           call x3zc1d ( vtmp, valp, valm, iords1, istps1, i, j
!     1                 , g31a, g31ai, g32b, g32bi, vpal, vmal )
!           call x3zc1d ( btmp, valp, valm, iordb1, istpb1, i, j
!     1                 , g31a, g31ai, g32b, g32bi, bpal, bmal )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
           do 72 k=kbeg-kone,kend+kone
            km1=k-kone
            kp1=k+kone
             dqm      = ( vtmp(k  ) - vtmp(km1) ) * dx3bi(k  )
             dqp      = ( vtmp(kp1) - vtmp(k  ) ) * dx3bi(kp1)
             dv (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm      = ( btmp(k  ) - btmp(km1) ) * dx3bi(k  )
             dqp      = ( btmp(kp1) - btmp(k  ) ) * dx3bi(kp1)
             db (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
72         continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           fact = dt * g31ai(i) * g32bi(j)
           do 74 k=kbeg,kend+kone
            km1=k-kone
            kp1=k+kone      !asif
             qv1    = vtmp(km1) + dx3a(km1) * dv (km1)
             qv2    = vtmp(k  ) - dx3a(k  ) * dv (k  )
             qb1    = btmp(km1) + dx3a(km1) * db (km1)
             qb2    = btmp(k  ) - dx3a(k  ) * db (k  )
!  
             xi     = valp(k) * fact
             q3     = sign ( haf, xi )
             vpal(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bpal(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
!
             xi     = valm(k) * fact
             q3     = sign ( haf, xi )
             vmal(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bmal(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
74         continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 80 k=kbeg,kend+kone
             q1       = sign ( one, bave(k) )
             bstar(k) = ( bpal(k) * srdpi(k) + bmal(k) * srdmi(k) &
                        + q1 * ( vpal(k) - vmal(k) ) ) &
                      / ( srdpi(k) + srdmi(k) )
80         continue
!
!      Evaluate transverse Lorentz force.
!
           q1 = g31ai(i) * g32bi(j)
           do 90 k=kbeg,kend
           kp1=k+kone
!JH
             if(lgeom .eq. 3) then
              if(x1a(i) .eq. 0.0) then
               st1(i,j,k) = 0.0D0
               go to 90
              endif
             endif
             st1(i,j,k) = st1(i,j,k) &
                        + ( bave (kp1) + bave (k) ) * q1 &
                        * ( bstar(kp1) - bstar(k) ) * dx3ai(k)
90         continue
95       continue
100      continue
       endif   !ldimen asif
!
110    continue
!-----------------------------------------------------------------------
!---- 2-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
!	asif: 3D only
        if (ldimen.eq.3) then
       do 165 j=jbeg,jend+1
!
!      By following the Alfven velocity in the 3-direction, evaluate
!  "bstar" from "b2" to estimate the 3-transverse Lorentz force.
!
         do 160 i=ibeg,iend
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 3-direction.
!
           do 120 k=kbeg,kend+kone
             km1=k-kone
             bave (k) = 0.5 * ( b3(i,j,k) + b3(i,j-1,k) )
             absb     = abs ( bave(k) )
             srdpi(k) = 1.0 / srd2(i,j,k  )
             srdmi(k) = 1.0 / srd2(i,j,km1)
             valp (k) =-absb * srdpi(k)
             valm (k) = absb * srdmi(k)
120        continue
!
!      Interpolate 1-D vectors of "u2" and "b2" in the 3-direction to
!  the footpoints of both Alfven characteristics.
!
           do 130 k=kbeg-2,kend+2
             vtmp(k) = u2(i,j,k) - vg2(j)
             btmp(k) = b2(i,j,k)
130        continue
!          call x3zc1d ( vtmp, valp, valm, iords2, istps2, i, j
!    1                 , g31b, g31bi, g32a, g32ai, vpal, vmal )
!          call x3zc1d ( btmp, valp, valm, iordb2, istpb2, i, j
!    1                 , g31b, g31bi, g32a, g32ai, bpal, bmal )
!
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
           do 132 k=kbeg-1,kend+1
             dqm      = ( vtmp(k  ) - vtmp(k-1) ) * dx3bi(k  )
             dqp      = ( vtmp(k+1) - vtmp(k  ) ) * dx3bi(k+1)
             dv (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm      = ( btmp(k  ) - btmp(k-1) ) * dx3bi(k  )
             dqp      = ( btmp(k+1) - btmp(k  ) ) * dx3bi(k+1)
             db (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
132         continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           fact = dt * g31bi(i) * g32ai(j)
           do 134 k=kbeg,kend+1
             qv1    = vtmp(k-1) + dx3a(k-1) * dv (k-1)
             qv2    = vtmp(k  ) - dx3a(k  ) * dv (k  )
             qb1    = btmp(k-1) + dx3a(k-1) * db (k-1)
             qb2    = btmp(k  ) - dx3a(k  ) * db (k  )
!  
             xi     = valp(k) * fact
             q3     = sign ( haf, xi )
             vpal(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (k-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bpal(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (k-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
!
             xi     = valm(k) * fact
             q3     = sign ( haf, xi )
             vmal(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (k-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bmal(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (k-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
134         continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 140 k=kbeg,kend+1
             q1       = sign ( one, bave(k) )
             bstar(k) = ( bpal(k) * srdpi(k) + bmal(k) * srdmi(k) &
                        + q1 * ( vpal(k) - vmal(k) ) ) &
                      / ( srdpi(k) + srdmi(k) )
140        continue
!
!      Evaluate transverse Lorentz force.
!
           q1 = g31bi(i) * g32ai(j)
           do 150 k=kbeg,kend
             st2(i,j,k) = ( bave (k+1) + bave (k) ) * q1 &
                        * ( bstar(k+1) - bstar(k) ) * dx3ai(k)
150        continue
160      continue
165    continue
!
         endif !ldimen  asif  7/31/03
!-----------------------------------------------------------------------
!
!      By following the Alfven velocity in the 1-direction, evaluate
!  "bstar" from "b2" to estimate the 1-transverse Lorentz force.
!
! *if def,IRIX.and.SGIMP
! C*$*ASSERT CONCURRENT CALL
! *endif IRIX.and.SGIMP
       do 210 k=kbeg,kend
         do 205 j=jbeg,jend+1
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 1-direction.
!
           do 170 i=ibeg,iend+1
             bave (i) = 0.5 * ( b1(i,j,k) + b1(i,j-1,k) )
             absb     = abs ( bave(i) )
             srdpi(i) = 1.0 / srd2(i  ,j,k)
             srdmi(i) = 1.0 / srd2(i-1,j,k)
             valp (i) =-absb * srdpi(i)
             valm (i) = absb * srdmi(i)
170        continue
!
!      Interpolate 1-D vectors of "u2" and "b2" in the 1-direction to
!  the footpoints of both Alfven characteristics.
!
           do 180 i=ibeg-2,iend+2
             vtmp(i) = u2(i,j,k) - vg2(j)
             btmp(i) = b2(i,j,k)
180        continue
!          call x1zc1d ( vtmp, valp, valm, iords2, istps2, j, k
!    1                 , vpal, vmal )
!          call x1zc1d ( btmp, valp, valm, iordb2, istpb2, j, k
!    1                 , bpal, bmal )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
           do 182 i=ibeg-1,iend+1
             dqm      = ( vtmp(i  ) - vtmp(i-1) ) * dx1bi(i  )
             dqp      = ( vtmp(i+1) - vtmp(i  ) ) * dx1bi(i+1)
             dv (i  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm      = ( btmp(i  ) - btmp(i-1) ) * dx1bi(i  )
             dqp      = ( btmp(i+1) - btmp(i  ) ) * dx1bi(i+1)
             db (i  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
 182       continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           do 184 i=ibeg,iend+1
             qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
             qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
             qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
!
             xi     = valp(i) * dt
             q3     = sign ( haf, xi )
             vpal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bpal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
!
             xi     = valm(i) * dt
             q3     = sign ( haf, xi )
             vmal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bmal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
184        continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 190 i=ibeg,iend+1
             q1       = sign ( one, bave(i) )
             bstar(i) = ( bpal(i) * srdpi(i) + bmal(i) * srdmi(i) &
                        + q1 * ( vpal(i) - vmal(i) ) ) * g2a(i) &
                      / ( srdpi(i) + srdmi(i) )
190        continue
!
!      Evaluate transverse Lorentz force.  Note that the metric
!  factor "g2a" has been absorbed into "bstar".
!
           do 200 i=ibeg,iend
             st2(i,j,k) = st2(i,j,k) &
                        + ( bave (i+1) + bave (i) ) * g2bi(i) &
                        * ( bstar(i+1) - bstar(i) ) * dx1ai(i)
200        continue
205      continue
210    continue
!
220    continue
!
!-----------------------------------------------------------------------
!---- 3-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
! *if def,IRIX.and.SGIMP
! C*$*ASSERT CONCURRENT CALL
! *endif IRIX.and.SGIMP
       do 330 k=kbeg,kend+kone   !asif 1->kone
        km1=k-kone
        kp1=k+kone
!
!      By following the Alfven velocity in the 1-direction, evaluate
!  "bstar" from "b3" to estimate the 1-transverse Lorentz force.
!
         do 270 j=jbeg,jend
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 1-direction.
!
           do 230 i=ibeg,iend+1
             bave (i) = 0.5 * ( b1(i,j,k) + b1(i,j,km1) )
             absb     = abs ( bave(i) )
             srdpi(i) = 1.0 / srd3(i  ,j,k)
             srdmi(i) = 1.0 / srd3(i-1,j,k)
             valp (i) =-absb * srdpi(i)
             valm (i) = absb * srdmi(i)
230        continue
!
!      Interpolate 1-D vectors of "u3" and "b3" in the 1-direction to
!  the footpoints of both Alfven characteristics.
!
           do 240 i=ibeg-2,iend+2
             vtmp(i) = u3(i,j,k) - vg3(k)
             btmp(i) = b3(i,j,k)
240        continue
!           call x1zc1d ( vtmp, valp, valm, iords3, istps3, j, k
!     1                 , vpal, vmal )
!           call x1zc1d ( btmp, valp, valm, iordb3, istpb3, j, k
!     1                 , bpal, bmal )
!
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
           do 242 i=ibeg-1,iend+1
             dqm      = ( vtmp(i  ) - vtmp(i-1) ) * dx1bi(i  )
             dqp      = ( vtmp(i+1) - vtmp(i  ) ) * dx1bi(i+1)
             dv (i  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm      = ( btmp(i  ) - btmp(i-1) ) * dx1bi(i  )
             dqp      = ( btmp(i+1) - btmp(i  ) ) * dx1bi(i+1)
             db (i  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
 242       continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           do 244 i=ibeg,iend+1
             qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
             qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
             qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
!
             xi     = valp(i) * dt
             q3     = sign ( haf, xi )
             vpal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bpal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
!
             xi     = valm(i) * dt
             q3     = sign ( haf, xi )
             vmal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bmal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
244        continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 250 i=ibeg,iend+1
             q1       = sign ( one, bave(i) )
             bstar(i) = ( bpal(i) * srdpi(i) + bmal(i) * srdmi(i) &
                        + q1 * ( vpal(i) - vmal(i) ) ) * g31a(i) &
                      / ( srdpi(i) + srdmi(i) )
250        continue
!
!      Evaluate transverse Lorentz force.  Note that the metric
!  factor "g31a" has been absorbed into "bstar".
!
           do 260 i=ibeg,iend
!JH         
            if(lgeom .eq. 3) then
             if(x1a(i) .eq. 0.0) then
              st3(i,j,k) = 0.0D0
              go to 260
             endif
            endif
             st3(i,j,k) = ( bave (i+1) + bave (i) ) * g31bi(i) &
                        * ( bstar(i+1) - bstar(i) ) * dx1ai(i)
260        continue
270      continue
!
!-----------------------------------------------------------------------
!
!      By following the Alfven velocity in the 2-direction, evaluate
!  "bstar" from "b3" to estimate the 2-transverse Lorentz force.
!
         do 320 i=ibeg,iend
!
!      Select an effective density and determine the Alfven speed for
!  each Alfven characteristic in the 2-direction.
!
           do 280 j=jbeg,jend+1
             bave (j) = 0.5 * ( b2(i,j,k) + b2(i,j,km1) )
             absb     = abs ( bave(j) )
             srdpi(j) = 1.0 / srd3(i,j  ,k)
             srdmi(j) = 1.0 / srd3(i,j-1,k)
             valp (j) =-absb * srdpi(j)
             valm (j) = absb * srdmi(j)
280        continue
!
!      Interpolate 1-D vectors of "u3" and "b3" in the 2-direction to
!  the footpoints of both Alfven characteristics.
!
           do 290 j=jbeg-2,jend+2
             vtmp(j) = u3(i,j,k) - vg3(k)
             btmp(j) = b3(i,j,k)
290        continue
!           call x2zc1d ( vtmp, valp, valm, iords3, istps3, k, i
!     1                 , g2b, g2bi, vpal, vmal )
!           call x2zc1d ( btmp, valp, valm, iordb3, istpb3, k, i
!     1                 , g2b, g2bi, bpal, bmal )
!
            do 292 j=jbeg-1,jend+1
             dqm      = ( vtmp(j  ) - vtmp(j-1) ) * dx2bi(j  )
             dqp      = ( vtmp(j+1) - vtmp(j  ) ) * dx2bi(j+1)
             dv (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
             dqm     = ( btmp(j  ) - btmp(j-1) ) * dx2bi(j  )
             dqp     = ( btmp(j+1) - btmp(j  ) ) * dx2bi(j+1)
             db (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
292         continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
           fact = dt * g2bi(i)
           do 294 j=jbeg,jend+1
             qv1    = vtmp(j-1) + dx2a(j-1) * dv (j-1)
             qv2    = vtmp(j  ) - dx2a(j  ) * dv (j  )
             qb1    = btmp(j-1) + dx2a(j-1) * db (j-1)
             qb2    = btmp(j  ) - dx2a(j  ) * db (j  )
!
             xi     = valp(j) * fact
             q3     = sign ( haf, xi )
             vpal(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bpal(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
!
             xi     = valm(j) * fact
             q3     = sign ( haf, xi )
             vmal(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) ) &
                    + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bmal(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) ) &
                    + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
 294       continue
!
!      Evaluate "bstar" by solving the characteristic equations.
!
           do 300 j=jbeg,jend+1
             q1       = sign ( one, bave(j) )
             bstar(j) = ( bpal(j) * srdpi(j) + bmal(j) * srdmi(j) &
                        + q1 * ( vpal(j) - vmal(j) ) ) * g32a(j) &
                      / ( srdpi(j) + srdmi(j) )
300        continue
!
!      Evaluate transverse Lorentz force.  Note that the metric
!  factor "g32a" has been absorbed into "bstar".
!
           do 310 j=jbeg,jend
!JH         
            if(lgeom .eq. 3) then
             if(x1a(i) .eq. 0.0) then
              st3(i,j,k) = 0.0D0
              go to 310
             endif
            endif
!JH
!JH         changed g2ai to g2bi since v3 lives at x1b 
!JH
            st3(i,j,k) = st3(i,j,k) &
                       + ( bave (j+1) + bave (j) ) * g2bi(i) * g32bi(j) &
                       * ( bstar(j+1) - bstar(j) ) * dx2ai(j)
310        continue
320      continue
!
330    continue
!
!-----------------------------------------------------------------------
!
!      Accelerate velocities and set boundaries.
!
       q1 = 0.5 * dt
       do 360 k=kbeg,kend
         do 350 j=jbeg,jend
           do 340 i=ibeg,iend
             w1(i,j,k) = u1(i,j,k) + q1 * st1(i,j,k) / srd1(i,j,k)**2
             w2(i,j,k) = u2(i,j,k) + q1 * st2(i,j,k) / srd2(i,j,k)**2
             w3(i,j,k) = u3(i,j,k) + q1 * st3(i,j,k) / srd3(i,j,k)**2
340        continue
350      continue
360    continue
!
!
!	asif: ZEUS-3d here calls bvalv*
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               L O R E N T Z               \\\\\\\\\!
!=======================================================================
!
