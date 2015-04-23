!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               H S M O C _ 1 D             \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                  University of California at San Diego
!
!=======================================================================
!
       subroutine hsmoc_1D( emf1, emf2, emf3 )
!
!    dac:zeus3d.mocemfs <-------------------------- MoC estimate of emfs
!                                                          october, 1992
!
!    written by: David Clarke
!    modified 1: Byung-Il Jun, July 1994
!                implemented John Hawley and Jim Stone's scheme to
!                fix pt. explosion of magnetic field in passive field.
!                Basically, this scheme mixes emfs computed with simple
!                upwinding(Evans and Hawley) and MoC.
!                The upwinded values are used to compute the wave
!                speeds for the characteristic cones for the MoC part.
!    modified 2: Robert Fiedler, February 1995
!                upgraded to ZEUS-3D version 3.4 -- improved cache
!                utilization and added parallelization directives for 
!                SGI multiprocessors.
!    modified 3: Mordecai-Mark Mac Low, December 1997 - March 1998
!                rewritten for ZEUS-MP without overlapping.  Calls to 
!                interpolation routines have been inlined.
!    modified 4: PSLi, December 1999
!                minor modications to prevent scratch arrays overwritten.
!
!    modified 5: JHayes, October 2005; created hsmoc_1D clone which
!                assumes symmetry about the J and K axes; after ZEUS3D.
!
!  PURPOSE:  Uses the Method of Characteristics (MoC, invented by Jim
!  Stone, John Hawley, Chuck Evans, and Michael Norman; see Stone and
!  Norman, ApJS, v80, p791) to evaluate the velocity and magnetic field
!  needed to estimate emfs that are properly upwinded in the character-
!  istic velocities for the set of equations describing transverse
!  Alfven waves.  This is *not* the full MHD characteristic problem, but
!  a subset which has been found (reference above) to yield good results
!  for the propagation of Alfven waves.  This routine differs from the
!  previous routines MOC1, MOC2, and MOC3 in version 3.1 in that the
!  Lorentz forces are computed *before* the emfs are estimated.  Thus,
!  the emfs now use the velocities which have been updated with all the
!  source terms, including the transverse Lorenz forces.
!
!  The characteristic equations governing the propagation of Alfven
!  waves in the 1-direction are (see ZEUS3D notes "Method of Character-
!  istics"):
!
!  "plus" characteristic (C+):
!
!    ( db1/dt + (v2 - a2) * db1/dx2 + (v3 - a3) * db1/dx3 ) / sqrt(d)
!  + ( dv1/dt + (v2 - a2) * dv1/dx2 + (v3 - a3) * dv1/dx3 )  =  S    (1)
!
!  "minus" characteristic (C-):
!
!    ( db1/dt + (v2 + a2) * db1/dx2 + (v3 + a3) * db1/dx3 ) / sqrt(d)
!  - ( dv1/dt + (v2 + a2) * dv1/dx2 + (v3 + a3) * dv1/dx3 )  = -S    (2)
!
!  where   a2 = b2/sqrt(d) is the Alfven velocity in the 2-direction
!          a3 = b3/sqrt(d) is the Alfven velocity in the 3-direction
!          g1, g2, g3 are the metric factors
!          S = b1 * ( b2/g2 * dg1/dx2 + b3/g3 * dg1/dx3 )
!
!  Equations (1) and (2) can be written in Lagrangian form:
!
!      1    D+/-           D+/-
!   ------- ----(b1)  +/-  ----(v1)  =  +/- S                        (3)
!   sqrt(d)  Dt             Dt
!
!  where the Lagrangian derivatives are given by
!
!  D+/-     d                   d                   d
!  ----  =  --  +  (v2 -/+ a2) ---  +  (v3 -/+ a3) ---               (4)
!   Dt      dt                 dx2                 dx3
!
!  Differencing equations (3) [e.g. D+(b1) = b* - b+; D-(b1) = b* - b-],
!  and then solving for the advanced time values of b* and v*, one gets:
!                             _                                _
!           sqrt (d+ * d-)   |     b+          b-               |
!  b* =  ------------------- |  -------- + --------- + v+ - v-  |    (5)
!        sqrt(d+) + sqrt(d-) |_ sqrt(d+)    sqrt(d-)           _|
!
!                 1
!  v* =  ------------------- [ v+*sqrt(d+) + v-*sqrt(d-) + b+ - b- ]
!        sqrt(d+) + sqrt(d-)                                         (6)
!
!     + S Dt
!
!  where b+(-), and v+(-) are the upwinded values of the magnetic field
!  and velocity interpolated to the time-averaged bases of C+ (C-), and
!  d+(-) are estimates of the density along each characteristic path
!  during the time interval Dt.
!
!  Equations (1) and (2) would suggest that when estimating "emf2" for
!  example, that the interpolated values for "v1" and "b1" be upwinded
!  in both the 2- and 3- components of the characteristic velocity.  It
!  turns out that this is impractical numerically, and so only the
!  "partial" characteristics are tracked.  While estimating "emf2", "v1"
!  and "b1" are upwinded only in the 3-component of the characteristic
!  velocities.  Conversely, while estimating "emf3", "v1" and "b1" are
!  upwinded only in the 2-component of the characteristic velocities.
!  Since updating "b1" requires both "emf2" and "emf3", the evolution of
!  "b1" will ultimately depend upon the full characteristics.  This
!  amounts to nothing more than directionally splitting the full MoC
!  algorithm.  The effects of such a directionally split implementation
!  are not fully known.  What is known is:
!
!  1) A non-directionally split implementation of the MoC algorithm is
!     not possible without either relocating the emfs to the zone
!     corners or the magnetic field components to the zone centres.
!     The former has been tried (change deck mocemf) and was found to
!     generate intolerable diffusion of magnetic field.  In addition,
!     the algorithm is not unconditionally stable.  The latter has not
!     been tried, but is dismissed on the grounds that div(B) will be
!     determined by truncation errors rather than machine round-off
!     errors.
!
!  2) A directionally split algorithm that is not also operator split
!     (ie, performing the Lorentz updates of the velocities separately
!     from the MoC estimation of the emfs) does not allow stable Alfven
!     wave propagation in 2-D.  Operator splitting the MoC algorithm so
!     that the transverse Lorentz forces are computed *before* the
!     magnetic field update does allow Alfven waves to propagate stably
!     in multi-dimensions but appears to introduce more diffusion for
!     sub-Alfvenic flow.  On the other hand, super-Alfvenic flow does
!     *not* appear to be more diffusive in the operator split scheme.
!
!  INPUT VARIABLES:
!
!  OUTPUT VARIABLES:
!    emf1      emf along 1-edge computed using MoC estimates of v2, b2,
!              v3, and b3.
!    emf2      emf along 2-edge computed using MoC estimates of v3, b3,
!              v1, and b1.
!    emf3      emf along 3-edge computed using MoC estimates of v1, b1,
!              v2, and b2.
!
!  LOCAL VARIABLES:
!
!    1-D variables
!    bave      spatially averaged magnetic field at edge
!    srdp      sqrt(density) along the plus  characteristic (C+)
!    srdm      sqrt(density) along the minus characteristic (C-)
!    vchp      characteristic velocity along C+ (v - va)
!    vchm      characteristic velocity along C- (v + va)
!    vpch      velocity interpolated to the time-centred footpoint of C+
!    vmch      velocity interpolated to the time-centred footpoint of C-
!    bpch      B-field  interpolated to the time-centred footpoint of C+
!    bmch      B-field  interpolated to the time-centred footpoint of C-
!    vsnm1     MoC estimate of v[n-1] used to evaluate emf[n], n=1,2,3
!    bsnm1     MoC estimate of b[n-1] used to evaluate emf[n], n=1,2,3
!
!    2-D variables
!    v3intj   upwinded v3 in 2-direction
!    b3intj   upwinded b3 in 2-direction
!    etc..
!
!    3-D variables
!    srd[n]    sqrt of spatially averaged density at [n]-face n=1,2,3
!    vsnp1     MoC estimate of v[n+1] used to evaluate emf[n], n=1,2,3
!    bsnp1     MoC estimate of b[n+1] used to evaluate emf[n], n=1,2,3
!    vsnp1     is reused as the vsnp1*bsnm1 term in emf[n]
!    bsnp1     is reused as the vsnm1*bsnp1 term in emf[n]
!
!  EXTERNALS:
!    BVALEMF1, BVALEMF2, BVALEMF3
!
! I have inlined these M-MML 8.3.98
!    X1ZC1D  , X2ZC1D  , X3ZC1D
!    X1INT1D , X2INT1D , X3INT1D
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
      integer  :: kone, km1,kp1,jone,jm1,jp1  !asif/jch
!
      real(rl) :: absb, sgnp, sgnm, q1, q2, src, &
                  qv1, qv2, qb1, qb2, q3 ,dqm ,dqp ,xi ,fact, &
                  dv(ijkn), db(ijkn)
!
      real(rl) :: bave (ijkn), srdp (ijkn), srdm (ijkn), &
                  srdpi(ijkn), srdmi(ijkn), vchp (ijkn), &
                  vchm (ijkn), vtmp (ijkn), btmp (ijkn), &
                  vpch (ijkn), vmch (ijkn), bpch (ijkn), &
                  bmch (ijkn), vsnm1(ijkn), bsnm1(ijkn), &
                  vave (ijkn), aave (ijkn)
!
      real(rl) :: vfl  (ijkn), vt   (ijkn), bt   (ijkn), &
                  vint (ijkn), bint (ijkn)
!
      real(rl) :: v3intj(jn,kn), b3intj(jn,kn), &
                  v2intk(jn,kn), b2intk(jn,kn), &
                  v1intk(kn,in), b1intk(kn,in), &
                  v3inti(kn,in), b3inti(kn,in), &
                  v2inti(in,jn), b2inti(in,jn), &
                  v1intj(in,jn), b1intj(in,jn)
!
      real(rl) :: emf1 (in,jn,kn), emf2 (in,jn,kn), emf3(in,jn,kn)
      real(rl) :: vsnp1(in,jn,kn), bsnp1(in,jn,kn)
!
!       equivalence   
!     .               ( bave    , w1da     ), 
!     .               ( srdp    , w1db     ),
!     1               ( srdm    , w1dc     ), 
!     .               ( srdpi   , w1dd     ),
!     1               ( srdmi   , w1de     ), 
!     .               ( vchp    , w1df     ),
!     1               ( vchm    , w1dg     ), 
!     .               ( vpch    , w1dh     ),
!     1               ( vmch    , w1di     ), 
!     .               ( bpch    , w1dj     ),
!     1               ( bmch    , w1dk     ),
!     1               ( vtmp    , w1dl     ),
!     1               ( btmp    , w1dm     ),
!     1               ( vave    , w1dn     ),
!     1               ( aave    , w1do     ),
!     1               ( vsnm1   , w1dp     ),
!     1               ( bsnm1   , w1dq     ),
!     1               ( vt      , w1dq     ),
!     1               ( bt      , w1dr     ),
!     1               ( vint    , w1ds     ),
!     1               ( bint    , w1dt     ),
!     1               ( vfl     , w1du     )
!
! there are no 2D scratch arrays like in ZEUS-3D, but these can all still 
! be equivalenced to each other. M-MML 4 Mar 98
!       equivalence   ( v3intj, v2intk, v1intk, v3inti, v2inti
!     2               , v1intj                          )
!     2             , ( b3intj, b2intk, b1intk, b3inti, b2inti
!     2               , b1intj                          )
!
!      Careful!  "wa3d" through "wc3d" are equivalenced in CT.
!                wa3d -> emf1    wb3d -> emf2    wc3d -> emf3
!
!      The worker arrays "we3d" and "wf3d" should still contain "srd2"
!  and "srd3" from LORENTZ. The worker array  "wd3d" will contain
!  "srd1", but "wd3d" is needed for "bsnp1" and "term2".  Thus "srd1"
!  will be recomputed once "srd3" is no longer needed.
!
!RAF We have more array space with wg3d, so don't recompute srd1.
!RAF SGIMP does not like equivalenced variables in the same loop nest,
!RAF so eliminate term1 and term2.
!
!PS
!       equivalence   
!                  ( srd1              , w3di     )
!     1             , ( srd2              , w3dj     )
!     1             , ( srd3              , w3df     ),
!     1               ( bsnp1             , w3dg     )
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!---- 1.  emf1 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!    BEGINNING OF I LOOPS
!
!
!      By following the characteristics of the flow in the 3-direction,
!  determine values for "v2" and "b2" ("vsnp1" and "bsnp1") to be used
!  in evaluating "emf1".
!
!     Compute upwinded b3 and v3 in the 2-direction and use them to
!     compute the wave speeds for the characteristic cones for the MoC.
!
!PS   Initialize vsnp1
       do k=1,kn
          do j=1,jn
             do i=1,in
                vsnp1(i,j,k)=0.
             enddo
          enddo
       enddo
!
       i = is  
       j = js
       k = ks
!
       kone = 0
       jone = 0
       km1  = k-kone
       kp1  = k+kone
       jm1  = j-jone
       jp1  = j+jone
!
       do 9 i=is,ie
        vfl   (j) = 0.5 * (v2(i,j,k) + v2(i,j,km1)) - vg2(j)
        vt    (j) = v3(i,j,k) - vg3(k)
        bt    (j) = b3(i,j,k)
7      continue
!
!           call x2int1d ( bt, vfl, iordb3, istpb3, k, i
!     1                  , g2b, g2bi, bint )
!           call x2int1d ( vt, vfl, iords3, istps3, k, i
!     1                  , g2b, g2bi, vint )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
        dqm      = ( vt(j  ) - vt(jm1) ) * dx2bi(j  )
        dqp      = ( vt(jp1) - vt(j  ) ) * dx2bi(jp1)
        dv (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
        dqm     = ( bt(j  ) - bt(jm1) ) * dx2bi(j  )
        dqp     = ( bt(jp1) - bt(j  ) ) * dx2bi(jp1)
        db (j  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
        fact = dt * g2bi(i)
        qv1    = vt(jm1) + dx2a(jm1) * dv (jm1)
        qv2    = vt(j  ) - dx2a(j  ) * dv (j  )
        qb1    = bt(jm1) + dx2a(jm1) * db (jm1)
        qb2    = bt(j  ) - dx2a(j  ) * db (j  )
!
        xi     = vfl(j) * fact
        q3     = sign ( haf, xi )
        vint(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (jm1) ) &
               + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
        bint(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (jm1) ) &
               + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
!
        emf1(i,j,k) = vint(j)
        emf2(i,j,k) = bint(j)
9      continue
!
      do 6 i=is,ie
       vsnp1(i,j,k) = (v2(i,j,k) - vg2(j)) * emf2(i,j,k)
       bsnp1(i,j,k) = b2(i,j,k) * emf1(i,j,k)
6     continue
!
!-----------------------------------------------------------------------
!
!      By following the characteristics of the flow in the 2-direction,
!  determine values for "v3" and "b3" ("vsnm1" and "bsnm1") to be used
!  in evaluating "emf1".
!
       src = 0.0
!
!       Compute upwinded b2 and v2 in the 3-direction and use them to
!       compute the wave speeds for the chracteristic cones for the MoC
!
      do 59 i=is,ie
       vfl   (k) = 0.5 * (v3(i,j,k) + v3(i,jm1,k)) - vg3(k)
       vt    (k) = v2(i,j,k) - vg2(j)
       bt    (k) = b2(i,j,k)
!           call x3int1d ( bt, vfl, iordb2, istpb2, i, j
!     1                 , g31b, g31bi, g32a, g32ai, bint )
!           call x3int1d ( vt, vfl, iords2, istps2, i, j
!     1                 , g31b, g31bi, g32a, g32ai, vint )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
       dqm      = ( vt(k  ) - vt(km1) ) * dx3bi(k  )
       dqp      = ( vt(kp1) - vt(k  ) ) * dx3bi(kp1)
       dv (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
       dqm      = ( bt(k  ) - bt(km1) ) * dx3bi(k  )
       dqp      = ( bt(kp1) - bt(k  ) ) * dx3bi(kp1)
       db (k  ) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
       fact = dt * g31bi(i) * g32ai(j)
       qv1    = vt(km1) + dx3a(km1) * dv (km1)
       qv2    = vt(k  ) - dx3a(k  ) * dv (k  )
       qb1    = bt(km1) + dx3a(km1) * db (km1)
       qb2    = bt(k  ) - dx3a(k  ) * db (k  )
!  
       xi     = vfl(k) * fact
       q3     = sign ( haf, xi )
       vint(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) ) &
              + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
       bint(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) ) &
              + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
!
       emf1(i,j,k) = vint(k)
       emf2(i,j,k) = bint(k)
59    continue
!
      km1 = k - kone
      do 56 i=is,ie
       vsnp1(i,j,k) = 0.5 * (vsnp1(i,j,k) &
                    + b3(i,j,k)*emf1(i,j,k) )
       bsnp1(i,j,k) = 0.5 * (bsnp1(i,j,k) &
                    + (v3(i,j,k) - vg3(k))*emf2(i,j,k) )
56    continue
!
!     END OF I LOOP
!
100    continue
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif
       call bvalemf1 ( vsnp1, bsnp1 )
!
!  Wait for communications to complete.
!
#ifdef MPI_USED
       if(nreq .ne. 0) &
           call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!
!      Compute "emf1" for all 1-edges, including the ghost zones.
!
       do 110 i=is-2,ie+2
        if(xvgrid) then
         emf1(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                     * dx1ah(i)
        else
         emf1(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                     * dx1a (i)
        endif
110    continue
!
!-----------------------------------------------------------------------
!---- 2.  emf2 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!     BEGINNING OF FIRST J LOOP
!
!      By following the characteristics of the flow in the 1-direction,
!  determine values for "v3" and "b3" ("vsnp1" and "bsnp1") to be used
!  in evaluating "emf2".
!
      src = 0.0
!
!      Compute upwinded b1 and v1 in the 3-direction and use them to
!      compute the wave speeds for the chracteristic cones for the MoC.
!
      do 139 i=is,ie+1
       vfl   (k) = 0.5*(v3(i,j,k) + v3(i-1,j,k)) - vg3(k)
       vt    (k) = v1(i,j,k) - vg1(i)
       bt    (k) = b1(i,j,k)
137   continue
!
!
!           call x3int1d ( bt, vfl, iordb1, istpb1, i, j
!     1                  , g31b, g31bi, g32a, g32ai, bint )
!           call x3int1d ( vt, vfl, iordb1, istpb1, i, j
!     1                  , g31b, g31bi, g32a, g32ai, vint )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
!
      dqm      = ( vt(k  ) - vt(km1) ) * dx3bi(k  )
      dqp      = ( vt(kp1) - vt(k  ) ) * dx3bi(kp1)
      dv (k  ) = max ( dqm * dqp, zro ) &
               * sign ( one, dqm + dqp ) &
               / max ( abs ( dqm + dqp ), tiny )
      dqm      = ( bt(k  ) - bt(km1) ) * dx3bi(k  )
      dqp      = ( bt(kp1) - bt(k  ) ) * dx3bi(kp1)
      db (k  ) = max ( dqm * dqp, zro ) &
               * sign ( one, dqm + dqp ) &
               / max ( abs ( dqm + dqp ), tiny )
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
      fact = dt * g31bi(i) * g32ai(j)
!
      qv1    = vt(km1) + dx3a(km1) * dv (km1)
      qv2    = vt(k  ) - dx3a(k  ) * dv (k  )
      qb1    = bt(km1) + dx3a(km1) * db (km1)
      qb2    = bt(k  ) - dx3a(k  ) * db (k  )
!  
      xi     = vfl(k) * fact
      q3     = sign ( haf, xi )
      vint(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) ) &
             + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
      bint(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) ) &
             + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
!
      v1intk(k,i) = vint(k)
      b1intk(k,i) = bint(k)
139      continue
!
!      Select an effective density and determine the characteristic
!  velocity using upwinded values for each characteristic in the
!  1-direction.
!
!
      do 140 i=is,ie+1
       vave (i) = v1intk(k,i)
       bave (i) = b1intk(k,i)
!            vave (i) = 0.5 * ( v1(i,j,k) + v1(i,j,km1) ) - vg1(i)
!            bave (i) = 0.5 * ( b1(i,j,k) + b1(i,j,km1) )
       absb     = abs ( bave(i) )
       aave (i) = 0.5 * absb * ( srd3(i,j,k) + srd3(i-1,j,k) ) &
                             / ( srd3(i,j,k) * srd3(i-1,j,k) )
       sgnp     = sign ( haf, vave(i) - aave(i) )
       sgnm     = sign ( haf, vave(i) + aave(i) )
       srdp (i) = ( 0.5 + sgnp ) * srd3(i-1,j,k) &
                + ( 0.5 - sgnp ) * srd3(i  ,j,k)
       srdm (i) = ( 0.5 + sgnm ) * srd3(i-1,j,k) &
                + ( 0.5 - sgnm ) * srd3(i  ,j,k)
       srdpi(i) = 1.0 / srdp(i)
       srdmi(i) = 1.0 / srdm(i)
       vchp (i) = vave(i) - absb * srdpi(i)
       vchm (i) = vave(i) + absb * srdmi(i)
140   continue
!
!      Interpolate 1-D vectors of "v3" and "b3" in the 1-direction to
!  the footpoints of both characteristics.
!
      do 150 i=is-2,ie+2
       vtmp(i) = v3(i,j,k) - vg3(k)
       btmp(i) = b3(i,j,k)
150   continue
!           call x1zc1d ( vtmp, vchp, vchm, iords3, istps3, j, k
!     1                 , vpch, vmch )
!           call x1zc1d ( btmp, vchp, vchm, iordb3, istpb3, j, k
!     1                 , bpch, bmch )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
      do 950 i=is-1,ie+1
       dqm      = ( vtmp(i  ) - vtmp(i-1) ) * dx1bi(i  )
       dqp      = ( vtmp(i+1) - vtmp(i  ) ) * dx1bi(i+1)
       dv (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
!
       dqm      = ( btmp(i  ) - btmp(i-1) ) * dx1bi(i  )
       dqp      = ( btmp(i+1) - btmp(i  ) ) * dx1bi(i+1)
       db (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
 950  continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
      do 960 i=is,ie+1
       qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
       qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
       qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
       qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
!
       xi     = vchp(i) * dt
       q3     = sign ( haf, xi )
       vpch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
       bpch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
!
       xi     = vchm(i) * dt
       q3     = sign ( haf, xi )
       vmch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
       bmch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
 960  continue
!
!      Evaluate "vsnp1" and "bsnp1" by solving the characteristic
!  equations.  The source term is non-zero for RTP coordinates since
!  dg31/dx1 = 1.0.
!
      do 160 i=is,ie+1
       q2           = sign ( one, bave(i) )
       if(lgeom .eq. 3) then
        src          = dt * dg31bd1(i) * g31ai(i) * bave(i) &
                     * ( b3  (i,j,k)    + b3  (i-1,j,k)    ) &
                     / ( srd3(i,j,k)**2 + srd3(i-1,j,k)**2 )
       endif ! lgeom
       vsnp1(i,j,k) = ( vpch (i) * srdp (i) + vmch (i) * srdm (i) &
                       + q2 * ( bpch(i) - bmch(i) ) ) &
                      / ( srdp (i) + srdm (i) ) + src
       bsnp1(i,j,k) = ( bpch (i) * srdpi(i) + bmch (i) * srdmi(i) &
                       + q2 * ( vpch (i) - vmch (i) ) ) &
                      / ( srdpi(i) + srdmi(i) )
       vsnp1(i,j,k) = vsnp1(i,j,k) * bave(i)
       bsnp1(i,j,k) = bsnp1(i,j,k) * vave(i)
160   continue
!
170   continue
!
!        END OF FIRST J LOOP
!
180    continue
!
!-----------------------------------------------------------------------
!
!
!       BEGINNING OF SECOND J LOOP
!
!
!
!      By following the characteristics of the flow in the 3-direction,
!  determine values for "v1" and "b1" ("vsnm1" and "bsnm1") to be used
!  in evaluating "emf2".
!
!       Compute upwinded b3 and v3 in the 1-direction and use them to
!       compute the wave speeds for the chracteristic cones for the MoC
!
      do 188 i=is-2,ie+2
       vfl   (i) = 0.5*(v1(i,j,k) + v1(i,j,km1)) - vg1(i)
       vt    (i) = v3(i,j,k) - vg3(k)
       bt    (i) = b3(i,j,k)
188   continue
!
!           call x1int1d ( bt, vfl, iordb1, istpb1, j, k
!     1                  , bint )
!           call x1int1d ( vt, vfl, iords1, istps1, j, k
!     1                  , vint )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
      do 988 i=is-1,ie+1
       dqm      = ( vt(i  ) - vt(i-1) ) * dx1bi(i  )
       dqp      = ( vt(i+1) - vt(i  ) ) * dx1bi(i+1)
       dv (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
       dqm      = ( bt(i  ) - bt(i-1) ) * dx1bi(i  )
       dqp      = ( bt(i+1) - bt(i  ) ) * dx1bi(i+1)
       db (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
 988  continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
      do 987 i=is,ie+1
       qv1    = vt(i-1) + dx1a(i-1) * dv (i-1)
       qv2    = vt(i  ) - dx1a(i  ) * dv (i  )
       qb1    = bt(i-1) + dx1a(i-1) * db (i-1)
       qb2    = bt(i  ) - dx1a(i  ) * db (i  )
!
       xi     = vfl(i) * dt
       q3     = sign ( haf, xi )
       vint(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
       bint(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
 987  continue
!
      do 187 i=is,ie+1
       v3inti(k,i) = vint(i)
       b3inti(k,i) = bint(i)
187   continue
189   continue
!
       do 186 i = is,ie+1
        vsnp1(i,j,k) = 0.5*( vsnp1(i,j,k) &
                     + b1(i,j,k)*v3inti(k,i) )
        bsnp1(i,j,k) = 0.5*( bsnp1(i,j,k) &
                     + ( v1(i,j,k) - vg1(i) )*b3inti(k,i))
186    continue
!
!       END OF SECOND J LOOP
!
230    continue
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif
       call bvalemf2 ( vsnp1, bsnp1 )
!
!  Wait for communications to complete.
!
#ifdef MPI_USED
       if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!      Compute "emf2" for all 2-edges, including the ghost zones.
!
      do 240 i=is-2,ie+3
       if(xvgrid) then
        emf2(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                    * dx2ah(j) * g2ah(i)
       else
        emf2(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                    * dx2a (j) * g2a (i)
       endif
240   continue
!
!-----------------------------------------------------------------------
!---- 3.  emf3 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      BEGINNING OF K LOOP
!
!
!      By following the characteristics of the flow in the 2-direction,
!  determine values for "v1" and "b1" ("vsnp1" and "bsnp1") to be used
!  in evaluating "emf3".
!
!      Compute upwinded b2 and v2 in the 1-direction and use them to
!      compute the wave speeds for the chracteristic cones for MoC.
!
      do 267 i=is-2,ie+2
       vfl   (i) = 0.5 * (v1(i,j,k) + v1(i,jm1,k)) - vg1(i)
       vt    (i) = v2(i,j,k) - vg2(j)
       bt    (i) = b2(i,j,k)
267   continue
!           call x1int1d ( bt, vfl, iordb2, istpb2, j, k
!     1                  , bint )
!           call x1int1d ( vt, vfl, iords2, istps2, j, k
!     1                  , vint )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
      do 1067 i=is-1,ie+1
       dqm      = ( vt(i  ) - vt(i-1) ) * dx1bi(i  )
       dqp      = ( vt(i+1) - vt(i  ) ) * dx1bi(i+1)
       dv (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
       dqm      = ( bt(i  ) - bt(i-1) ) * dx1bi(i  )
       dqp      = ( bt(i+1) - bt(i  ) ) * dx1bi(i+1)
       db (i  ) = max ( dqm * dqp, zro ) &
                * sign ( one, dqm + dqp ) &
                / max ( abs ( dqm + dqp ), tiny )
 1067 continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
      do 1068 i=is,ie+1
       qv1    = vt(i-1) + dx1a(i-1) * dv (i-1)
       qv2    = vt(i  ) - dx1a(i  ) * dv (i  )
       qb1    = bt(i-1) + dx1a(i-1) * db (i-1)
       qb2    = bt(i  ) - dx1a(i  ) * db (i  )
!
       xi     = vfl(i) * dt
       q3     = sign ( haf, xi )
       vint(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
       bint(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
1068  continue
!
      do 268 i=is,ie+1
       v2inti(i,j) = vint(i)
       b2inti(i,j) = bint(i)
268   continue
269   continue
!
      do 266 i = is,ie+1
       vsnp1(i,j,k) = (v1(i,j,k) - vg1(i))*b2inti(i,j)
       bsnp1(i,j,k) = b1(i,j,k)*v2inti(i,j)
266   continue
!
!-----------------------------------------------------------------------
!
!      By following the characteristics of the flow in the 1-direction,
!  determine values for "v2" and "b2" ("vsnm1" and "bsnm1") to be used
!  in evaluating "emf3".
!
       src = 0.0
!
!      Compute upwinded b1 and v1 in the 2-direction and use them to
!      compute the wave speeds for the chracteristic cones for Moc
!
       do 319 i=is,ie+1
           vfl   (j) = 0.5 * (v2(i,j,k) + v2(i-1,j,k)) - vg2(j)
           vt    (j) = v1(i,j,k) - vg1(i)
           bt    (j) = b1(i,j,k)
!         call x2int1d ( bt, vfl, iordb1, istpb1, k, i
!     1                , g2b, g2bi, bint )
!         call x2int1d ( vt, vfl, iords1, istps1, k, i
!     1                , g2b, g2bi, vint )
!
           dqm      = ( vt(j  ) - vt(jm1) ) * dx2bi(j  )
           dqp      = ( vt(jp1) - vt(j  ) ) * dx2bi(jp1)
           dv (j  ) = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
           dqm     = ( bt(j  ) - bt(jm1) ) * dx2bi(j  )
           dqp     = ( bt(jp1) - bt(j  ) ) * dx2bi(jp1)
           db (j  ) = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
         fact = dt * g2bi(i)
           qv1    = vt(jm1) + dx2a(jm1) * dv (jm1)
           qv2    = vt(j  ) - dx2a(j  ) * dv (j  )
           qb1    = bt(jm1) + dx2a(jm1) * db (jm1)
           qb2    = bt(j  ) - dx2a(j  ) * db (j  )
!
           xi     = vfl(j) * fact
           q3     = sign ( haf, xi )
           vint(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (jm1) ) &
                  + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
           bint(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (jm1) ) &
                  + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
!
           b1intj(i,j) = bint(j)
           v1intj(i,j) = vint(j)
319    continue
!
!      Select an effective density and determine the characteristic
!  velocity using upwinded values for each characteristic in the
!  1-direction.
!
         do 320 i=is,ie+1
           vave (i) = v1intj(i,j)
           bave (i) = b1intj(i,j)
!          vave (i) = 0.5 * ( v1(i,j,k) + v1(i,jm1,k) ) - vg1(i)
!          bave (i) = 0.5 * ( b1(i,j,k) + b1(i,jm1,k) )
           absb     = abs ( bave(i) )
           aave (i) = 0.5 * absb * ( srd2(i,j,k) + srd2(i-1,j,k) ) &
                                 / ( srd2(i,j,k) * srd2(i-1,j,k) )
           sgnp     = sign ( haf, vave(i) - aave(i) )
           sgnm     = sign ( haf, vave(i) + aave(i) )
           srdp (i) = ( 0.5 + sgnp ) * srd2(i-1,j,k) &
                    + ( 0.5 - sgnp ) * srd2(i  ,j,k)
           srdm (i) = ( 0.5 + sgnm ) * srd2(i-1,j,k) &
                    + ( 0.5 - sgnm ) * srd2(i  ,j,k)
           srdpi(i) = 1.0 / srdp(i)
           srdmi(i) = 1.0 / srdm(i)
           vchp (i) = vave(i) - absb * srdpi(i)
           vchm (i) = vave(i) + absb * srdmi(i)
320      continue
!
!      Interpolate 1-D vectors of "v2" and "b2" in the 1-direction to
!  the footpoints of both characteristics.
!
         do 330 i=is-2,ie+2
           vtmp(i) = v2(i,j,k) - vg2(j)
           btmp(i) = b2(i,j,k)
330      continue
!         call x1zc1d ( vtmp, vchp, vchm, iords2, istps2, j, k
!     1               , vpch, vmch )
!         call x1zc1d ( btmp, vchp, vchm, iordb2, istpb2, j, k
!     1               , bpch, bmch )
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
         do 1130 i=is-1,ie+1
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
 1130    continue
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
         do 1140 i=is,ie+1
           qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
           qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
           qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
           qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
!          
           xi     = vchp(i) * dt
           q3     = sign ( haf, xi )
           vpch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                  + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
           bpch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                  + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
!
           xi     = vchm(i) * dt
           q3     = sign ( haf, xi )
           vmch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) ) &
                  + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
           bmch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) ) &
                  + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
1140    continue
!
!      Evaluate "vsnm1" and "bsnm1" by solving the characteristic
!  equations.  The source term is non-zero for RTP coordinates since
!  dg2/dx1 = 1.0.  Compute the two terms in "emf3".
!
         do 340 i=is,ie+1
           q2           = sign ( one, bave(i) )
           if(lgeom .eq. 3) then
            src          = dt * dg2bd1(i) * g2ai(i) * bave(i) &
                         * ( b2  (i,j,k)    + b2  (i-1,j,k)    ) &
                         / ( srd2(i,j,k)**2 + srd2(i-1,j,k)**2 )
           endif ! lgeom
           vsnm1(i    ) = ( vpch (i) * srdp (i) + vmch (i) * srdm (i) &
                          + q2 * ( bpch(i) - bmch(i) ) ) &
                        / ( srdp (i) + srdm (i) ) + src
           bsnm1(i    ) = ( bpch (i)* srdpi(i) + bmch (i)* srdmi(i) &
                          + q2 * ( vpch(i) - vmch(i) ) ) &
                        / ( srdpi(i) + srdmi(i) )
!
           vsnm1(i    ) = vsnm1(i) * bave(i)
           bsnm1(i    ) = bsnm1(i) * vave(i)
!
           vsnp1(i,j,k) = 0.5*(  vsnp1(i,j,k) + bsnm1(i) )
           bsnp1(i,j,k) = 0.5*(  vsnm1(i)     + bsnp1(i,j,k) )
340      continue
350    continue
!
!       END OF K LOOP
!
360    continue
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif
       call bvalemf3 ( vsnp1, bsnp1 )
!
!  Wait for communications to complete.
!
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!      Compute "emf3" for all 3-edges, including the ghost zones.
!
      do 370 i=is-2,ie+3
       if(xvgrid) then
        emf3(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                    * dx3ah(k) * g31ah(i) * g32ah(j)
       else
        emf3(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) ) &
                    * dx3a (k) * g31a (i) * g32a (j)
       endif
370   continue
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               H S M O C _ 1 D             \\\\\\\\\\
!
!=======================================================================
