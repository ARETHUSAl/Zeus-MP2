!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 M O M X 1                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine momx1 (ibeg,iend,jbeg,jend,kbeg,kend &
                        ,s1,s2,s3,mflx,atwid1,vtwid,sflx,dq)
!
!    dac:zeus3d.momx1 <--------------- transports momenta in 1-direction
!    from jms:zeus2d.momx1, mln:zeus04.momz                    may, 1990
!
!    written by: David Clarke
!    modified 1: November, 1992 by David Clarke; momenta are now updated
!                between and including i=is,ie, j=js,je, and k=ks,ke to
!                allow for proper treatment of periodic boundaries.
!    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
!                for ZEUS-MP
!
!  PURPOSE:  Transports the three components of the momentum density in
!  the 1-direction using the consistent transport algorithm, including
!  the effects of grid compression.  The transported fluxes are thus
!  given by the mass fluxes times the time centred area of the control
!  volume faces times the interpolated velocities.  Interpolations are
!  performed in-line.
!
!  INPUT VARIABLES:
!    mflx    mass flux in 1-direction (computed in TRANX1)
!    s1      momentum density in 1-direction
!    s2      momentum density in 2-direction
!    s3      momentum density in 3-direction
!
! BOUNDARY VALUES USED:
!
!    Macro defined  var   ii    oi    ij    oj    ik    ok
!    -------------  ---  ----  ----  ----  ----  ----  ----
!                  mflx  is-1  ie+1  js-1        ks-1
!                    u1  is-2  ie+2  js-1        ks-1
!                    u2  is-2  ie+2
!                    u3  is-2  ie+2
!
!  OUTPUT VARIABLES:
!    s1      momentum density in 1-direction updated in the 1-direction
!    s2      momentum density in 2-direction updated in the 1-direction
!    s3      momentum density in 3-direction updated in the 1-direction
!
!  LOCAL VARIABLES:
!    vel     velocity used for upwinding in interpolation routine
!    vtwid   interpolated velocity
!    sflx    momentum fluxes
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use scratch
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, jm1
!
      real(rl) :: dqm, dqp, xi, q1, vel
!
      real(rl) :: atwid1(ijkn), atwid2(ijkn), atwid3(ijkn)
      real(rl) :: vtwid (ijkn), sflx  (ijkn), &
                  dq    (ijkn)
!
      real(rl) :: mflx(in,jn,kn), s1(in,jn,kn), &
                  s2  (in,jn,kn), s3(in,jn,kn)
!
!---------------------------- TRANSPORT S1 -----------------------------
!
!      Compute time-centred area factors.
!
       do 10 i=ibeg-1,iend+1
        if(xvgrid) then
         atwid1(i) = 0.5 * g2bh(i) * g31bh(i)
         atwid2(i) = 0.5 * g2a (i) * g2ah (i) * g31ah(i)
         atwid3(i) = 0.5 * g31a(i) * g2ah (i) * g31ah(i)
        else
         atwid1(i) = 0.5 * g2b (i) * g31b (i)
         atwid2(i) = 0.5 * g2a (i) * g2a  (i) * g31a (i)
         atwid3(i) = 0.5 * g31a(i) * g2a  (i) * g31a (i)
        endif
10     continue
      do 590 k=kbeg,kend
       do 580 j=jbeg,jend
        if(ldimen .eq. 1) then
         jm1 = j
        else
         jm1 = j - 1
        endif
!
!      Interpolate "v1" at the zone centres.
!
!       call x1fc3d ( v1, vel1, js, ks, je, ke, iords1, vtwid1 )
!
!  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
!
           do 70 i=ibeg-1,iend+1
             dqm        = (v1 (i  ,j,k) - v1 (i-1,j,k)) * dx1ai(i-1)
             dqp        = (v1 (i+1,j,k) - v1 (i  ,j,k)) * dx1ai(i  )
             dq(i)      = max ( dqm * dqp, zro ) &
                        * sign ( one, dqm + dqp ) &
                        / max ( abs ( dqm + dqp ), tiny )
70         continue
!
!  2.  Choose time averaged, upwinded interface value.
!
           do 100 i=ibeg-1,iend
!
!      Construct an i-average of "v1-vg1" to be used for interpolation.
!
             vel        = 0.5 * ( v1(i  ,j,k) - vg1(i  ) &
                                + v1(i+1,j,k) - vg1(i+1) )
             xi        = vel         * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v1 (i  ,j,k) &
                       + ( dx1b(i  ) - xi ) * dq  (i      ) ) &
                       + ( 0.5 - q1 ) * ( v1 (i+1,j,k) &
                       - ( dx1b(i+1) + xi ) * dq  (i+1    ) )
!
!      Construct the 1-momentum flux at the zone centres and perform
!  1-momentum advection.  Note that the timestep "dt" is hidden in the
!  mass flux.
!
             sflx (i    ) = ( mflx (i,j,k) + mflx (i+1,j,k) ) &
                          * vtwid (i    ) * atwid1(i)
100        continue
           do 170 i=ibeg,iend
            if(xvgrid) then
             s1   (i,j,k) = ( s1(i,j,k) * dvl1b(i) &
                          - sflx (i    ) + sflx (i-1    ) ) &
                          * dvl1bni(i)
            else
             s1   (i,j,k) = ( s1(i,j,k) * dvl1b(i) &
                          - sflx (i    ) + sflx (i-1    ) ) &
                          * dvl1bi(i)
            endif
170        continue
!         if(ldimen .gt. 1) then
!
!---------------------------- TRANSPORT S2 -----------------------------
!
!      Interpolate "v2" at the 1-interfaces.
!
!       call x1zc3d ( v2, vel1, js, ks, je, ke, iords2, istps2
!     1             , vtwid1, p      )
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
           do 270 i=ibeg-1,iend+1
             dqm       = (v2 (i  ,j,k) - v2 (i-1,j,k)) * dx1bi(i  )
             dqp       = (v2 (i+1,j,k) - v2 (i  ,j,k)) * dx1bi(i+1)
             dq(i)     = max ( dqm * dqp, zro ) &
                       * sign ( one, dqm + dqp ) &
                       / max ( abs ( dqm + dqp ), tiny )
270        continue
!
!     2.  Choose time averaged, upwinded interface value.
!
           do 300 i=ibeg,iend+1
             vel       = 0.5 * ( v1(i,jm1,k) + v1(i,j,k) )
             xi        = ( vel         - vg1(i) ) * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v2 (i-1,j,k) &
                       + ( dx1a(i-1) - xi ) * dq   (i-1) ) &
                       + ( 0.5 - q1 ) * ( v2 (i  ,j,k) &
                       - ( dx1a(i  ) + xi ) * dq   (i  ) )
!
!      Construct the 2-momentum flux at the 1-interfaces and perform
!  2-momentum advection.  Note that the timestep "dt" is hidden in the
!  mass flux.
!
             sflx (i    ) = ( mflx (i,jm1,k) + mflx (i,j,k) ) &
                          * vtwid (i    ) * atwid2(i)
300        continue
           do 370 i=ibeg,iend
            if(xvgrid) then
             s2(i,j,k) = ( s2(i,j,k) * dvl1a(i) &
                         - sflx (i+1    ) + sflx (i    ) ) * dvl1ani(i)
            else
             s2(i,j,k) = ( s2(i,j,k) * dvl1a(i) &
                         - sflx (i+1    ) + sflx (i    ) ) * dvl1ai(i)
            endif
370        continue
!         endif ! ldimen
!
!JH         if(ldimen .eq. 3) then
!
!---------------------------- TRANSPORT S3 -----------------------------
!
!      Interpolate "v3" at the 1-interfaces.
!
!       call x1zc3d ( v3, vel1, js, ks, je, ke, iords3, istps3
!     1             , vtwid1, p      )
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
           do 470 i=ibeg-1,iend+1
             dqm       = (v3 (i  ,j,k) - v3 (i-1,j,k)) * dx1bi(i  )
             dqp       = (v3 (i+1,j,k) - v3 (i  ,j,k)) * dx1bi(i+1)
             dq(i)     = max ( dqm * dqp, zro ) &
                       * sign ( one, dqm + dqp ) &
                       / max ( abs ( dqm + dqp ), tiny )
470        continue
!
!     2.  Choose time averaged, upwinded interface value.
!
           do 500 i=ibeg,iend+1
!
!      Construct a k-average of "v1" to be used for interpolation.
!
            if(ldimen .eq. 3) then
             vel       = 0.5 * ( v1(i,j,k-1) + v1(i,j,k) )
            else
             vel       = v1(i,j,k)
            endif
             xi        = ( vel         - vg1(i) ) * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v3 (i-1,j,k) &
                       + ( dx1a(i-1) - xi ) * dq   (i-1) ) &
                       + ( 0.5 - q1 ) * ( v3 (i  ,j,k) &
                       - ( dx1a(i  ) + xi ) * dq   (i  ) )
!
!      Construct the 3-momentum flux at the 1-interfaces and perform
!  3-momentum advection.  Note that the timestep "dt" is hidden in the
!  mass flux.
!
            if(ldimen .eq. 3) then
             sflx (i    ) = ( mflx (i,j,k-1) + mflx (i,j,k) ) &
                          * vtwid (i    ) * g32b(j) * atwid3(i)
            else
             sflx (i    ) = 2.0*( mflx (i,j,k) ) &
                          * vtwid (i    ) * g32b(j) * atwid3(i)
            endif
500        continue
           do 570 i=ibeg,iend
            if(xvgrid) then
             s3(i,j,k) = ( s3(i,j,k) * dvl1a(i) &
                         - sflx (i+1    ) + sflx (i    ) ) * dvl1ani(i)
            else
             s3(i,j,k) = ( s3(i,j,k) * dvl1a(i) &
                         - sflx (i+1    ) + sflx (i    ) ) * dvl1ai(i)
            endif
570        continue
!         else ! ldimen
!          do i = ibeg, iend
!           s3(i,j,k) = 0.0D0
!          enddo
!         endif ! ldimen
580    continue
590   continue
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 M O M X 1                 \\\\\\\\\\
!
!=======================================================================
!
!
