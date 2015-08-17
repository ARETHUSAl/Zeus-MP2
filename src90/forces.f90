!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                F O R C E S                \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine forces (ibeg,iend,jbeg,jend,kbeg,kend &
                         ,u1,u2,u3,w1,w2,w3)
!
!
! Computes Pressure, Gravity, and Pseudo-Rotational Forces.
!   and includes magnetic pressure forces (longitudinal Lorentz forces)
!
! Arrays u1 , u2 , u3  hold the old velocity values, while
!        w1 , w2 , w3  receive the updated values.
!
! BOUNDARY VALUES USED:
!
!  Macro defined  var   ii    oi    ij    oj    ik    ok
!  -------------  ---  ----  ----  ----  ----  ----  ----
!                  d   is-1        js-1        ks-1
!  ZRP             d   is-1        js-1        ks-1  ke+1
!  RTP             d   is-1        js-1  je+1  ks-1  ke+1
!                  e   is-1        js-1        ks-1
!  TOTAL_ENERGY    u1  is-1  ie+1  js-1        ks-1
!  RTP             u2  is-1              je+1
!  TOTAL_ENERGY    u2  is-1        js-1  je+1  ks-1
!  TOTAL_ENERGY    u3  is-1        js-1        ks-1  ke+1
!  ZRP             u3              js-1              ke+1
!  RTP             u3  is-1        js-1              ke+1
!  RAD             er  is-1        js-1        ks-1
!  MHD             b1  is-1  ie+1  js-1  je+1  ks-1  ke+1
!  MHD             b2  is-1  ie+1  js-1  je+1  ks-1  ke+1
!  MHD             b3  is-1  ie+1  js-1  je+1  ks-1  ke+1
! Written by RAF; modified 3/13/97 (JCH); 
!     modified 2 Mar 1998 to add MHD (M-MML)
!
! Rewritten for F90 by John Hayes, God knows when, exactly.
!
! Last modified: May 2006, by John Hayes: corrected k indexing error
! in expression for derdx3 in v2 update.
! Modified: 09/01/2006, by John Hayes; implemented Sean Matt's
! corrections to components of point-mass potential in RTP geometry
! Modified: 12/20/2006, by John Hayes; corrected typo in formula
! for corrections to v2 due to rotational pseudoforces.  This bug
! was inherited from ZEUS-3D!!
! 
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
      use gravmod
      use cons
      use opac
      use mpiyes
      use mpipar
      use chem
      implicit NONE
!
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend,km1,kp1
      integer  :: jp1, jm1
      integer  :: i,j,k
!
      real(rl) :: u1 (in,jn,kn),u2 (in,jn,kn),u3 (in,jn,kn), &
                  w1 (in,jn,kn),w2 (in,jn,kn),w3 (in,jn,kn)
!
      real(rl) :: poo(in),pmo(in),pom(in)
      real(rl) :: st,rhoi
!
      real(rl) :: s2oo(in),s2po(in)
!
      real(rl) :: s3oo(in),s3mo(in),s3op(in),s3mp(in)
!
      real(rl) :: r2,r2i,prpq,pzpq,prpt,pzpt
      real(rl) :: d1b2oo(in), d1b3oo(in), d2b1oo(in), d2b3oo(in), &
                  d3b1oo(in), d3b2oo(in)
      real(rl) :: d1b2po(in), d1b3op(in), d2b1po(in), d2b3op(in), &
                  d3b1po(in), d3b2op(in)
!
      real(rl) :: rdmcp, ros_mfp, flx_lim, big_r, derdx1, derdx2, &
                  derdx3, rdmcm, grad_er
!
!-----------------------------------------------------------------------
!
! Main loop for pressure, magnetic pressure, self-gravity, and
!     rotational pseudo-forces.  
!
! Active zones for all velocity components are (is:ie, js:je, ks:ke).
!
! In its most general configuration, this loop needs 1 layer of 
! boundary data for each of d, e, u1, u2, u3, and gp.
!
! If GRAV or GRAV_FFT is not defined, we don't need gp.
! If ISO  is     defined, we don't need e .
!
! For CARTESIAN grids and TOTAL_ENERGY not defined, we don't need u1,u2,u3
! boundary values.  All we do need in this case is the inner boundary 
! data for d (and e).
!
! Of course, TOTAL_ENERGY and ISO are mutually exclusive. 
!
! **> at the moment, TOTAL_ENERGY and MHD also can't be used together
!       (M-MML 4 Mar 98)
!
!
       do 40 k=kbeg,kend
        if(ldimen .eq. 3) then
         km1 = k-1
         kp1 = k+1
        else
         km1 = ks
         kp1 = ks
        endif
        do 30 j=jbeg,jend
         if(ldimen .eq. 1) then
          jm1 = js
          jp1 = js
         else
          jm1 = j-1
          jp1 = j+1
         endif
         do 10 i=ibeg-1,iend
!
! Compute thermal pressure
!
           poo(i) = p(i  ,j  ,k  )
           pmo(i) = p(i  ,jm1,k  )
           pom(i) = p(i  ,j  ,km1)
!.......................................................................
!   
          if(lgeom .ne. 1) then
!
!      Construct momentum densities from velocities.
!
           s3oo(i)         = u3(i  ,j  ,k  )*0.5*g31b(i  )*g32b(j  ) &
                           * ( d(i  ,j  ,km1) + d(i  ,j  ,k  ) )
           s3op(i)         = u3(i  ,j  ,kp1)*0.5*g31b(i  )*g32b(j  ) &
                           * ( d(i  ,j  ,k  ) + d(i  ,j  ,kp1) )
           s3mo(i)         = u3(i  ,jm1,k  )*0.5*g31b(i  )*g32b(jm1) &
                           * ( d(i  ,jm1,km1) + d(i  ,jm1,k  ) )
           s3mp(i)         = u3(i  ,jm1,kp1)*0.5*g31b(i  )*g32b(jm1) &
                           * ( d(i  ,jm1,k  ) + d(i  ,jm1,kp1) )
          endif ! lgeom 
          if(lgeom .eq. 3) then
           s2oo(i)         = u2(i  ,j  ,k  )*0.5*g2b (i  ) &
                           * ( d(i  ,jm1,k  ) + d(i  ,j  ,k  ) )
           s2po(i)         = u2(i  ,jp1,k  )*0.5*g2b (i  ) &
                           * ( d(i  ,j  ,k  ) + d(i  ,jp1,k  ) )
          endif ! lgeom
 10      continue
!.......................................................................
         if(xmhd) then
!
!        Compute differences in the squares of the magnetic field
!        components for the longitudinal Lorentz forces.
!
          do 11 i=ibeg-1,iend
           d1b2oo(i) = ( ( g2b (i  ) * b2(i  ,j  ,k  ) )**2 &
                       - ( g2b (i-1) * b2(i-1,j  ,k  ) )**2 ) &
                        * g2ai (i) * g2ai (i)
           d1b2po(i) = ( ( g2b (i  ) * b2(i  ,jp1,k  ) )**2 &
                       - ( g2b (i-1) * b2(i-1,jp1,k  ) )**2 ) &
                        * g2ai (i) * g2ai (i)
           d1b3oo(i) = ( ( g31b(i  ) * b3(i  ,j  ,k  ) )**2 &
                       - ( g31b(i-1) * b3(i-1,j  ,k  ) )**2 ) &
                        * g31ai(i) * g31ai(i)
           d1b3op(i) = ( ( g31b(i  ) * b3(i  ,j  ,kp1) )**2 &
                       - ( g31b(i-1) * b3(i-1,j  ,kp1) )**2 ) &
                        * g31ai(i) * g31ai(i)
           d2b3oo(i) = ( ( g32b(j  ) * b3(i  ,j  ,k  ) )**2 &
                       - ( g32b(jm1) * b3(i  ,jm1,k  ) )**2 ) &
                        * g32ai(j) * g32ai(j)
           d2b3op(i) = ( ( g32b(j  ) * b3(i  ,j  ,kp1) )**2 &
                       - ( g32b(jm1) * b3(i  ,jm1,kp1) )**2 ) &
                        * g32ai(j) * g32ai(j)
           d2b1oo(i) = b1(i  ,j  ,k  ) * b1(i  ,j  ,k  ) &
                     - b1(i  ,jm1,k  ) * b1(i  ,jm1,k  )
           d2b1po(i) = b1(i+1,j  ,k  ) * b1(i+1,j  ,k  ) &
                     - b1(i+1,jm1,k  ) * b1(i+1,jm1,k  )
           d3b1oo(i) = b1(i  ,j  ,k  ) * b1(i  ,j  ,k  ) &
                     - b1(i  ,j  ,km1) * b1(i  ,j  ,km1)
           d3b1po(i) = b1(i+1,j  ,k  ) * b1(i+1,j  ,k  ) &
                     - b1(i+1,j  ,km1) * b1(i+1,j  ,km1)
           d3b2oo(i) = b2(i  ,j  ,k  ) * b2(i  ,j  ,k  ) &
                     - b2(i  ,j  ,km1) * b2(i  ,j  ,km1)
           d3b2op(i) = b2(i  ,jp1,k  ) * b2(i  ,jp1,k  ) &
                     - b2(i  ,jp1,km1) * b2(i  ,jp1,km1)
11        continue
         endif ! xmhd
         do 20 i=ibeg,iend
!.......................................................................
!
! Perform an explicit update for v1
!
          rhoi    = 2.0 / ( d(i-1,j,k) + d(i,j,k) )
!
!  1.  pressure gradient
!
          st      = - rhoi &
                  * ( poo(i)    - poo(i-1)    ) * dx1bi(i)
!
!  2.  gravitational potential gradient
!
          if(xgrav .or. xgrvfft) then
           if (.not. xhse) then
            if(xsphgrv) then
             if(x1a(i) .ne. 0.0) then
              st      = st - guniv*intm(i)*x1ai(i)*x1ai(i)
             endif
            else
             st      = st &
                     + ( gp(i,j,k) - gp(i-1,j,k) ) * dx1bi(i)
            endif
           endif
           if (xhse) then
            if (nhy .eq. 0 .and. irestart .ne. 1) then
              gp1(i,j,k) =   rhoi &
                         * ( poo(i) - poo(i-1)) * dx1bi(i)
            endif
            st         = st + gp1(i,j,k)
           endif ! xhse
          endif  ! grav
!
!  3.  rotational pseudo-forces
!
          if(lgeom .eq. 3) then
           st      = st &
                   + ( s2oo(i  )     + s2po(i  ) &
                     + s2oo(i-1)     + s2po(i-1)   ) &
                   * ( u2(i  ,j,k) + u2(i  ,jp1,k  ) &
                     + u2(i-1,j,k) + u2(i-1,jp1,k  ) ) / 16.0 &
                   * rhoi * dg2bd1 (i) * g2ai(i)    * g2ai(i)
           st      = st &
                   + ( s3oo(i  )     + s3op(i  ) &
                     + s3oo(i-1)     + s3op(i-1)   ) &
                   * ( u3(i  ,j,k) + u3(i  ,j  ,kp1) &
                     + u3(i-1,j,k) + u3(i-1,j  ,kp1) ) / 16.0 &
                   * rhoi * dg31bd1(i) * g31ai(i)**2 * g32bi(j)
          endif ! lgeom
          if(xptmass) then
!
!  4.  gravitational point mass 
!
           if(lgeom .eq. 1) then
             r2        = ( x1a(i) - x1ptm )**2 &
                       + ( x2b(j) - x2ptm )**2 &
                       + ( x3b(k) - x3ptm )**2
             r2i       = ( x1a(i) - x1ptm ) &
                       / ( r2    * sqrt(r2   ) + tiny )
           endif ! CARTESIAN
           if(lgeom .eq. 2) then
             r2        = ( x1a(i) - x1ptm )**2 &
                       + ( x2b(j) - x2ptm )**2 &
                       + 2.0 * x2b(j) * x2ptm &
                       * ( 1.0 - cos(x3b(k) - x3ptm) )
             r2i       = ( x1a(i) - x1ptm ) &
                       / ( r2    * sqrt(r2   ) + tiny )
           endif ! CYLINDRICAL
           if(lgeom .eq. 3) then
!JH
!JH  correction due to Sean Matt -- enforce condition that radial
!JH  force due to point mass depend only on radius.
!JH
            if(x1ptm .eq. 0.0) then
             r2i = x1ai(i)**2
            else !x1ptm
             prpq      = x1a(i) * sin(x2b(j))
             pzpq      = x1a(i) * cos(x2b(j))
             prpt      = x1ptm  * sin(x2ptm)
             pzpt      = x1ptm  * cos(x2ptm)
!
             r2        = ( pzpq - pzpt )**2 &
                       + ( prpq - prpt )**2 &
                       + 2.0 * prpq * prpt &
                       * ( 1.0 - cos(x3b(k) - x3ptm) )
             r2i       = ( ( pzpq - pzpt ) * cos(x2b(j)) &
                       + ( prpq - prpt * cos(x3b(k) - x3ptm)) &
                       * sin(x2b(j))) / ( r2    * sqrt(r2   ) + tiny )
            endif ! x1ptm
           endif ! SPHERICAL
            st        = st - guniv * ptmass * r2i
          endif ! POINT MASS GRAVITY
!
!  5.  Radiation force
!
          if(lrad .eq. 1) then
           derdx1 = ( er(i  ,j  ,k  ) - er(i-1,j  ,k  ) ) &
                  * dx1bi(i)
           if(ldimen .eq. 1) then
             derdx2 = 0.0
           else
            derdx2 = ( er(i-1,j  ,k  ) - er(i-1,j-1,k  ) ) &
                   * dx2bi(j) &
                   + ( er(i  ,j  ,k  ) - er(i  ,j-1,k  ) ) &
                   * dx2bi(j)
            derdx2 = derdx2 * 0.5 * g2ai(i)
           endif ! ldimen = 1
           if(ldimen .lt. 3) then
            derdx3 = 0.0
           else
            derdx3 = ( er(i-1,j  ,k  ) - er(i-1,j  ,km1) ) &
                   * dx3bi(k) &
                   + ( er(i  ,j  ,k  ) - er(i  ,j  ,km1) ) &
                   * dx3bi(k)
            derdx3 = derdx3 * 0.5 * g31ai(i) * g32bi(j)
           endif ! ldimen < 3
           grad_er    = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i-1,j,k)
           rdmcp   = kapr(i  ,j,k)
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp &
                   / max(er(i  ,j,k)+er(i-1,j,k), tiny)
!
!
!          LP flux limiter
!
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
           st         = st - rhoi * flx_lim * derdx1
          endif ! lrad
!
!  6. Magnetic pressure (longitudinal Lorentz force)
!
          if(xmhd) then
             st = st - (d1b2oo(i) + d1b2po(i) + d1b3oo(i) + d1b3op(i)) &
                       * rhoi * dx1bi(i) * 0.25
          endif ! xmhd
!
! Include ionizing radiation pressure as dv_ioniz computed in coolchem
!
          w1(i,j,k) = u1(i,j,k) + dt * st ! + dv_ioniz(i,j,k)
!.......................................................................
!
! Perform an explicit update for v2 
!
         if(ldimen .ne. 1) then
          rhoi    = 2.0 / ( d(i,j-1,k) + d(i,j,k) )
!
!  1.  pressure gradient
!
          if(poo(i  ) .lt. 0.0D0) &
           write(2,"('Pressure = ',1pd12.4,' at i,j,k = ',3i4, &
                     ' on PID ',i3,'; nhy = ',i6)")poo(i),i,j,k,myid,nhy
          if((pmo(i  ) .lt. 0.0D0) .and. (j .eq. js)) &
           write(2,"('Pressure = ',1pd12.4,' at i,j,k = ',3i4, &
                     ' on PID ',i3,'; nhy = ',i6)")pmo(i),i,j-1,k,myid, &
                       nhy
          st      = - rhoi &
                  * ( poo(i)     - pmo(i)     ) * dx2bi(j) &
                  * g2bi(i)
!  2.  gravitational potential gradient
!
          if(xgrav .or. xgrvfft) then
           if(.not. xhse) then
            if(.not. xsphgrv) then
             st      = st &
                     + ( gp(i,j,k) - gp(i,j-1,k) ) * dx2bi(j) &
                     * g2bi(i)
            endif
           endif  
           if (xhse) then
            if (nhy .eq. 0 .and. irestart .ne. 1) then
              gp2(i,j,k) =  rhoi &
                         * (poo(i) - pmo(i)) * dx2bi(j) &
                         *  g2bi(i)
            endif
            st      = st + gp2(i,j,k)
           endif  ! xhse
          endif   ! grav
!
!  3.  rotational pseudo-force
!
          if(lgeom .ne. 1) then
           if(x2a(j) .ne. 0.0) then
             st      = st &
                     + ( s3oo(i)     + s3op(i) &
                       + s3mo(i)     + s3mp(i)      ) &
                     * ( u3(i,j  ,k) + u3(i,j  ,kp1) &
                       + u3(i,j-1,k) + u3(i,j-1,kp1) ) / 16.0 &
                     * rhoi       * g2bi(i) &
                     * dg32ad2(j) * g31bi(i) * g32ai(j) * g32ai(j)
!JH  -- changed dg32bd2 to dg32ad2 in line above
           endif
          endif ! lgeom
         else ! ldimen
          st = 0.0
         endif ! ldimen
          if(xptmass) then
!
!  4.  gravitational point mass
!
             if(lgeom .eq. 1) then
              r2        = ( x1b(i) - x1ptm )**2 &
                        + ( x2a(j) - x2ptm )**2 &
                        + ( x3b(k) - x3ptm )**2
              r2i       = ( x2a(j) - x2ptm ) &
                        / ( r2    * sqrt(r2   ) + tiny )
             endif ! CARTESIAN
             if(lgeom .eq. 2) then
              r2        = ( x1b(i) - x1ptm )**2 &
                        + ( x2a(j) - x2ptm )**2 &
                        + 2.0 * x2a(j) * x2ptm &
                        * ( 1.0 - cos(x3b(k) - x3ptm) )
              r2i       = ( x2a(j) - x2ptm * cos(x3b(k) - x3ptm)) &
                        / ( r2    * sqrt(r2   ) + tiny )
             endif ! CYLINDRICAL
             if(lgeom .eq. 3) then
!JH
!JH  correction due to Sean Matt -- enforces condition that theta
!JH  component of force be zero if point mass located at origin.
!JH
              if(x1ptm .eq. 0.0) then
               r2i       = 0.0
              else ! x1ptm
               prpq      = x1b(i) * sin(x2a(j))
               pzpq      = x1b(i) * cos(x2a(j))
               prpt      = x1ptm  * sin(x2ptm)
               pzpt      = x1ptm  * cos(x2ptm)
!
               r2        = ( pzpq - pzpt )**2 &
                         + ( prpq - prpt )**2 &
                         + 2.0 * prpq * prpt &
                         * ( 1.0 - cos(x3b(k) - x3ptm) )
               r2i       = (-( pzpq - pzpt ) * sin(x2a(j)) &
                         + ( prpq - prpt * cos(x3b(k) - x3ptm)) &
                         * cos(x2a(j))) / ( r2    * sqrt(r2   ) + tiny )
              endif ! x1ptm
             endif ! SPHERICAL
             st        = st - guniv * ptmass * r2i
          endif ! xptmass
!
!  5.  Radiation force 
!
          if(ldimen .ne. 1) then
          if(lrad .eq. 1) then
           derdx2 = ( er     (i  ,j  ,k  ) - er     (i  ,j-1,k  ) ) &
                  * dx2bi(j) * g2bi(i)
           derdx1 = ( er     (i  ,j-1,k  ) - er     (i-1,j-1,k  ) ) &
                  * dx1bi(i) &
                  + ( er     (i  ,j  ,k  ) - er     (i-1,j  ,k  ) ) &
                  * dx1bi(i)
           derdx1 = derdx1 * 0.5
           if(ldimen .eq. 2) then
            derdx3 = 0.0
           else
            derdx3 = ( er     (i  ,j-1,k  ) - er     (i  ,j-1,km1) ) &
                   * dx3bi(k) &
                   + ( er     (i  ,j  ,k  ) - er     (i  ,j  ,km1) ) &
                   * dx3bi(k)
            derdx3 = derdx3 * 0.5 * g31bi(i) * g32ai(j)
           endif ! ldimen
           grad_er = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i,j-1,k)     
           rdmcp   = kapr(i,j  ,k)
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp &
                   / max(er(i,j  ,k)+er(i,j-1,k), tiny)
!
!   LP flux limiter
!
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
           st        = st - rhoi * flx_lim * derdx2
          endif ! lrad
!
!  6. Magnetic pressure (longitudinal Lorentz force)
!
          if(xmhd) then
           st = st - (d2b3oo(i) + d2b3op(i) + d2b1oo(i) + d2b1po(i)) &
              * rhoi * dx2bi(j) * g2bi(i) * 0.25
          endif ! xmhd
         endif ! ldimen
!
         w2(i,j,k) = u2(i,j,k) + dt * st
!
!.......................................................................
!
! Perform an explicit update for v3
!
         st = 0.0
         if(ldimen .ne. 1) then
          rhoi = 2.0 / ( d(i,j,km1) + d(i,j,k) )
!
!  1.  pressure gradient
!
          st = - rhoi &
             * ( poo(i)    - pom(i)      ) * dx3bi(k) &
             * g31bi(i) * g32bi(j)
!
!  2.  gravitational potential gradient
!
          if(xgrav .or. xgrvfft) then
           if(.not. xhse) then
            if(.not. xsphgrv) then
             st      = st &
                     + ( gp(i,j,k) - gp(i,j,km1) ) * dx3bi(k) &
                     * g31bi(i) * g32bi(j)
            endif
           endif
           if (xhse) then
            if (nhy .eq. 0 .and. irestart .ne. 1) then
              gp3(i,j,k) =  rhoi &
                         * (poo(i)  - pom(i)) * dx3bi(k) &
                         * g31bi(i) * g32bi(j)
            endif
            st      = st + gp3(i,j,k)
           endif ! xhse
          endif  ! grav
         endif   ! ldimen
!
!  4.  gravitational point mass
!
          if(xptmass) then
           if(lgeom .eq. 1) then
             r2        = ( x1b(i) - x1ptm )**2 &
                       + ( x2b(j) - x2ptm )**2 &
                       + ( x3a(k) - x3ptm )**2
             r2i       = ( x3a(k) - x3ptm ) &
                       / ( r2    * sqrt(r2   ) + tiny )
           endif ! CARTESIAN
           if(lgeom .eq. 2) then
             r2        = ( x1b(i) - x1ptm )**2 &
                       + ( x2b(j) - x2ptm )**2 &
                       + 2.0 * x2b(j) * x2ptm &
                       * ( 1.0 - cos(x3a(k) - x3ptm) )
             r2i       = x2ptm * sin(x3a(k) - x3ptm) &
                       / ( r2    * sqrt(r2   ) + tiny )
           endif ! CYLINDRICAL
           if(lgeom .eq. 3) then
!JH
!JH  correction due to Sean Matt -- enforces condition that phi
!JH  component of force be zero if point mass located at origin.
!JH
            if(x1ptm .eq. 0.0) then
             r2i       = 0.0
            else ! x1ptm
             prpq      = x1b(i) * sin(x2b(j))
             pzpq      = x1b(i) * cos(x2b(j))
             prpt      = x1ptm  * sin(x2ptm)
             pzpt      = x1ptm  * cos(x2ptm)
!
             r2        = ( pzpq - pzpt )**2 &
                       + ( prpq - prpt )**2 &
                       + 2.0 * prpq * prpt &
                       * ( 1.0 - cos(x3a(k) - x3ptm) )
             r2i       = prpt * sin(x3a(k) - x3ptm) &
                       / ( r2    * sqrt(r2   ) + tiny )
            endif ! x1ptm
           endif ! SPHERICAL
             st        = st - guniv * ptmass * r2i
          endif ! xptmass
!
!  5. Radiation forces
!
         if(ldimen .ne. 1) then
          if(lrad .eq. 1) then
           derdx3 = ( er     (i  ,j  ,k  ) - er     (i  ,j  ,km1) ) &
                  * dx3bi(k) * g31bi(i) * g32bi(j)
           derdx1 = ( er     (i  ,j  ,km1) - er     (i-1,j  ,km1) ) &
                  * dx1bi(i) &
                  + ( er     (i  ,j  ,k  ) - er     (i-1,j  ,k  ) ) &
                  * dx1bi(i)
           derdx1 = derdx1 * 0.5
           derdx2 = ( er     (i  ,j  ,km1) - er     (i  ,j-1,km1) ) &
                  * dx2bi(j) &
                  + ( er     (i  ,j  ,k  ) - er     (i  ,j-1,k  ) ) &
                  * dx2bi(j)
           derdx2 = derdx2 * 0.5 * g2bi(i)  
           grad_er = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i,j,km1)
           rdmcp   = kapr(i,j,k  )       
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp &
                  / max(er(i,j,k  )+er(i,j,k-1), tiny)
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
           st        = st - rhoi * flx_lim * derdx3
          endif
!
          if(xmhd) then
!
!  6. Magnetic pressure (longitudinal Lorentz force)
!
           st = st - (d3b1oo(i) + d3b1po(i) + d3b2oo(i) + d3b2op(i)) &
              * rhoi * dx3bi(k) * g31bi(i) * g32bi(j) * 0.25
          endif ! xmhd
!
          endif ! ldimen
          w3(i,j,k) = u3(i,j,k) + dt * st
20      continue
30     continue
40    continue
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                F O R C E S                \\\\\\\\\!
!=======================================================================
!
!
