!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 A V I S C                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine avisc (ibeg,iend,jbeg,jend,kbeg,kend,dvdxmn &
                        ,w1,w2,w3,u1,u2,u3,s1,s2,s3)
!
!    mln:zeus3d.viscous <------------- artificial viscosity source terms
!                                                        ?????????, 19??
!
!    written by: Mike Norman
!    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
!    modified 2: February, 1990 by David Clarke; incorporated into
!                ZEUS3D
!    modified 3: June, 1992 by David Clarke; expunged "ISMIN", thereby
!                decreasing execution time by 30%.
!    modified 4: Oct., 1994 by Robert Fiedler to run in parallel on SGIs
!    modified 5: Totally rewritten 2/27/96 by RAF for ZEUS-MP.
!    modified 6: October 2005 by John Hayes; returned linear viscosity
!                terms, after ZEUS-3D.
!
!  PURPOSE: Computes the artificial viscosity source terms in the
!  momentum and energy equations.  i.e., it computes
!
!             dv / dt = -DIV(Q) / rho             for v1, v2, and v3
!      and    de / dt = -Q * delta(v) / delta(x)
!
!  This routine uses the von Neumann-Richtmyer form of the artificial
!  viscosity.  This means that geometric terms are not included in
!  DIV(Q).
!
!
!  LOCAL VARIABLES:
!    qqs        diagonal elements of viscous tensor.  Thus, this is a
!               linear treatment of the artificial viscosity.
!    dvelb      v1(i+1,j,k)-v1(i,j,k) for i-sweep
!               v2(i,j+1,k)-v2(i,j,k) for j-sweep
!               v3(i,j,k+1)-v3(i,j,k) for k-sweep
!    dvela      min ( zro, dvelb ) - ensures that only compressional
!               waves (shocks) are affected.
!    dvdxmn     min ( ( delta(v1) / delta(x1) ),
!                     ( delta(v2) / delta(x2) ),
!                     ( delta(v3) / delta(x3) ) )
!    w1,w2,w3   velocity values prior to viscosity update.
!    u1,u2,u3   velocity values after    viscosity update.
!    s1,s2,s3   updated momentum densities for transport step.
!
! BOUNDARY VALUES USED:
!
!  Macro defined  var   ii    oi    ij    oj    ik    ok
!  -------------  ---  ----  ----  ----  ----  ----  ----
!                  d   is-1        js-1        ks-1
!                  u1  is-1  ie+1
!                  u2              js-1  je+1
!                  u3                          ks-1  ke+1
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
      implicit none
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, km1, &
                  jm1
      real(rl) :: q1, q3, g3i, q2, q11p, q11m, q22p, q22m, q33p, q33m, &
                  qt
      real(rl) :: qqs,qqsm, qqsp, dvelas,dvelasm, dvelasp, dvdxmn, &
                  cs000, csm00, cs0m0, cs00m, dvelbs, dvelbsm
      real(rl) :: w1(in,jn,kn),w2(in,jn,kn),w3(in,jn,kn), &
                  u1(in,jn,kn),u2(in,jn,kn),u3(in,jn,kn), &
                  s1(in,jn,kn),s2(in,jn,kn),s3(in,jn,kn)
!
!
!-----------------------------------------------------------------------
!      Start artificial viscosity update.
!-----------------------------------------------------------------------
!
      if(xsubav) then
       q1 = avisc_dt * qcon
       q2 = avisc_dt * qlin
      else
       q1 = dt * qcon
       q2 = dt * qlin
      endif
!
      do 160 k=kbeg,kend
       if(ldimen .eq. 3) then
        km1 = k - 1
       else
        km1 = k
       endif
       do 150 j=jbeg,jend
        if(ldimen .gt. 1) then
         jm1 = j - 1
        else
         jm1 = j
        endif
        do 140 i=ibeg,iend
         if(qlin .ne. 0.0) then
          cs000 = sqrt(gamm1*gamma*e(i  ,j  ,k  )/d(i  ,j  ,k  ))
          csm00 = sqrt(gamm1*gamma*e(i-1,j  ,k  )/d(i-1,j  ,k  ))
          cs0m0 = sqrt(gamm1*gamma*e(i  ,jm1,k  )/d(i  ,jm1,k  ))
          cs00m = sqrt(gamm1*gamma*e(i  ,j  ,km1)/d(i  ,j  ,km1))
         endif
!
!      Do v1.
!
         dvelasm = min ( zro, w1(i  ,j,k) - w1  (i-1,j,k) )
         dvelas  = min ( zro, w1(i+1,j,k) - w1  (i  ,j,k) )
         if(qlin .eq. 0.0) then
          qqsm    = q1 * d(i-1,j,k) * dvelasm * dvelasm
          qqs     = q1 * d(i  ,j,k) * dvelas  * dvelas
          q3      = dvelas  * dx1ai(i)
         else ! qlin
          dvelbsm = w1(i  ,j,k) - w1  (i-1,j,k)
          dvelbs  = w1(i+1,j,k) - w1  (i  ,j,k)
          qqsm    = d(i-1,j,k)*dvelbsm &
                  * (q1 * dvelasm - q2 * csm00)
          qqs     = d(i  ,j,k)*dvelbs &
                  * (q1 * dvelas  - q2 * cs000)
          q3      = dvelbs  * dx1ai(i)
         endif ! qlin
!
         dvdxmn  = min ( dvdxmn, q3 )
         if(xiso .eqv. .false.) then
          if(xtotnrg .eqv. .false.) then
           e  (i,j,k) = e(i,j,k) - q3 * qqs
          else
           dvelasp = min ( zro, w1(i+2,j,k) - w1  (i+1,j,k) )
           qqsp    = q1 * d(i+1,j,k) * dvelasp * dvelasp
           q11p    = 0.5D0*(qqs + qqsp)
           q11m    = 0.5D0*(qqs + qqsm)
           qt      = (g2a(i+1)*g31a(i+1) * q11p * w1(i+1,j,k) &
                   -  g2a(i  )*g31a(i  ) * q11m * w1(i  ,j,k) ) &
                   *                                dvl1ai(i)
           e(i,j,k) = e(i,j,k) - qt
          endif
         endif
         u1 (i,j,k) = w1(i,j,k) &
                    - ( qqs  - qqsm   ) * dx1bi(i) &
                    * 2.0 / ( d(i-1,j,k) + d(i,j,k) )
!
!      Do v2.
!
         if(ldimen .gt. 1) then
          dvelasm = min ( zro, w2(i,j  ,k) - w2  (i,j-1,k) )
          dvelas  = min ( zro, w2(i,j+1,k) - w2  (i,j  ,k) )
          if(qlin .eq. 0.0) then
           qqsm    = q1 * d(i,j-1,k) * dvelasm  * dvelasm
           qqs     = q1 * d(i,j  ,k) * dvelas   * dvelas
           q3      = dvelas   * dx2ai(j) * g2bi(i)
          else ! qlin
           dvelbsm = w2(i,j  ,k) - w2  (i,j-1,k)
           dvelbs  = w2(i,j+1,k) - w2  (i,j  ,k)
           qqsm    = d(i,j-1,k) * dvelbsm &
                   * (q1 * dvelasm - q2 * cs0m0)
           qqs     = d(i,j  ,k) * dvelbs &
                   * (q1 * dvelas  - q2 * cs000)
           q3      = dvelbs   * dx2ai(j) * g2bi(i)
          endif ! qlin
          dvdxmn  = min ( dvdxmn, q3 )
          if(xiso .eqv. .false.) then
           if(xtotnrg .eqv. .false.) then
            e(i,j,k) = e(i,j,k) - q3 * qqs
           else
            dvelasp = min ( zro, w2(i,j+2,k) - w2  (i,j+1,k) )
            qqsp    = q1 * d(i,j+1,k) * dvelasp * dvelasp
            q22p    = 0.5D0*(qqs + qqsp)
            q22m    = 0.5D0*(qqs + qqsm)
            qt      = ( g32a(j+1) * q22p * w2(i,j+1,k) &
                      - g32a(j  ) * q22m * w2(i,j  ,k) ) &
                    *   g2bi(i)             * dvl2ai(j)
            e(i,j,k) = e(i,j,k) - qt
           endif ! xtotnrg
          endif
          u2(i,j,k) = w2(i,j,k) &
                    - ( qqs   - qqsm    ) * dx2bi(j) * g2bi(i) &
                    * 2.0 / ( d(i,j-1,k) + d(i,j,k) )
         else ! ldimen
          u2(i,j,k)  = w2(i,j,k)
         endif ! ldimen
!
!      Do v3.
!
         IF(LDIMEN .EQ. 3) THEN
         g3i     = g31bi(i) * g32bi(j)
         dvelasm = min ( zro, w3(i,j,k  ) - w3  (i,j,k-1) )
         dvelas  = min ( zro, w3(i,j,k+1) - w3  (i,j,k  ) )
         if(qlin .eq. 0.0) then
          qqsm    = q1 * d(i,j,k-1) * dvelasm  * dvelasm
          qqs     = q1 * d(i,j,k  ) * dvelas   * dvelas
          q3      = dvelas   * dx3ai(k) * g3i
         else ! qlin
          dvelbsm = w3(i,j,k  ) - w3(i,j,k-1)
          dvelbs  = w3(i,j,k+1) - w3(i,j,k  )
          qqsm    = d(i,j,k-1) * dvelbsm &
                  * (q1 * dvelasm - q2 * cs00m)
          qqs     = d(i,j,k  ) * dvelbs &
                  * (q1 * dvelas  - q2 * cs000)
          q3      = dvelbs   * dx3ai(k) * g3i
         endif ! qlin
         dvdxmn  = min ( dvdxmn, q3 )
         if(xtotnrg) q33(i,j,k) = qcon*dvelas*dvelas
         if((xiso .eqv. .false.) .and. &
            (xtotnrg .eqv. .false.)) then
!             if(x1dflag(i) .eqv. .false.) e(i,j,k) = e(i,j,k) - q3*qqs
              e(i,j,k) = e(i,j,k) - q3*qqs
         endif
         u3(i,j,k) = w3(i,j,k) &
                   - ( qqs   - qqsm    ) * dx3bi(k) * g3i &
                   *  2.0 / ( d(i,j,k-1) + d(i,j,k) )
         ELSE
         u3(i,j,k) = w3(i,j,k)
         ENDIF
!
! Save the updated momentum densities for the transport step, since 
! everything should be in the cache already.
!
         s1(i,j,k) = u1(i,j,k) * 0.5 * (d(i-1,j  ,k  ) + d(i,j,k))
         s2(i,j,k) = u2(i,j,k) * 0.5 * (d(i  ,jm1,k  ) + d(i,j,k)) &
                   * g2b(i)
         s3(i,j,k) = u3(i,j,k) * 0.5 * (d(i  ,j  ,km1) + d(i,j,k)) &
                   * g31b(i) * g32b(j)
140     continue
150    continue
160   continue
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 A V I S C                 \\\\\\\\\!
!=======================================================================
!
!
