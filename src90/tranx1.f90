!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                T R A N X 1                \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine tranx1 (ibeg,iend,jbeg,jend ,kbeg ,kend , dlo ,den, &
                   eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd ,deod, &
                   ero ,ern ,abo ,abn)
!
!    dac:zeus3d.tranx1 <----- transports zone-centred variables along x1
!    from jms:zeus2d.tranx1, mln:zeus04.tranz                  may, 1990
!
!    written by: David Clarke
!    modified 1: June 1992, by David Clarke; added the total energy
!                option originally designed by Byung-IL Jun.
!    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
!                for ZEUS-MP.
!    modified 3: Dec. 19, 1996 by Robert Fiedler; added radiation
!    modified 4: October 2005 by John Hayes; corrected transposed subscripts
!                in "dxo" array
!
!  PURPOSE:  Transports all zone centred variables in the 1-direction
!  only.  Currently transported are:
!
!                      mass   density
!                      energy density
!
!  The consistent advection algorithm, in which mass fluxes are used to
!  construct the fluxes of all variables across the interfaces, is used
!  for all hydrodynamical variables.  Thus, the mass fluxes are passed
!  to MOMX1 on order to transport the momenta as well.  The magnetic
!  field components are treated separately from the mass fluxes in CT.
!  Interpolations are done in-line.
!
!  INPUT VARIABLES: 
!    ibeg,iend,jbeg,jend,kbeg,kend  index ranges to cover.
!    dlo      mass            density at previous substep.
!    eod      specific energy density at previous substep; equals
!             (e+p)/d  if TOTAL_ENERGY is defined.
!
! BOUNDARY VALUES USED:
!
!    Macro defined  var   ii    oi    ij    oj    ik    ok
!    -------------  ---  ----  ----  ----  ----  ----  ----
!                    d   is-3  ie+2  js-1        ks-1
!                   e/d  is-2  ie+2
!                    u1  is-1  ie+1  js-1        ks-1
!    TOTAL_ENERGY    u2                    je+1
!    TOTAL_ENERGY    u3                                ke+1
!
!  OUTPUT VARIABLES:
!    den      updated mass            density.
!    edn      updated specific energy density.
!    mflx     mass flux (in the 1-direction)
!
!  LOCAL VARIABLES:
!    atwid    effective cross sectional area of the 1-interfaces
!    etwid    interpolated specific energy densities (e/d) at all
!             1-interfaces
!    eflx     energy density flux across all 1-interfaces  (reuse etwid)
!    dtwid    interpolated mass densities at all 1-interfaces
!    dflx     mass density flux across all 1-interfaces    (reuse dtwid)
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
      use bndry
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, n, kp1
!
      real(rl) :: dqm, dqp, xi,q1
      real(rl) :: atwid(ijkn), mflux(ijkn)
      real(rl) :: dtwid(ijkn), dd   (ijkn)
      real(rl) :: etwid(ijkn), deod (ijkn)
      real(rl), allocatable :: rtwid(:), dero (:)
      real(rl), allocatable :: xtwid(:,:), dxo (:,:)
!
      real(rl) :: mflx(in,jn,kn)
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn)
      real(rl) :: eod(in,jn,kn), edn(in,jn,kn)
      real(rl),optional :: ero(in,jn,kn), ern(in,jn,kn)
      real(rl),optional :: abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
 
      if (present(ero)) then
        allocate(rtwid(ijkn))
        allocate(dero (ijkn))
        allocate(xtwid(ijkn,nspec))
        allocate(dxo  (ijkn,nspec))
      endif
!-----------------------------------------------------------------------
!
! Compute time-centered area factors.
!
      do 10 i=ibeg-1,iend+1
       if(xvgrid) then
        atwid (i)  =       g2ah(i) * g31ah(i)
       else
        atwid (i)  =       g2a (i) * g31a (i)
       endif
10    continue
!
! Transport all zone-centered quantities in the 1 direction.
! Note that transporting v1 in MOMX1 will require the mass flux at 
! x1a(is-1).  To get it from the field variables, we need d at is-3.
! We also need mflx at js-1 and ks-1 for i=is,ie+1 for v2 and v3.
! Extend loop indices to compute mass fluxes beyond inner borders.
! Be careful to assign values to mflx, (and den, edn) only
! within the range (ibeg:iend,jbeg:jend,kbeg:kend) when these
! indices are not on the borders, so that they can't get overwritten
! when this routine is called with various index ranges.
!
      do 100 k=kbeg-1,kend
       do 90 j=jbeg-1,jend
!
!   Interpolate to obtain zone-centered quantities at zone faces.
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
        if (ibeg .eq. is) then  !  Need d(is-3) from neighbor.
         i         = is - 2
         dqm       = (dlo(i  ,j,k) - diib (j,k,3)) * dx1bi(i  )
         dqp       = (dlo(i+1,j,k) - dlo(i  ,j,k)) * dx1bi(i+1)
         dd(i)     = max ( dqm * dqp, zro ) &
                   * sign ( one, dqm + dqp ) &
                   / max ( abs ( dqm + dqp ), tiny )
         if(xiso .eqv. .false.) then
          deod  (i) = zro  ! Not valid, but we don't use it.
         endif
         if(lrad .ne.0) dero  (i) = zro
        endif ! ibeg
        do 30 i=max(ibeg-2,is-1),iend+1
         dqm       = (dlo(i  ,j,k) - dlo(i-1,j,k)) * dx1bi(i  )
         dqp       = (dlo(i+1,j,k) - dlo(i  ,j,k)) * dx1bi(i+1)
         dd(i)     = max ( dqm * dqp, zro ) &
                   * sign ( one, dqm + dqp ) &
                   / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm       = (abo(i  ,j,k,n)-abo(i-1,j,k,n))*dx1bi(i  )
           dqp       = (abo(i+1,j,k,n)-abo(i  ,j,k,n))*dx1bi(i+1)
           dxo(i,n) = max ( dqm * dqp, zro ) &
                     * sign ( one, dqm + dqp ) &
                     / max ( abs ( dqm + dqp ), tiny )
          enddo
         endif ! nspec
         if(xiso .eqv. .false.) then
          dqm       = (eod(i  ,j,k) - eod(i-1,j,k)) * dx1bi(i  )
          dqp       = (eod(i+1,j,k) - eod(i  ,j,k)) * dx1bi(i+1)
          deod(i)   = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
         endif ! xiso
         if(lrad .ne. 0) then
          dqm       = (ero(i  ,j,k) - ero(i-1,j,k)) * dx1bi(i  )
          dqp       = (ero(i+1,j,k) - ero(i  ,j,k)) * dx1bi(i+1)
          dero(i)   = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
         endif ! lrad
30      continue
!
!     2.  Choose time averaged, upwinded interface values.
!
        do 40 i=ibeg-1,iend+1
         xi        = ( v1  (i,j,k) - vg1(i) ) * dt
         q1        = sign ( haf, xi )
         dtwid (i) = ( 0.5 + q1 ) * ( dlo(i-1,j,k) &
                    + ( dx1a(i-1) - xi ) * dd   (i-1) ) &
                    + ( 0.5 - q1 ) * ( dlo(i  ,j,k) &
                    - ( dx1a(i  ) + xi ) * dd   (i  ) )
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(i,n) = ( 0.5 + q1 )*(abo(i-1,j,k,n) &
                      + ( dx1a(i-1) - xi )*dxo(i-1,n) ) &
                      + ( 0.5 - q1 )*(abo(i  ,j,k,n) &
                      - ( dx1a(i  ) + xi )*dxo(i  ,n) )
          enddo
         endif
         if(xiso .eqv. .false.) then
          etwid (i) = ( 0.5 + q1 ) * ( eod(i-1,j,k) &
                    + ( dx1a(i-1) - xi ) * deod (i-1) ) &
                    + ( 0.5 - q1 ) * ( eod(i  ,j,k) &
                    - ( dx1a(i  ) + xi ) * deod (i  ) )
         endif ! xiso
         if(lrad .ne. 0) then
          rtwid (i) = ( 0.5 + q1 ) * ( ero(i-1,j,k) &
                    + ( dx1a(i-1) - xi ) * dero (i-1) ) &
                    + ( 0.5 - q1 ) * ( ero(i  ,j,k) &
                    - ( dx1a(i  ) + xi ) * dero (i  ) )
         endif ! lrad
40      continue
!
!  For the purposes of consistent advection, construct the mass
!  flux across each 1-interface.  The mass flux will be used to create
!  the fluxes of all variables, including the momenta which are updated
!  in MOMX1.
!
        do 50 i=ibeg-1,iend+1
         mflux (i    ) = dtwid (i    ) * ( v1(i,j,k) - vg1(i) ) * dt
         dtwid (i    ) = mflux (i    ) * atwid (i    )
         if(nspec .gt. 1)then
          do n = 1, nspec
           xtwid(i,n) = xtwid(i,n) * dtwid(i)
          enddo
         endif
         if(xiso .eqv. .false.) then
          etwid (i    ) = dtwid (i    ) * etwid (i    )
         endif
         if(lrad .ne. 0)rtwid (i    ) = dtwid (i    ) * rtwid (i    )
50      continue
!
!  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
!  only for zones next to the inner borders.
!
        if ( (j.eq.js-1 .and. k.ge.kbeg)   .or. &
             (j.ge.jbeg .and. k.eq.ks-1)   .or. &
             (j.eq.js-1 .and. k.eq.ks-1) ) then
         if (ibeg.eq.is) mflx(is-1,j,k) = mflux (is-1)
         do 60 i=ibeg,iend
           mflx(i,j,k) = mflux (i      )
60       continue
         if (iend.eq.ie) mflx(ie+1,j,k) = mflux (ie+1)
        endif
!
!  Perform mass density and energy density advection.  Note that
!  the timestep "dt" is hidden the fluxes "dflx" and "eflx".
!  Do only zones inside (ibeg:iend,jbeg:jend,kbeg:kend).
!
        if (j.ge.jbeg .and. k.ge.kbeg) then
         if (ibeg.eq.is) mflx(is-1,j,k) = mflux (is-1)
          do 80 i=ibeg,iend
           mflx(i,j,k)= mflux (i      )
           if(xvgrid) then
            den(i,j,k) = ( dlo(i,j,k) * dvl1a(i) &
                       - dtwid(i+1    )+dtwid (i    ) )*dvl1ani(i)
           else
            den(i,j,k) = ( dlo(i,j,k) * dvl1a(i) &
                       - dtwid(i+1    )+dtwid (i    ) )*dvl1ai(i)
           endif
           if(nspec .gt. 1) then
            do n = 1, nspec
             if(xvgrid) then
              abn(i,j,k,n) = (dlo(i,j,k)*abo(i,j,k,n)*dvl1a(i) &
                           -  xtwid(i+1,n) + xtwid(i,n))*dvl1ani(i) &
                           / den(i,j,k)
             else
              abn(i,j,k,n) = (dlo(i,j,k)*abo(i,j,k,n)*dvl1a(i) &
                           -  xtwid(i+1,n) + xtwid(i,n))*dvl1ai(i) &
                           / den(i,j,k)
             endif
            enddo ! n
           endif ! nspec
           if(xiso .eqv. .false.) then
            if(xvgrid) then
             e(i,j,k) = ( e  (i,j,k) * dvl1a(i) &
                        - etwid(i+1) + etwid (i) ) * dvl1ani(i)
            else
             e(i,j,k) = ( e  (i,j,k) * dvl1a(i) &
                      -   etwid(i+1) + etwid (i) ) * dvl1ai(i)
            endif
!
! Compute e/d for the next substep.
!
            if(ldimen .eq. 3) then
             kp1 = k+1
            else
             kp1 = ks
            endif
            if(xtotnrg .eqv. .false.) then
             edn(i,j,k) = e(i,j,k) / den(i,j,k)
            else ! xtotnrg
             edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k) &
                        - gamm1 * ( (v1(i,j,k)+v1(i+1,j  ,k  ))**2 &
                                  + (v2(i,j,k)+v2(i  ,j+1,k  ))**2 &
                                  + (v3(i,j,k)+v3(i  ,j  ,kp1))**2 ) &
                                * 0.125
            endif ! xtotnrg
           endif ! xiso
           if(lrad .ne. 0) then
            if(xvgrid) then
             ern(i,j,k) = ( ero(i,j,k) * dlo(i,j,k) * dvl1a(i) &
                        -   rtwid(i+1) + rtwid(i) ) * dvl1ani(i)
            else
             ern(i,j,k) = ( ero(i,j,k) * dlo(i,j,k) * dvl1a(i) &
                        -   rtwid(i+1) + rtwid(i) ) * dvl1ai(i)
            endif ! xvgrid
!
! Work with er/d for the next sweep.
!
            ern(i,j,k) = ern(i,j,k) / den(i,j,k)
           endif ! lrad
80        continue
          if (iend.eq.ie) mflx(ie+1,j,k) = mflux (ie+1)
         endif
90      continue
100    continue
!
      if (allocated(rtwid))  deallocate(rtwid)
      if (allocated(xtwid))  deallocate(xtwid)
      if (allocated(dero ))  deallocate(dero )
      if (allocated(dxo  ))  deallocate( dxo )
      return
end subroutine tranx1
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                T R A N X 1                \\\\\\\\\!
!=======================================================================
!
!
