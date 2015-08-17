!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                T R A N X 2                \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine tranx2 (ibeg,iend,jbeg,jend,kbeg,kend,dlo,den,eod,edn &
                   ,mflx,atwid,ero,ern,abo,abn)
!
!    dac:zeus3d.tranx2 <----- transports zone-centred variables along x2
!    from jms:zeus2d.tranx2, mln:zeus04.tranz                  may, 1990
!
!    written by: David Clarke
!    modified 1: June 1992, by David Clarke; added the total energy
!                option originally designed by Byung-IL Jun.
!    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
!                for ZEUS-MP.
!    modified 3: Aug. 8, 1996 by Robert Fiedler; unrolled i-loops.
!    modified 4: Dec. 23, 1996 by Robert Fiedler; added radiation.
!    modified 5: Sep. 8, 2003 by John Hayes; re-rolled i-loops.
!
!  PURPOSE:  Transports all zone centred variables in the 2-direction
!  only.  Currently transported are:
!
!                      mass   density
!                      energy density
!
!  The consistent advection algorithm, in which mass fluxes are used to
!  construct the fluxes of all variables across the interfaces, is used
!  for all hydrodynamical variables.  Thus, the mass fluxes are passed
!  to MOMX2 on order to transport the momenta as well.  The magnetic
!  field components are treated separately from the mass fluxes in CT.
!  Interpolations are done in-line.
!
!  INPUT VARIABLES: 
!    ibeg,iend,jbeg,jend,kbeg,kend  index ranges to cover.
!    dlo      mass            density at previous substep.
!    eod      specific energy density at previous substep; equals
!             (e+p)/d  if TOTAL_ENERGY is defined.
!    ero      specific radiation energy density at previous substep.
!
! BOUNDARY VALUES USED:
!
!    Macro defined  var   ii    oi    ij    oj    ik    ok
!    -------------  ---  ----  ----  ----  ----  ----  ----
!                    d   is-1        js-3  je+2  ks-1
!                   e/d              js-2  je+2
!    TOTAL_ENERGY    u1        ie+1                         
!                    u2  is-1        js-1  je+1  ks-1
!    TOTAL_ENERGY    u3                                ke+1
!
!  OUTPUT VARIABLES:
!    den      updated mass            density.
!    edn      updated specific energy density.
!    ern      updated specific radiation energy density.
!    mflx     mass flux (in the 2-direction)
!
!  LOCAL VARIABLES:
!    atwid    effective cross sectional area of the 1-interfaces
!    etwid    interpolated specific energy densities (e/d) at all
!             2-interfaces
!    rtwid    interpolated specific radiation energy densities at all
!             2-interfaces
!    eflx     energy density flux across all 2-interfaces  (reuse etwid)
!    dtwid    interpolated mass densities at all 2-interfaces
!    dflx     mass density flux across all 2-interfaces    (reuse dtwid)
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
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, kp1, l
      integer  :: kstart, n, m
!
      real(rl) :: dqm,dqp, xi,q1, atwid(ijkn)
!
      real(rl) :: mflux(ijkn)
      real(rl) :: dtwid(ijkn), dd  (ijkn)
      real(rl) :: etwid(ijkn), deod(ijkn)
!
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn)
      real(rl) :: eod(in,jn,kn), edn(in,jn,kn)
      real(rl),optional :: ero(in,jn,kn), ern(in,jn,kn)
      real(rl),optional :: abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
!
      real(rl), allocatable :: rtwid(:), dero(:)
      real(rl), allocatable :: dxo(:,:), xtwid(:,:)
      if (present(ero)) then 
        allocate(rtwid(ijkn))
        allocate(dero (ijkn))
        allocate(dxo  (ijkn,nspec))
        allocate(xtwid(ijkn,nspec))
      endif
!-----------------------------------------------------------------------
!
! Compute time-centered area factors.
!
      do 10 i=ibeg-1,iend+1
        atwid (i)         = g31b(i) * dx1a(i) * dvl1ai(i)
10    continue
!
! Transport all zone-centered quantities in the 2 direction.
!
      if(ldimen .ne. 3) then
       kstart = ks
      else
       kstart = kbeg-1
      endif
      do 2100 k=kstart,kend
       if(ldimen .eq. 3) then
        kp1 = k+1
       else
        kp1 = ks
       endif
       if ( (k.ge.kbeg) .or. (k.eq.ks-1) ) then
       do 1090 i=ibeg-1,iend
!
!   Interpolate to obtain zone-centered quantities at zone faces.
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
        if (jbeg .eq. js) then  !  Need d(js-3) from neighbor.
         j         = js - 2
         dqm       = (dlo(i,j  ,k) - dijb (i,k,3)) * dx2bi(j  )
         dqp       = (dlo(i,j+1,k) - dlo(i,j  ,k)) * dx2bi(j+1)
         dd  (j) = max ( dqm * dqp, zro ) &
                   * sign ( one, dqm + dqp ) &
                   / max ( abs ( dqm + dqp ), tiny )
         if(xiso .eqv. .false.) then
          deod(j) = zro  ! Not valid, but we don't use it.
         endif
         if(lrad .ne. 0) then
           dero(j) = zro
         endif
        endif
        do 1030 j=max(jbeg-2,js-1),jend+1
         dqm       = (dlo(i,j  ,k) - dlo(i,j-1,k)) * dx2bi(j  )
         dqp       = (dlo(i,j+1,k) - dlo(i,j  ,k)) * dx2bi(j+1)
         dd  (j) = max ( dqm * dqp, zro ) &
                   * sign ( one, dqm + dqp ) &
                   / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm        = (abo(i,j  ,k,n)-abo(i,j-1,k,n))*dx2bi(j  )
           dqp        = (abo(i,j+1,k,n)-abo(i,j  ,k,n))*dx2bi(j+1)
           dxo(j,n) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
          enddo ! n
         endif
         if(xiso .eqv. .false.) then
          dqm       = (eod(i,j  ,k) - eod(i,j-1,k)) * dx2bi(j  )
          dqp       = (eod(i,j+1,k) - eod(i,j  ,k)) * dx2bi(j+1)
          deod(j) = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
         endif ! xiso
         if(lrad .ne. 0) then
          dqm       = (ero(i,j  ,k)-ero(i,j-1,k))*dx2bi(j  )
          dqp       = (ero(i,j+1,k)-ero(i,j  ,k))*dx2bi(j+1)
          dero(j) = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
         endif
1030    continue
!
!     2.  Choose time averaged, upwinded interface value.
!
!  For the purposes of consistent advection, construct the mass
!  flux across each 1-interface.  The mass flux will be used to create
!  the fluxes of all variables, including the momenta which are updated
!  in MOMX1.
!
        do 1040 j=jbeg-1,jend+1
         xi          = ( v2  (i,j  ,k) - vg2(j  ) ) * dt * g2bi(i)
         q1          = sign ( haf, xi )
         dtwid (j) = ( 0.5 + q1 ) * ( dlo(i,j-1,k) &
                     + ( dx2a(j-1) - xi ) * dd   (j-1) ) &
                     + ( 0.5 - q1 ) * ( dlo(i,j  ,k) &
                     - ( dx2a(j  ) + xi ) * dd   (j  ) )
!
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(j,n) = (0.5 + q1)*(abo(i,j-1,k,n) &
                        + (dx2a(j-1) - xi ) * dxo(j-1,n) ) &
                        + (0.5 - q1)*(abo(i,j  ,k,n) &
                        - (dx2a(j  ) + xi ) * dxo(j,n  ) )
          enddo
         endif ! nspec
!
         mflux (j) = dtwid (j) * ( v2(i,j  ,k) - vg2(j  ) ) * dt
!
         if(xvgrid) then
          dtwid (j) = mflux (j) * atwid (i) * g32ah(j  )
         else
          dtwid (j) = mflux (j) * atwid (i) * g32a (j  )
         endif
!
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(j,n) = xtwid(j,n)*dtwid(j)
          enddo
         endif
!
         if(xiso .eqv. .false.) then
          etwid (j)= ( 0.5 + q1 ) * ( eod(i,j-1,k) &
                     + ( dx2a(j-1) - xi ) * deod (j-1) ) &
                     + ( 0.5 - q1 ) * ( eod(i,j  ,k) &
                     - ( dx2a(j  ) + xi ) * deod (j  ) )
!
          etwid (j)= dtwid (j) * etwid (j)
         endif ! xiso
         if(lrad .ne. 0) then
          rtwid (j)= ( 0.5 + q1 ) * ( ero(i,j-1,k) &
                     + ( dx2a(j-1) - xi ) * dero (j-1) ) &
                     + ( 0.5 - q1 ) * ( ero(i,j,k) &
                     - ( dx2a(j  ) + xi ) * dero (j) )
          rtwid (j)= dtwid(j)*rtwid(j)
         endif
1040    continue
!
!  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
!  only for zones next to the inner orders.
!
        if( (i.eq.is-1 .and. k.ge.kbeg) .or. &
            (i.ge.ibeg .and. k.eq.ks-1) .or. &
            (i.eq.is-1 .and. k.eq.ks-1) ) then
         if (jbeg .eq. js) mflx(i,js-1,k) = mflux (js-1)
         do 1060 j=jbeg,jend
          mflx(i,j,k) = mflux(j)
1060     continue
         if (jend .eq. je) mflx(i,je+1,k) = mflux (je+1)
        endif
!
!  Perform mass density and energy density advection.  Note that
!  the timestep "dt" is hidden the fluxes "dflx" and "eflx".
!
        if (k.ge.kbeg .and. i.ge.ibeg) then
         if (jbeg .eq. js) mflx(i,js-1,k) = mflux (js-1)
         if (jend .eq. je) mflx(i,je+1,k) = mflux (je+1)
         do 1080 j=jbeg,jend
          mflx(i,j,k) = mflux (j)
!
          if(xvgrid) then
!
            den(i,j,k) = ( dlo(i,j,k) * dvl2a(j) &
                       - dtwid (j+1)  + dtwid (j)  ) * dvl2ani(j)
!
          else
!
            den(i,j,k) = ( dlo(i,j,k) * dvl2a(j) &
                       - dtwid (j+1)  + dtwid (j)  ) * dvl2ai(j)
!
          endif ! xvgrid
!
          if(nspec .gt. 1) then
           do n = 1, nspec
            if(xvgrid) then
!
              abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl2a(j) &
                           - xtwid(j+1,n)+xtwid(j,n))*dvl2ani(j) &
                           / den(i,j,k)
!
            else ! xvgrid
!
              abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl2a(j) &
                           - xtwid(j+1,n)+xtwid(j,n))*dvl2ai(j) &
                           / den(i,j,k)
            endif ! xvgrid
           enddo ! n
          endif ! nspec
!
          if(xiso .eqv. .false.) then
           if(xvgrid) then
!
             e  (i,j,k) = ( e  (i,j,k) * dvl2a(j) &
                        - etwid (j+1)  + etwid (j)  ) * dvl2ani(j)
!
           else
!
             e  (i,j,k) = ( e  (i,j,k) * dvl2a(j) &
                        - etwid (j+1)  + etwid (j)  ) * dvl2ai(j)
           endif
!
! Compute e/d for the next substep.
!
           if(xtotnrg .eqv. .false.) then
            edn(i,j,k) = e(i,j,k) / den(i,j,k)
           else ! xtotnrg
            edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k) &
                       - gamm1 * ( ( v1(i,j,k) + v1(i+1,j  ,k  ) )**2 &
                                 + ( v2(i,j,k) + v2(i  ,j+1,k  ) )**2 &
                                 + ( v3(i,j,k) + v3(i  ,j  ,kp1) )**2 &
                                 ) * 0.125
            if(xgrav .or. xgrvfft) then
             edn(i,j,k) = edn(i,j,k) &
                        + gamm1 * gp(i,j,k)
            endif
           endif ! xtotnrg
          endif ! xiso
!
          if(lrad .ne. 0) then
           if(xvgrid) then
            ern(i  ,j,k) = ( ero(i  ,j,k) * dlo(i  ,j,k) * dvl2a(j) &
                         - rtwid (j+1) + rtwid (j) ) * dvl2ani(j)
           else
            ern(i  ,j,k) = ( ero(i  ,j,k) * dlo(i  ,j,k) * dvl2a(j) &
                         - rtwid (j+1) + rtwid (j) ) * dvl2ai(j)
           endif
!
! Work with er/d for the next sweep.
!
           ern(i  ,j,k) = ern(i  ,j,k) / den(i  ,j,k)
          endif
1080     continue
        endif
1090   continue
!
       endif  !  k=ks-1 || k>=kbeg
2100   continue
!
       
      if (allocated(rtwid))  deallocate(rtwid)
      if (allocated(dero ))  deallocate(dero )
      if (allocated(dxo  ))  deallocate(dxo  )
      if (allocated(xtwid))  deallocate(xtwid)
      return
end subroutine tranx2
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                T R A N X 2                \\\\\\\\\!
!=======================================================================
!
!
