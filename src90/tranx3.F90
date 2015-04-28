!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                T R A N X 3                \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine tranx3 (ibeg,iend,jbeg,jend,kbeg,kend,dlo,den,eod,edn &
                  ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod,ero,ern,abo,abn)
!
!    dac:zeus3d.tranx3 <----- transports zone-centred variables along x3
!    from jms:zeus2d.tranx2, mln:zeus04.tranz                  may, 1990
!
!    written by: David Clarke
!    modified 1: June 1992, by David Clarke; added the total energy
!                option originally designed by Byung-IL Jun.
!    modified 2: Feb. 22, 1996 by Robert Fiedler; completely rewritten
!                for ZEUS-MP.
!    modified 3: Aug. 6,  1996 by Robert Fiedler; rearranged for
!                maximum efficiency, and unrolled i-loops.
!    modified 3: Dec. 18,  1996 by Robert Fiedler; added radiation.
!    modified 4: by John Hayes; rewrote for F90 and removed CPP logic
!    modified 5: by John Hayes; restored metric factors multiplying "xi"
!    modified 6: by John Hayes; corrected transposed array subscripts in
!                "dxo" array
!
!  PURPOSE:  Transports all zone centred variables in the 3-direction
!  only.  Currently transported are:
!
!                      mass   density
!                      energy density
!
!  The consistent advection algorithm, in which mass fluxes are used to
!  construct the fluxes of all variables across the interfaces, is used
!  for all hydrodynamical variables.  Thus, the mass fluxes are passed
!  to MOMX3 in order to transport the momenta as well.  The magnetic
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
!                    d   is-1        js-1        ks-3  ke+2
!                   e/d                          ks-2  ke+2
!    TOTAL_ENERGY    u1        ie+1
!    TOTAL_ENERGY    u2                    je+1
!                    u3  is-1        js-1        ks-1  ke+1
!
!  OUTPUT VARIABLES:
!    den      updated mass            density.
!    edn      updated specific energy density.
!    mflx     mass flux (in the 3-direction)
!
!  LOCAL VARIABLES:
!    atwid    effective cross sectional area of the 1-interfaces
!    atwid2   effective cross sectional area of the 2-interfaces
!    etwid    interpolated specific energy densities (e/d) at all
!             3-interfaces
!    eflx     energy density flux across all 3-interfaces (reuse etwid)
!    dtwid    interpolated mass densities at all 3-interfaces
!    dflx     mass density flux across all 3-interfaces   (reuse dtwid)
!
!  EXTERNALS:
!
      use real_prec
      use config
      use param
      use root
      use field
      use bndry
      use grid
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, n
!
      real(rl) :: dqm, dqp, xi, q1
!
      real(rl) :: atwid2(ijkn),  atwid(ijkn)
      real(rl) :: mflux (ijkn,1)
      real(rl) :: dtwid (ijkn,1), dd  (ijkn,1)
      real(rl) :: etwid (ijkn,1), deod(ijkn,1)
      real(rl), allocatable :: rtwid (:,:), dero(:,:)
      real(rl), allocatable :: xtwid(:,:),dxo(:,:)
!
      real(rl) :: mflx(in,jn,kn)
      real(rl) :: dlo (in,jn,kn), den(in,jn,kn)
      real(rl) :: eod (in,jn,kn), edn(in,jn,kn)
      real(rl), optional :: ero (in,jn,kn), ern(in,jn,kn)
      real(rl), optional :: abo (in,jn,kn,nspec), abn(in,jn,kn,nspec)

      if (present(ero)) then
        allocate(rtwid (ijkn,1))
        allocate(dero  (ijkn,1))
        allocate(xtwid(ijkn,nspec))
        allocate(dxo  (ijkn,nspec))
      endif
!-----------------------------------------------------------------------
!
! Compute time-centered area factors.
!
      do 10 i=ibeg-1,iend+1
       atwid (i)           = g2b(i) * dx1a(i) * dvl1ai(i)
10    continue
      do 20 j=jbeg-1,jend+1
       atwid2(j)           = dx2a(j) * dvl2ai(j)
20    continue
!
! Transport all zone-centered quantities in the 3 direction.
!
      do 2100 j=jbeg-1,jend
       if ( (j.ge.jbeg) .or. (j.eq.js-1) ) then
!.......................................................................
!
! Split off the i=ibeg-1 iteration to ease unrolling.
!
        i = ibeg - 1
!
!   Interpolate to obtain zone-centered quantities at zone faces.
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
        if (kbeg .eq. ks) then  !  Need d(ks-3) from neighbor.
         k        = ks - 2
         dqm      = (dlo(i  ,j,k  ) - dikb (i  ,j,3)) * dx3bi(k  )
         dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
         dd(k,1)  = max ( dqm * dqp, zro ) &
                  * sign ( one, dqm + dqp ) &
                  / max ( abs ( dqm + dqp ), tiny )
        endif
        do 30 k=max(kbeg-2,ks-1),kend+1
         dqm      = (dlo(i  ,j,k  ) - dlo(i  ,j,k-1)) * dx3bi(k  )
         dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
         dd(k,1)  = max ( dqm * dqp, zro ) &
                  * sign ( one, dqm + dqp ) &
                  / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm         = (abo(i,j,k  ,n)-abo(i,j,k-1,n))*dx3bi(k  )
           dqp         = (abo(i,j,k+1,n)-abo(i,j,k  ,n))*dx3bi(k+1)
           dxo(k,n) = max ( dqm * dqp, zro ) &
                      * sign ( one, dqm + dqp ) &
                      / max ( abs ( dqm + dqp ), tiny )
          enddo
         endif
30      continue
!
!     2.  Choose time averaged, upwinded interface value.
!
!  For the purposes of consistent advection, construct the mass
!  flux across each 1-interface.  The mass flux will be used to create
!  the fluxes of all variables, including the momenta which are updated
!  in MOMX1.
!
        do 40 k=kbeg-1,kend+1
         xi          = ( v3  (i  ,j,k  ) - vg3(k  ) ) * dt &
                     * g31bi(i) * g32bi(j)
         q1          = sign ( haf, xi )
         dtwid (k,1) = ( 0.5 + q1 ) * ( dlo(i  ,j,k-1) &
                     + ( dx3a(k-1) - xi ) * dd   (k-1,1) ) &
                     + ( 0.5 - q1 ) * ( dlo(i  ,j,k  ) &
                     - ( dx3a(k  ) + xi ) * dd   (k  ,1) )
!
         mflux (k,1) = dtwid (k,1) * ( v3(i  ,j,k) - vg3(k) ) * dt
40      continue
!
!  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
!  only for zones along the inner borders.
!
        if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
         do 60 k=kbeg,kend
          mflx(i  ,j,k) = mflux (k   ,1)
60       continue
         if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
!
!.......................................................................
!
         do 1090 i=ibeg,iend
!
!   Interpolate to obtain zone-centered quantities at zone faces.
!
!     1.  Evaluate monotonised, van Leer differences across the zone.
!
          if (kbeg .eq. ks) then  !  Need d(ks-3) from neighbor.
           k        = ks - 2
           dqm      = (dlo(i  ,j,k  ) - dikb (i  ,j,3)) * dx3bi(k  )
           dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
           dd(k,1)  = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
           if(xiso .eqv. .false.) then
            deod(k,1) = zro  ! Not valid, but we don't use it.
           endif
           if(lrad .ne. 0) dero(k,1) = zro
          endif
          do 1030 k=max(kbeg-2,ks-1),kend+1
           dqm      = (dlo(i  ,j,k  ) - dlo(i  ,j,k-1)) * dx3bi(k  )
           dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
           dd(k,1)  = max ( dqm * dqp, zro ) &
                    * sign ( one, dqm + dqp ) &
                    / max ( abs ( dqm + dqp ), tiny )
           if(nspec .gt. 1) then
            do n = 1, nspec
             dqm      = (abo(i,j,k  ,n)-abo(i,j,k-1,n))*dx3bi(k  )
             dqp      = (abo(i,j,k+1,n)-abo(i,j,k  ,n))*dx3bi(k+1)
             dxo(k,n) = max ( dqm * dqp, zro ) &
                        * sign ( one, dqm + dqp ) &
                        / max ( abs ( dqm + dqp ), tiny )
            enddo
           endif
           if(xiso .eqv. .false.) then
            dqm      = (eod(i  ,j,k  ) - eod(i  ,j,k-1)) * dx3bi(k  )
            dqp      = (eod(i  ,j,k+1) - eod(i  ,j,k  )) * dx3bi(k+1)
            deod(k,1)= max ( dqm * dqp, zro ) &
                     * sign ( one, dqm + dqp ) &
                     / max ( abs ( dqm + dqp ), tiny )
           endif
           if(lrad .ne. 0) then
            dqm      = (ero(i  ,j,k  ) - ero(i  ,j,k-1)) * dx3bi(k  )
            dqp      = (ero(i  ,j,k+1) - ero(i  ,j,k  )) * dx3bi(k+1)
            dero(k,1)= max ( dqm * dqp, zro ) &
                     * sign ( one, dqm + dqp ) &
                     / max ( abs ( dqm + dqp ), tiny )
           endif
1030      continue
!
!     2.  Choose time averaged, upwinded interface value.
!
!  For the purposes of consistent advection, construct the mass
!  flux across each 1-interface.  The mass flux will be used to create
!  the fluxes of all variables, including the momenta which are updated
!  in MOMX1.
!
          do 1040 k=kbeg-1,kend+1
           xi          = ( v3  (i  ,j,k  ) - vg3(k  ) ) * dt &
                       * g31bi(i) * g32bi(j)
           q1          = sign ( haf, xi )
           dtwid (k,1) = ( 0.5 + q1 ) * ( dlo(i  ,j,k-1) &
                       + ( dx3a(k-1) - xi ) * dd   (k-1,1) ) &
                       + ( 0.5 - q1 ) * ( dlo(i  ,j,k  ) &
                       - ( dx3a(k  ) + xi ) * dd   (k  ,1) )
!
           mflux (k,1) = dtwid (k,1) * ( v3(i  ,j,k) - vg3(k) ) * dt
           dtwid (k,1) = mflux (k,1) * atwid (i  ) * atwid2(j)
           if(nspec .gt. 1) then
            do n = 1, nspec
             xtwid(k,n) = ( 0.5 + q1 ) * (abo(i,j,k-1,n) &
                          + ( dx3a(k-1) - xi ) * dxo(k-1,n) ) &
                          + ( 0.5 - q1 ) * (abo(i,j,k  ,n) &
                          - ( dx3a(k  ) + xi ) * dxo(k  ,n) )
            enddo
            do n = 1, nspec
             xtwid(k,n) = xtwid(k,n)*dtwid(k,1)
            enddo
           endif ! nspec
           if(xiso .eqv. .false.) then
            etwid (k,1) = ( 0.5 + q1 ) * ( eod(i  ,j,k-1) &
                        + ( dx3a(k-1) - xi ) * deod (k-1,1) ) &
                        + ( 0.5 - q1 ) * ( eod(i  ,j,k  ) &
                        - ( dx3a(k  ) + xi ) * deod (k  ,1) )
!
            etwid (k,1) = dtwid (k,1) * etwid (k,1)
           endif
           if(lrad .ne. 0) then
            rtwid (k,1) = ( 0.5 + q1 ) * ( ero(i  ,j,k-1) &
                        + ( dx3a(k-1) - xi ) * dero (k-1,1) ) &
                        + ( 0.5 - q1 ) * ( ero(i  ,j,k  ) &
                        - ( dx3a(k  ) + xi ) * dero (k  ,1) )
!
            rtwid (k,1) = dtwid (k,1) * rtwid (k,1)
           endif
1040      continue
!
!  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
!  only for zones along the inner borders.
!
          if (j.eq.js-1) then
           if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
           do 1060 k=kbeg,kend
            mflx(i  ,j,k) = mflux (k   ,1)
1060       continue
           if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
          endif
!
!  Perform mass density and energy density advection.  Note that
!  the timestep "dt" is hidden the fluxes "dtwid" and "etwid".
!
          if (j.ge.jbeg) then
           if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
           do 1080 k=kbeg,kend
            mflx(i  ,j,k)= mflux (k   ,1)
            if(xvgrid) then
             den(i,j,k) = ( dlo(i  ,j,k) * dvl3a(k) &
                        - dtwid(k+1,1)+dtwid (k,1)) * dvl3ani(k)
            else
             den(i,j,k) = ( dlo(i  ,j,k) * dvl3a(k) &
                        - dtwid(k+1,1)+dtwid (k,1)) * dvl3ai(k)
            endif
            if(nspec .gt. 1) then
             do n = 1, nspec
              if(xvgrid) then
               abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl3a(k) &
                            -  xtwid(k+1,n)+xtwid(k,n)) * &
                               dvl3ani(k)/den(i,j,k)
              else
               abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl3a(k) &
                            -  xtwid(k+1,n)+xtwid(k,n)) * &
                               dvl3ai(k)/den(i,j,k)
              endif ! xvgrid
             enddo
            endif ! nspec
            if(xiso .eqv. .false.) then
             if(xvgrid) then
              e(i,j,k) = ( e  (i  ,j,k) * dvl3a(k) &
                       - etwid (k+1,1) + etwid (k,1) ) * dvl3ani(k)
             else
              e(i,j,k) = ( e  (i  ,j,k) * dvl3a(k) &
                       - etwid (k+1,1) + etwid (k,1) ) * dvl3ai(k)
             endif
!
! Compute e/d for the next substep.
!
             if(xtotnrg .eqv. .false.) then
              edn(i,j,k) =         e(i,j,k) / den(i,j,k)
             else !xtotnrg
              edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k) - gamm1 &
                         * ( &
                           ( v1(i  ,j,k) + v1(i+1,j  ,k  ) )**2 &
                         + ( v2(i  ,j,k) + v2(i  ,j+1,k  ) )**2 &
                         + ( v3(i  ,j,k) + v3(i  ,j  ,k+1) )**2 &
                            ) * 0.125
             endif ! xtotnrg
            endif ! xiso
            if(lrad .ne. 0) then
             if(xvgrid) then
              ern(i,j,k) = (ero(i,j,k) * dlo(i  ,j,k) * dvl3a(k) &
                           - rtwid(k+1,1)+rtwid (k,1) ) * dvl3ani(k)
             else
              ern(i,j,k) = (ero(i,j,k) * dlo(i  ,j,k) * dvl3a(k) &
                         -   rtwid(k+1,1)+rtwid (k,1) ) * dvl3ai(k)
            endif
!
! Work with er/d for the next sweep.
!
            ern(i  ,j,k) = ern(i  ,j,k) / den(i  ,j,k)
           endif ! lrad
1080      continue
          if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
         endif
1090    continue
       endif  !  j>=jbeg || j=js-1
2100  continue
!
      if (allocated(rtwid))  deallocate(rtwid)
      if (allocated(dero ))  deallocate(dero )
      if (allocated(xtwid))  deallocate(xtwid)
      if (allocated(dxo  ))  deallocate(dxo  )
      return
end subroutine tranx3
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                T R A N X 3                \\\\\\\\\\
!
!=======================================================================
!
!
