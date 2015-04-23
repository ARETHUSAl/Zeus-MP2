!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                  N U D T                  \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine nudt
!
!    mln:zeus3d.nudt <------------------------ mhd time step calculation
!                                                              may, 1986
!
!    written by: Mike Norman
!    modified 1: May, 1986 by David Clarke; adapted for mhd.
!    modified 2: October, 1987 by Mike Norman; reworked for covariant
!                formulation.
!    modified 3: June, 1988 by Jim Stone; incorporated into ZEUS2D.
!    modified 4: February, 1990 by David Clarke; incorporated into
!                ZEUS3D.
!    modified 5: September, 1990 by David Clarke; moved magnetic fields
!                to face-centres.
!    modified 6: March 5, 1996 by Robert Fiedler; completely rewritten
!                for ZEUS-MP.
!    modified 7: Aug. 23, 1996 by Robert Fiedler; minor change to 
!                substep counter.
!    modified 8: Dec. 19, 1996 by Robert Fiedler; added radiation 
!                diffusion time step.
!    modified 9: Feb. 14, 1997 by RAF; corrected for non-hydro and
!                lower dimensionality options.
!    modified 10: kluge for Garching T3E (search on M-MML) 26 Feb 98
!    modified 11: December, 1999 by PSLi; add in dt calculation in case
!                 of subcycle of artificial viscosity.
!
!  PURPOSE:  Computes the new timestep for explicit calculations from
!  the values of the field variables updated by the source and transport
!  steps.
!
!  In explicit calculations, the timestep is given by:
!
!     dt = courno * sqrt [ min ( dtcs**2 + dtv1**2 + dtv2**2 + dtv3**2
!                              + dtal**2 + dtqq**2 + dtrd**2 ) ]
!
!  where the variable names are described below.  The timestep can be
!  reduced in size by any amount, but can be larger than the old timstep
!  by no more than a factor of 1.26.
!
!  LOCAL VARIABLES:
!
!  i-sweep
!  dr*i      inverse of distance across zone in 1-, 2-, or 3-direction
!  drimax    maximum of dr1i, dr2i, and dr3i
!  dt**i2i   square of the inverse time step of the ** physical process
!            gathered during the i-sweep.  Possible values of ** are:
!                cs = sound speed
!                v1 = fluid motion in x1 direction
!                v2 = fluid motion in x2 direction
!                v3 = fluid motion in x2 direction
!                al = Alfven speed
!                qq = artificial viscosity
!                rd = radiation diffusion
!            The first five are vectors; the next to last is a scalar
!            which has been computed in ARTIFICIALVISC (passed in root).
!            Last is a scalar computed in rad_pfl (passed in radexp).
!  dttoi2i   sum of dt**i2i (without artificial viscosity contribution)
!
!  j-sweep
!  dttoi2j   vector of maximum values of dttoi2i from each i-sweep.
!            This vector is filled during a j-sweep.
!  imaxj     i-index of each dttoi2j
!  dt**i2j   values of above at zone i=imaxj
!
!  k-sweep
!  dttoi2k   vector of maximum values of dttoi2j from each j-sweep.
!  imaxk     i-index of each dttoi2k.
!  jmaxk     j-index of each dttoi2k.
!  dt**i2j   values of above at zone i=imaxk, j=jmaxk
!
!  grand maximum inverse time
!  imin      i-index where dttoi2k is a maximum
!  jmin      j-index where dttoi2k is a maximum
!  kmin      k-index where dttoi2k is a maximum
!  dt**      time step of the ** physical process at i=imax, j=jmax,
!            and k=kmax
!  dtnew     new timestep, which is limited to be no greater than 1.26
!            times the previous timestep.
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use field
      use bndry
      use root
      use scratch
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, &
                  imax, jmax, kmax, &
                  imin, jmin, kmin, &
                  k0  , k1  , k2  , j0
!
      real(rl) :: dtcsm, dtv1m, dtv2m, dtv3m, dtalm, dttoi2m
!
!-----------------------------------------------------------------------
!
! Find the minimum time step required by the Courant condition for
! this tile.  We will first compute 1/dt**2 required by each of the
! various physical processes in the calculation, and save the maximum
! value of their sum at each zone.  In the process, we will compute
! the updated velocity.
!
       dttoi2m = 0.0
       dtcsm   = 0.0
       dtv1m   = 0.0
       dtv2m   = 0.0
       dtv3m   = 0.0
       dtalm   = 0.0
       imin    = is
       jmin    = js
       kmin    = ks
       if (nx2z .eq. 1) then
!
! For one-dimensional problems only
!
         j0 = js
       else
         j0 = js + 1
       endif
       if (nx3z .eq. 1) then
!
! For one-dimensional or two-dimensional problems only
!
         k0 = ks
       else
         k0 = ks + 1
       endif
!
! Divide the computational volume into three equal pieces.  We must
! have at least 3 active zones in the 3-direction.
!
      if(ldimen .eq. 3) then
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
      else
       k1 = ks
       k2 = ks
      endif
!
      do i = 1, 6
       bvstat(i,1) = 0
       bvstat(i,3) = 0
       bvstat(i,4) = 0
       bvstat(i,5) = 0
       bvstat(i,7) = 0
       bvstat(i,8) = 0
      enddo
!......................................................................
!
! i boundaries
!
!    1) Post sends and receives. 
!
       nreq = 0
       nsub = nsub + 1
      if(xtotnrg) then
       imax = ie - 1
       if (nx2z .eq. 1) then
!
! One-dimensional problems only
!
         jmax = js
       else
         jmax = je - 1
       endif
       kmax = ke - 1
       call bvald  (1,1,0,0,0,0,d   )
       call bvalv1 (0,1,0,0,0,0,w3da)
      else ! xtotnrg
       imax = ie
!JH       jmax = je
!JH       kmax = ke
       if(ldimen .gt. 1) then
        jmax = je
       else
        jmax = js
       endif
       if(ldimen .gt. 2) then
        kmax = ke
       else
        kmax = ks
       endif
       call bvald  (1,0,0,0,0,0,d   )
      endif ! xtotnrg
!
!    2) Do first 1/3 of the interior points.
!
       call newdt (is+1,imax,j0  ,jmax,k0  ,k1 &
                      ,  imin,jmin,kmin &
                      ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
!
!       subroutine newdt (ibeg,iend,jbeg,jend,kbeg,kend
!     &                ,  imin,jmin,kmin
!     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
!
!
!    3) Wait for communications to complete.
!
#ifdef MPI_USED
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.
!
      if(ldimen .gt. 1) then
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d   )
        call bvalv2 (0,0,0,1,0,0,w3db)
       else ! xtotnrg
        call bvald  (0,0,1,0,0,0,d   )
       endif ! xtotnrg
      endif ! ldimen
!
!    2) Do middle 1/3 of the interior points, plus some on borders.
!
       call newdt (is  ,is  ,j0  ,jmax,k0  ,k1 &
                      ,  imin,jmin,kmin &
                      ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      if(xtotnrg) then
       call newdt (ie  ,ie  ,j0  ,jmax,k0  ,k1 &
                      ,  imin,jmin,kmin &
                      ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      endif ! xtotnrg
      if (ldimen .eq. 3) then
         call newdt (is  ,ie  ,j0  ,jmax,k1+1,k2 &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      endif
!
!    3) Wait for communications to complete.
!
#ifdef MPI_USED
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! k boundaries
!
!    1) Post sends and receives.
!
      if(ldimen .eq. 3) then
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,0,0,1,1,d   )
        call bvalv3 (0,0,0,0,0,1,w3dc)
       else ! xtotnrg
        call bvald  (0,0,0,0,1,0,d   )
       endif ! xtotnrg
      endif ! ldimen
!
!    2) Do last 1/3 of the interior points, plus some on borders.
!
!JH       if (nx2z .gt. 1) then
       if (ldimen .gt. 1) then
         call newdt (is  ,ie  ,js  ,js  ,k0  ,k2 &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
!        if(ldimen .eq. 2) then
!         call newdt (is  ,ie  ,je  ,je  ,k0  ,k2
!     &                  ,  imin,jmin,kmin
!     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
!        endif ! ldimen
        if(xtotnrg) then
         call newdt (is  ,ie  ,je  ,je  ,k0  ,k2 &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        endif ! xtotnrg
       endif
!JH       if (nx3z .gt. 1) then
       if (ldimen .gt. 2) then
         call newdt (is  ,ie  ,js  ,je  ,k2+1,kmax &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
       endif
!
!      Mark the velocity boundary values out of date.
!
       do 10 i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
10     continue
!
!    3) Wait for communications to complete. 
!
#ifdef MPI_USED
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! Finally, do the remaining border zones.
!
       if (ke .gt. ks) then
         call newdt (is  ,ie  ,js  ,je  ,ks  ,ks &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        if(xtotnrg) then
         call newdt (is  ,ie  ,js  ,je  ,ke  ,ke &
                        ,  imin,jmin,kmin &
                        ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        endif ! xtotnrg
       endif
!
!-----------------------------------------------------------------------
!
! Compute preliminary new time step.
!
#ifdef MPI_USED
!
! Now find the smallest dtnew among all tiles, and send the result
! to all in buf_out.  We need the info for the print statement below
! if the time step is too small, so use the MPI_MINLOC operation
! to pass values AND RANKS.
!
      if(xhydro) then
       buf_in(1) = dttoi2m
       buf_in(2) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                         , MPI_2DOUBLE_PRECISION &
                         , MPI_MAXLOC, comm3d, ierr)
       dttoi2m   =   buf_out(1)
!
       buf_in(1) = dtqqi2
       buf_in(3) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                         , MPI_2DOUBLE_PRECISION &
                         , MPI_MAXLOC, comm3d, ierr)
       dtqqi2    =   buf_out(1)
      endif
      if(lrad .ne. 0) then
       buf_in(1) = dtnri2
       buf_in(4) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                         , MPI_2DOUBLE_PRECISION &
                         , MPI_MAXLOC, comm3d, ierr)
       dtnri2    =   buf_out(1)
!
       buf_in(1) = dtimrdi2
       buf_in(5) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                         , MPI_2DOUBLE_PRECISION &
                         , MPI_MAXLOC, comm3d, ierr)
       dtimrdi2  =   buf_out(1)
      endif
#endif /* MPI */
!
      if(lrad .ne. 0) then
       if(xforce .eqv. .false.) then
        dtnew = courno / (sqrt( dtnri2 + dtimrdi2 ) + tiny )
       else ! xforce
        dtnew = courno / ( sqrt ( dttoi2m + dtqqi2 + dtnri2 + &
                                  dtimrdi2) + tiny )
       endif ! xforce
      endif ! lrad
!
      if(lrad .eq. 0) then
!PS
       if(xsubav) then
        dtnew  = courno / ( sqrt ( dttoi2m ) + tiny )
       else
        dtnew  = courno / ( sqrt ( dttoi2m + dtqqi2 ) + tiny )
       endif
!
      endif ! lrad
!
!#ifdef MPI_USED
!c
!c Now find the smallest dtnew among all tiles, and send the result
!c to all in buf_out.  We need the info for the print statement below 
!c if the time step is too small, so use the MPI_MINLOC operation
!c to pass values AND RANKS.
!c
!       buf_in(1) = dtnew
!       buf_in(2) = real( myid )
!c  M-MML: on advice of Jakob Pichlmeier (Cray Munich) changed MPI_2FLOAT
!c  to MPI_2DOUBLE_PRECISION...  26.2.98
!       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
!     &                   , MPI_2DOUBLE_PRECISION
!     &                   , MPI_MINLOC, comm3d, ierr)
!       dtnew  =   buf_out(1)
!#endif /* MPI */
       dt     =   min ( dtnew, 1.26*dt )
       if(time .ne. tlim) then
        if ((time+dt) .gt. tlim) dt = tlim-time
       endif
!
       if (dt .le. dtmin) then
!
! Determine which tile requires this short time step and report its
! coordinates and other info.  Convert 1/dt**2 to dt for each 
! physical process.
!
#ifdef MPI_USED
         if ( myid .eq. int( buf_out(2) ) ) then
!
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2010) coords(1), coords(2), coords(3) &
                         , imin, jmin, kmin, nhy, dt, dtmin, dtcs &
                         , dtv1, dtv2, dtv3, dtqq, dtal
         endif
         if ( myid .eq. int( buf_out(3) ) ) then
!
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2011) coords(1), coords(2), coords(3) &
                         , imin, jmin, kmin, nhy, dt, dtmin, dtcs &
                         , dtv1, dtv2, dtv3, dtqq, dtal
         endif
#endif /* MPI */
#ifndef MPI_USED
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2010) coords(1), coords(2), coords(3) &
                         , imin, jmin, kmin, nhy, dt, dtmin, dtcs &
                         , dtv1, dtv2, dtv3, dtqq, dtal
#endif /* NO MPI */
         nwarn = nwarn + 1
       endif
       return
!
!-----------------------------------------------------------------------
!----------------------- Write format statements -----------------------
!-----------------------------------------------------------------------
!
2010   format('NUDT    : **** VELOCITY **** Hot zone on tile (',i4,',' &
             ,i4,',',i4,')',/ &
             ,'NUDT    : at i=',i4,' j=',i4,' k=',i4,' (dt < dtmin)' &
             ,', nhy  = ',i6,',',/ &
             ,'NUDT    : dt   = ',1pe12.5,', dtmin= ',1e12.5,', dtcs = ' &
             ,1e12.5,',',/ &
             ,'NUDT    : dtv1 = ',1e12.5,', dtv2 = ',1e12.5,', dtv3 = ' &
             ,1e12.5,',',/ &
             ,'NUDT    : dtqq = ',1e12.5,', dtal = ' &
             ,1e12.5,'.')
!
2011   format('NUDT    : **** VISCOSITY **** Hot zone on tile (',i4,',' &
             ,i4,',',i4,')',/ &
             ,'NUDT    : at i=',i4,' j=',i4,' k=',i4,' (dt < dtmin)' &
             ,', nhy  = ',i6,',',/ &
             ,'NUDT    : dt   = ',1pe12.5,', dtmin= ',1e12.5,', dtcs = ' &
             ,1e12.5,',',/ &
             ,'NUDT    : dtv1 = ',1e12.5,', dtv2 = ',1e12.5,', dtv3 = ' &
             ,1e12.5,',',/ &
             ,'NUDT    : dtqq = ',1e12.5,', dtal = ' &
             ,1e12.5,'.')
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                  N U D T                  \\\\\\\\\\
!
!=======================================================================
!
!
