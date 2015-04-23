!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              F O R C E S _ D              \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine forces_d &
       (v1old, v2old, v3old, v1new, v2new, v3new)
!
!  PURPOSE
!    Driver for forces module.  Updates the velocities due to the
!    following body forces:
!
!    1) Thermal pressure gradient
!    2) Self-gravity
!    3) Rotational Pseudo-Forces
!    4) Gravitational point mass
!    5) Radiation force
!    6) Magnetic pressure
!
!  AUTHOR
!    Robert A. Fiedler
!
!  LAST MODIFIED
!    02/17/97, MHD on 3 Mar 98 by M-M Mac Low
!
!  INPUT
!    v*old    velocity components before acceleration
!
!  OUTPUT
!    v*new    velocity components after update
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    BVALD   , BVALE   , BVALER
!    FORCES
!
!.......................................................................
!
      use real_prec
      use config
      use param
      use grid
      use root
      use field
      use bndry
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
      real(rl) :: v1old(in,jn,kn), v2old(in,jn,kn), v3old(in,jn,kn), &
                  v1new(in,jn,kn), v2new(in,jn,kn), v3new(in,jn,kn)
!
      real(rl) :: tg
      integer  :: sf
!
      integer :: i, j, k, k1, k2, imax, jmax, kmax
!
!-----------------------------------------------------------------------
!      Not doing forces, so bail
!-----------------------------------------------------------------------
!
      if(xforce .eqv. .false.) then
!
! Just copy the old velocities to the new for all active zones.
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             v1new(i,j,k) = v1old(i,j,k)
             v2new(i,j,k) = v2old(i,j,k)
             v3new(i,j,k) = v3old(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
!
       do i = 1, 6
        bvstat(i,3) = 0
        bvstat(i,4) = 0
        bvstat(i,5) = 0
       enddo
!
       return
      endif ! xforce
!
!-----------------------------------------------------------------------
      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111
!-----------------------------------------------------------------------
!
! Divide the computational volume into three equal pieces.  We must
! have at least 5 active zones in the 3-direction.
!
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
!
!......................................................................
!
! i boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald  (1,0,0,0,0,0,d )
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
      if(lgeom .eq. 3) then
       jmax = je - 1
       call bvalv2 (1,0,0,0,0,0,v2old)
       call bvalv3 (1,0,0,0,0,0,v3old)
      else
       jmax = je
      endif ! RTP
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
!
!    2) Do first portion of the interior points.
!
      call forces (is+1,imax,js+1,jmax,ks+1,k1 &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
#ifdef MPI_USED
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.  We need d(j+1) only for RTP.
!
      nreq = 0
      nsub = nsub + 1
!
      if(lgeom .eq. 3) then
       call bvald  (0,0,1,1,0,0,d )
      else
       call bvald  (0,0,1,0,0,0,d )
      endif ! lgeom
!
      if(xiso .eqv. .false.) call bvale  (0,0,1,0,0,0,e )
      if(lgeom .eq. 3) then
       call bvalv2 (0,0,0,1,0,0,v2old)
      endif
      if(lgeom .ne. 1) then
       call bvalv3 (0,0,1,0,0,0,v3old)
      endif
      if(lrad .ne. 0) call bvalers (0,0,1,0,0,0,er)
!
!    2) Do middle portion of the interior points, plus some on borders.
!
      call forces (is  ,is  ,js+1,jmax,ks+1,k1 &
                  ,v1old,v2old,v3old,v1new,v2new,v3new)
      call forces (is  ,ie  ,js+1,jmax,k1+1,k2 &
                  ,v1old,v2old,v3old,v1new,v2new,v3new)
#ifdef MPI_USED
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! k boundaries
!
!    1) Post sends and receives.  We need d(k+1) for ZRP or RTP.
!
       nreq = 0
       nsub = nsub + 1
!
      if(lgeom .ne. 1) then
       call bvald(0,0,0,0,1,1,d )
      else
       call bvald(0,0,0,0,1,0,d )
      endif
!
      if(xiso .eqv. .false.) call bvale  (0,0,0,0,1,0,e )
      if(lgeom .ne. 1) then
       kmax = ke - 1
       call bvalv3 (0,0,0,0,1,1,v3old)
      else ! lgeom
       kmax = ke
      endif ! lgeom
      if(lrad .ne. 0) call bvalers (0,0,0,0,1,0,er)
!
!    2) Do last portion of the interior points, plus some on borders.
!
      call forces (is  ,ie  ,js  ,js  ,ks+1,k2 &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
      if(lgeom .eq. 3) then
       call forces (is  ,ie  ,je  ,je  ,ks+1,k2 &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! rtp
!
      call forces (is  ,ie  ,js  ,je  ,k2+1,kmax &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
!      Mark the velocity boundary values out of date.
!
       do 10 i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
10     continue
#ifdef MPI_USED
!
!    3) Wait for communications to complete. 
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! Finally, do the remaining border zones.
!
      call forces (is  ,ie  ,js  ,je  ,ks  ,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
      if(lgeom .ne. 1) then
       call forces (is  ,ie  ,js  ,je  ,ke  ,ke &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! xtotnrg or zrp/rtp
!
      go to 999
!
!......................................................................
!     2D TRANSPORT
!......................................................................
222   continue
!
! i boundaries
!
!    1) Post sends and receives.  We never use d(i+1) in any of the
!       source term substeps, so we need to pass only 1 "m" layer
!       of density boundary values.  We also need only one
!       "m" layer of energy values.  If we need the velocity, we
!       need just one layer, but both "m" and "p".  This should
!       minimize communication for any physics/geometry.  By
!       exchanging the i, j, and k boundary data in three separate 
!       steps, we ensure that the corner and edge data are correctly
!       updated.
!
      nreq = 0
      nsub = nsub + 1
      if(leos .ne. 2) call bvald  (1,0,0,0,0,0,d )
!
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
!
      if(lgeom .eq. 3) then
       jmax = je - 1
       call bvalv2 (1,0,0,0,0,0,v2old)
       call bvalv3 (1,0,0,0,0,0,v3old)
      else
       jmax = je
      endif ! rtp
!
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
!
!    2) Do first portion of the interior points.
!
      call forces (is+1,imax,js+1,jmax,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.  We need d(j+1) only for RTP.
!
      nreq = 0
      nsub = nsub + 1
!
      if(lgeom .eq. 3) then
       call bvald  (0,0,1,1,0,0,d )
      else
       call bvald  (0,0,1,0,0,0,d )
      endif ! lgeom
      if(xiso .eqv. .false.) call bvale  (0,0,1,0,0,0,e )
!
      if(lgeom .eq. 3) then
       call bvalv2 (0,0,0,1,0,0,v2old)
      endif
!
      if(lgeom .ne. 1) then
       call bvalv3 (0,0,1,0,0,0,v3old)
      endif
!
      if(lrad .ne. 0) call bvalers (0,0,1,0,0,0,er)
!
!    2) Do middle portion of the interior points, plus some on borders.
!
      call forces (is  ,is  ,js+1,jmax,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
!
!    3) Wait for communications to complete.
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!
!    2) Do last portion of the interior points, plus some on borders.
!
      call forces (is  ,ie  ,js  ,js  ,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
      if(lgeom .eq. 3) then
       call forces (is  ,ie  ,je  ,je  ,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! rtp
!
!      Mark the velocity boundary values out of date.
!
       do i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v2
       enddo
!
      go to 999
!......................................................................
!     1D TRANSPORT
!......................................................................
111   continue
!
! i boundaries
!
!    1) Post sends and receives.  We never use d(i+1) in any of the
!       source term substeps, so we need to pass only 1 "m" layer
!       of density boundary values.  We also need only one
!       "m" layer of energy values.  If we need the velocity, we
!       need just one layer, but both "m" and "p".  This should
!       minimize communication for any physics/geometry.  By
!       exchanging the i, j, and k boundary data in three separate 
!       steps, we ensure that the corner and edge data are correctly
!       updated.
!
      nreq = 0
      nsub = nsub + 1
      call bvald  (1,0,0,0,0,0,d )
!
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
!
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
!
!    2) Do first portion of the interior points.
!
      call forces (is+1,imax,js,js,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
!       subroutine forces (ibeg,iend,jbeg,jend,kbeg,kend
!     &                   ,u1,u2,u3,w1,w2,w3)
!c
!c Arrays u1 , u2 , u3  hold the old velocity values, while
!c        w1 , w2 , w3  receive the updated values.
!
!
!    3) Wait for communications to complete.
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
!
      call forces (is  ,is  ,js,js,ks,ks &
                   ,v1old,v2old,v3old,v1new,v2new,v3new)
!
!      Mark the velocity boundary values out of date.
!
       do i = 1,6
         bvstat(i,3) = 0      !  v1
       enddo
!
999   return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////             F O R C E S _ D               \\\\\\\\\\
!
!=======================================================================
!
