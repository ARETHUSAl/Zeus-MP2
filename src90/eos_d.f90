!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 E O S _ D                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine eos_d
!
!     Driver for non-ideal equations of state (leos > 1).  The logic
!     is cloned from the forces_d routine written by R. A. Fiedler.
!     Note that the EOS must be computed in one layer of inner ghost
!     zones on each face so that pressure gradients at the inner
!     boundaries can be computed by FORCES.
!
!     Written by J. Hayes, 5-2003
!
      use real_prec
      use config
      use param
      use grid
      use root
      use field
      use bndry
      use scratch
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer :: i, j, k, k1, k2, imax, jmax, kmax
!
!-----------------------------------------------------------------------
!
      if(.not. xtotnrg) then
       imax = ie
       jmax = je
       kmax = ke
      else
       imax = ie-2
       jmax = je-2
       kmax = ke-2
      endif
      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     3D Grids
!-----------------------------------------------------------------------
!
! Divide the computational volume into three equal pieces.  We must
! have at least 5 active zones in the 3-direction.
!
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
!
!
!    2) Do first portion of the interior points.
!
!
! i boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
       call bvalv2(1,1,0,0,0,0,v2)
       call bvalv3(1,1,0,0,0,0,v3)
      endif
      if(leos .ne. 1) call bvalt(1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
!
      call eos(is+1,imax,js+1,jmax,ks+1,k1)
!
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,1,0,0,0,d )
      call bvale(0,0,1,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,1,1,0,0,v1)
       call bvalv2(0,0,1,1,0,0,v2)
       call bvalv3(0,0,1,1,0,0,v3)
      endif
      if(leos .ne. 1) call bvalt(0,0,1,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,1,0,0,0,abun)
!
!    2) Do middle portion of the interior points, plus some on borders.
!
      call eos(is  ,is  ,js+1,jmax,ks+1,k1)
      call eos(is  ,ie  ,js+1,jmax,k1+1,k2)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! k boundaries
!
!    1) Post sends and receives.  We need d(k+1) for ZRP or RTP.
!
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,0,0,1,0,d )
      call bvale(0,0,0,0,1,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,0,0,1,1,v1)
       call bvalv2(0,0,0,0,1,1,v2)
       call bvalv3(0,0,0,0,1,1,v3)
      endif
      if(leos .ne. 1) call bvalt(0,0,0,0,1,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,0,0,1,0,abun)
!
!    2) Do last portion of the interior points, plus some on borders.
!
      call eos(is  ,ie  ,js  ,js  ,ks+1,k2)
      call eos(is  ,ie  ,js  ,je  ,k2+1,kmax)
!
!      Mark the boundary values out of date.
!
      do 10 i = 1,6
       bvstat(i,7) = 0      !  v1
       bvstat(i,8) = 0      !  v2
10    continue
!
!    3) Wait for communications to complete. 
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border zones.
!
      call eos(is-1,is  ,js  ,je  ,ks  ,ke) ! 1-face + ghost layer
      call eos(is  ,ie  ,js-1,js  ,ks  ,ke) ! 2-face + ghost layer
      call eos(is  ,ie  ,js  ,je  ,ks-1,ks) ! 3-face + ghost layer
      if(xtotnrg) then
       call eos(ie-1,ie  ,js  ,je  ,ks  ,ke-2)
       call eos(is  ,ie-2,je-1,je  ,ks  ,ke-2)
       call eos(is  ,ie  ,js  ,je  ,ke-1,ke  )
      endif
      go to 999
!
!-----------------------------------------------------------------------
!     2D Grids
!-----------------------------------------------------------------------
!
222   continue
      if(xtotnrg) then
       do i = 1, 6
        bvstat(i,3) = 0
        bvstat(i,4) = 0
       enddo
      endif
!
! i boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
       call bvalv2(1,1,0,0,0,0,v2)
      endif
      if(leos .ne. 1) call bvalt(1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
!
!    2) Do first portion of the interior points.
!
      call eos(is+1,imax,js+1,jmax,ks,ks)
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,1,0,0,0,d )
      call bvale(0,0,1,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,1,1,0,0,v1)
       call bvalv2(0,0,1,1,0,0,v2)
      endif
      if(leos .ne. 1) call bvalt    (0,0,1,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,1,0,0,0,abun)
!
      call eos(is  ,is  ,js+1,jmax,ks,ks)
      if(xtotnrg) call eos(ie-1,ie  ,js+1,jmax,ks,ks)
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
      call eos(is  ,ie  ,js  ,js  ,ks,ks)
!
      if(xtotnrg) call eos(is  ,ie  ,je-1,je  ,ks,ks)
!
      call eos(is-1,is-1,js  ,je  ,ks,ks) ! 1-dir inner ghost layer
      call eos(is  ,ie  ,js-1,js-1,ks,ks) ! 2-dir inner ghost layer
!
      if(.false.) then
       call eos(ie+1,ie+1,js  ,je  ,ks,ks) ! 1-dir outer ghost layer
       call eos(is  ,ie  ,je+1,je+1,ks,ks) ! 2-dir outer ghost layer
      endif
!
!      Mark the boundary values out of date.
!
       do i = 1,6
         bvstat(i,7) = 0      !  T
         bvstat(i,8) = 0      !  abun
       enddo
!
      go to 999
!-----------------------------------------------------------------------
!     1D Grids
!-----------------------------------------------------------------------
!
111   continue
!
! i boundaries
!
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
      endif
      if(leos .ne. 1) call bvalt    (1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
!
!    2) Do first portion of the interior points.
!
      call eos(is+1,imax,js,js,ks,ks)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
!     do first i zone plus inner ghost zone
!
      call eos(is-1,is,js,js,ks,ks)
      if(xtotnrg) call eos(ie-1,ie,js,js,ks,ks)
!
!      Mark the boundary values out of date.
!
      do i = 1,6
       bvstat(i,7) = 0 
       bvstat(i,8) = 0 
      enddo
!
999   return
      end
