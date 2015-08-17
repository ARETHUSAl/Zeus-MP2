!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              L O R E N T Z _ D            \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine lorentz_d (w1, w2, w3, u1, u2, u3)
!
!  PURPOSE
!    Driver for lorentz update. Updates the velocities due to the
!    transverse magnetic forces (the longitudinal forces are handled 
!    in FORCES).  Adapted from FORCES_D; but LORENTZ needs two slabs 
!    of velocity both above and below
!
!  AUTHOR
!    Mordecai-Mark Mac Low
!
!  LAST MODIFIED
!    3 Mar 98 by M-M Mac Low
!
!    April 29, 2006 by John Hayes:  consolidated BVAL* and LORENTZ
!    calls so that all BVAL* updates are completed before calling
!    subroutine LORENTZ.  LORENTZ is therefore called one time with
!    all zones being processed as a group.  The original form of this
!    routine attempted to subdivide the LORENTZ update and perform
!    it in stages, alternating calls to LORENTZ with calls to the
!    BVAL* routines on each tile face.  Unfortunately, this does not work
!    properly in parallel.  It may be possible to correct the problem,
!    but as of this date it remains unsolved.  This
!    driver routine, in its simplified form, does produce proper
!    results (when compare to a 1-processor run) when used in a
!    parallel calculation.  -- JCH 04/29/2006
!
!    May 15, 2006 by John Hayes: changed velocity references to point
!    to local scratch arrays w[123] and u[123].  Only density is accessed
!    via module FIELD.  Calls to BVALV[123] and LORENTZ are affected.
!    These changes, in combination with those made in CT/CT_2D/CT_1D,
!    correct an error in which two variables were assigned the same
!    space in memory.
!
!  INPUT
!    w[1,2,3]    velocity components before acceleration
!
!  OUTPUT
!    u[1,2,3]    velocity components after update
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    BVALD   , BVALE  
!    LORENTZ
!
!.......................................................................
!
      use real_prec
      use config
      use param
      use grid
      use root
      use field, ONLY : d
      use bndry
      use scratch
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, k1, k2
      real(rl) :: u1(in,jn,kn), u2(in,jn,kn), u3(in,jn,kn)
      real(rl) :: w1(in,jn,kn), w2(in,jn,kn), w3(in,jn,kn)
!
      if(xforce .eqv. .false.) then
!
! Just copy the old velocities to the new for all active zones.
!
      do k=ks,ke
       do j=js,je
         do i=is,ie
          u1(i,j,k) = w1(i,j,k)
          u2(i,j,k) = w2(i,j,k)
          u3(i,j,k) = w3(i,j,k)
         enddo ! i
        enddo ! j
       enddo ! k
       return
      endif ! xforce
!-----------------------------------------------------------------------
!
! Divide the computational volume into three equal pieces.  We must
! have at least 5 active zones in the 3-direction.
!
      k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
      k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
!
!
!     asif: accomodate 2D MHD runs
!
      if(ldimen.eq.2) goto 222
!
!     JCH: accomodate 1D MHD runs
!
      if(ldimen.eq.1) goto 111
!......................................................................
!
! i boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
!
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
! j boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,3,1,0,0,d )
       call bvalv1 (0,0,3,3,0,0,w1)
       call bvalv2 (0,0,3,3,0,0,w2)
       call bvalv3 (0,0,3,3,0,0,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
! k boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,0,0,3,1,d )
       call bvalv1 (0,0,0,0,3,3,w1)
       call bvalv2 (0,0,0,0,3,3,w2)
       call bvalv3 (0,0,0,0,3,3,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
       call lorentz (is,ie,js,je,ks,ke,w1,w2,w3,u1,u2,u3)
!
       goto 999
222     continue
!	asif: 2D lorentz
!......................................................................
!
! i boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
! j boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,3,1,0,0,d )
       call bvalv1 (0,0,3,3,0,0,w1)
       call bvalv2 (0,0,3,3,0,0,w2)
       call bvalv3 (0,0,3,3,0,0,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
       call lorentz_2D(is,ie,js,je,ks,ks,w1,w2,w3,u1,u2,u3)
!
       goto 999
!
111     continue
!
!      1-D Lorentz driver (JCH)
!......................................................................
!
! i boundaries
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
       call lorentz_1D(is,ie,js,js,ks,ks,w1,w2,w3,u1,u2,u3)
!
999    continue
!
!      Mark the velocities as out of date.
!
      do  i = 1,6
       bvstat(i,1) = 0      !  d
       bvstat(i,3) = 0      !  w1
       bvstat(i,4) = 0      !  w2
       bvstat(i,5) = 0      !  w3
      enddo
!
      return
!
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////           L O R E N T Z _ D               \\\\\\\\\!
!=======================================================================
