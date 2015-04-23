!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               A V I S C _ D
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine avisc_d &
                  (v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!  PURPOSE
!    Driver for artificial viscosity module.  Updates velocity
!    components and material energy.
!
!  AUTHOR
!    Robert A. Fiedler
!
!  LAST MODIFIED
!    01/21/97
!
!  INPUT
!    v*old     Velocity components before viscous update.
!
!  OUPUT
!    v*new     Velocity components after  viscous update.
!    s*        Momentum components for use in the transport step.
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    BVALD
!    AVISC
!
!.......................................................................
!
      use real_prec
      use config
      use param
      use root
      use field
      use bndry
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
!
      real(rl) :: v1old(in,jn,kn),v2old(in,jn,kn),v3old(in,jn,kn), &
                  v1new(in,jn,kn),v2new(in,jn,kn),v3new(in,jn,kn), &
                  s1   (in,jn,kn),s2   (in,jn,kn),s3   (in,jn,kn)
!
      integer  :: i, k1, k2, imax, jmax, kmax
      real(rl) :: dvdxmn
!
!-----------------------------------------------------------------------
!
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
!
! Save the largest velocity gradient for computing the viscous
! time step.
!
      dvdxmn = 0.0
!
      nsub = nsub + 1
      if(xtotnrg) then
       call bvald  (1,1,0,0,0,0,d    )
       call bvalv1 (3,3,0,0,0,0,v1old)
      else
       call bvald  (1,0,0,0,0,0,d    )
       call bvalv1 (1,1,0,0,0,0,v1old)
      endif
!
      call avisc (is+1,ie-2,js+1,je-2,ks+1,k1  ,dvdxmn &
                 ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
#ifdef MPI_USED
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d    )
        call bvalv2 (0,0,3,3,0,0,v2old)
       else
        call bvald  (0,0,1,0,0,0,d    )
        call bvalv2 (0,0,1,1,0,0,v2old)
       endif
!
       call avisc (is  ,is  ,js+1,je-2,ks+1,k1  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (ie-1,ie  ,js+1,je-2,ks+1,k1  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js+1,je-2,k1+1,k2  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald (0,0,0,0,1,1,d    )
        call bvalv3(0,0,0,0,3,3,v3old)
       else
        call bvald (0,0,0,0,1,0,d    )
        call bvalv3(0,0,0,0,1,1,v3old)
       endif
!
       call avisc (is  ,ie  ,js  ,js  ,ks+1,k2  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,je-1,je  ,ks+1,k2  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js  ,je  ,k2+1,ke-2,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!      Mark the boundary values out of date.
!
       do 20 i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
20     continue
#ifdef MPI_USED
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
       call avisc (is  ,ie  ,js  ,je  ,ks  ,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js  ,je  ,ke-1,ke  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
      go to 999
!-----------------------------------------------------------------------
222   continue
!
! Save the largest velocity gradient for computing the viscous
! time step.
!
       dvdxmn = 0.0
!
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (1,1,0,0,0,0,d    )
        call bvalv1 (3,3,0,0,0,0,v1old)
       else
        call bvald  (1,0,0,0,0,0,d    )
        call bvalv1 (1,1,0,0,0,0,v1old)
       endif
!
       call avisc (is+1,ie-2,js+1,je-2,ks,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d    )
        call bvalv2 (0,0,3,3,0,0,v2old)
       else
        call bvald  (0,0,1,0,0,0,d    )
        call bvalv2 (0,0,1,1,0,0,v2old)
       endif
!
       call avisc (is  ,is  ,js+1,je-2,ks,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (ie-1,ie  ,js+1,je-2,ks,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!      Mark the boundary values out of date.
!
       do i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
       enddo
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!......................................................................
       call avisc (is  ,ie  ,js  ,js  ,ks  ,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,je-1,je  ,ks  ,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
      go to 999
!-----------------------------------------------------------------------
111   continue
!
! Save the largest velocity gradient for computing the viscous
! time step.
!
       dvdxmn = 0.0
!
! i boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (1,1,0,0,0,0,d    )
        call bvalv1 (3,3,0,0,0,0,v1old)
       else
        call bvald  (1,0,0,0,0,0,d    )
        call bvalv1 (1,1,0,0,0,0,v1old)
       endif
!
!    2) Do first portion of the interior points.
!
       call avisc (is+1,is+3,js,js,ks,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!       subroutine avisc (ibeg,iend,jbeg,jend,kbeg,kend,dvdxmn
!     &                  ,w1,w2,w3,u1,u2,u3,s1,s2,s3)
!
!c    w1,w2,w3   velocity values prior to viscosity update.
!c    u1,u2,u3   velocity values after    viscosity update.
!c    s1,s2,s3   updated momentum densities for transport step.
!
!    3) Wait for communications to complete.
!
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
!
!      Mark the boundary values out of date.
!
       do i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
       enddo
!......................................................................
!
! Do the remaining border zones.
!
       call avisc (is  ,is  ,js  ,js  ,ks  ,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is+4,ie  ,js  ,js  ,ks  ,ks  ,dvdxmn &
                  ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!......................................................................
999   continue
!......................................................................
!
!  Compute the viscous timestep.  Note that the minimum dv/dx is found
!  since it is less than 0.  Thus the minimum dv/dx gives the maximum
!  absolute value.  We'll do a global min on this later, in nudt.
!
       dtqqi2 = ( 4.0 * qcon * dvdxmn )**2
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               A V I S C _ D               \\\\\\\\\\
!
!=======================================================================
!
