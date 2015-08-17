!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 P D V _ D                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine pdv_d &
                  (dcopy, eod)
!
!  PURPOSE
!    Compute the compressional work term (PdV) in the material
!    energy equation (for no radiation).
!
!  AUTHOR
!    Robert A. Fiedler
!
!  LAST MODIFIED
!    01/21/97
!
!  INPUT
!    (none)
!
!  OUTPUT
!    dcopy     A copy of the density.
!    eod       e/d for the transport step.
!
!  We need just 1 "p" layer of updated boundary data for v1, v2, and v3,
!  but none for d and e.  
!
!  Routine pdv also saves the density and e/d in (dcopy,eod) for use
!  in the transport step.  This is why it is being called even when
!  TOTAL_ENERGY is defined.
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    PDV
!.......................................................................
!
      use real_prec
      use config
      use param
      use root
      use field
      use bndry
      use grid
      use mpiyes
      use mpipar
!
      implicit NONE
!
      real(rl) :: dcopy(in,jn,kn), eod(in,jn,kn)
!
      integer  :: i, k1, k2, imax, jmax, kmax
      integer  :: j, k
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
! i boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)
!
!    2) Do first portion of the interior points.
!
!      Copy density boundary values.
!
       do k=ks-2,ks-1
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
!
       call pdv (is,ie-1,js,je-1,ks  ,k1 &
                ,dcopy ,eod )
!
!       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
!     &                ,dlo,eod)
!c
!c    dlo   Mass density (copy)                      for transport step.
!c    eod   Specific energy density e/d (or (e+p)/d] for transport step.
!
!
!      Copy density boundary values.
!
       do k=ks,ke
         do j=js-2,js-1
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
!
!    3) Wait for communications to complete.
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv2 (0,0,0,1,0,0,v2)
!
!    2) Do middle 1/3 of the interior points, and some on borders.
!
       call pdv (ie  ,ie  ,js,je-1,ks  ,k1 &
                ,dcopy ,eod )
       call pdv (is  ,ie  ,js,je-1,k1+1,k2 &
                ,dcopy ,eod )
!
!      Copy density boundary values.
!
       do k=ks,ke
         do j=je+1,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
       do k=ks,ke
         do j=js,je
           do i=is-2,is-1
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
!
!    3) Wait for communications to complete.
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! k boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv3 (0,0,0,0,0,1,v3)
!
!    2) Do last 1/3 of the interior points, and some on borders.
!
       call pdv (is  ,ie  ,je  ,je  ,ks  ,k2 &
                ,dcopy ,eod )
       call pdv (is  ,ie  ,js  ,je  ,k2+1,ke-1 &
                ,dcopy ,eod )
!
!      Copy density boundary values.
!
       do k=ke+1,ke+2
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
       do k=ks,ke
         do j=js,je
           do i=ie+1,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
!
!    3) Wait for communications to complete. 
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border.
!
       call pdv (is  ,ie  ,js  ,je  ,ke  ,ke &
                ,dcopy ,eod )
!
      go to 999
!======================================================================
!     2D TRANSPORT
!======================================================================
!
222   continue
!
! i boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)
!
!    2) Do first portion of the interior points.
!
!      Copy density values.
!
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,ks) = d(i,j,ks)
           enddo ! i
         enddo ! j
!
       call pdv (is,ie-1,js,je-1,ks  ,ks &
                ,dcopy ,eod )
!
!       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
!     &                ,dlo,eod)
!c
!c    dlo   Mass density (copy)                      for transport step.
!c    eod   Specific energy density e/d (or (e+p)/d] for transport step.
!
!    3) Wait for communications to complete.
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv2 (0,0,0,1,0,0,v2)
!
       call pdv (ie  ,ie  ,js,je-1,ks  ,ks &
                ,dcopy ,eod )
!
!    3) Wait for communications to complete.
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!
! Finally, do the remaining border.
!
       call pdv (is  ,ie  ,je  ,je  ,ks  ,ks &
                ,dcopy ,eod )
!
!      Mark the boundary values out of date.
!
      go to 999
!======================================================================
!     1D TRANSPORT
!======================================================================
!
111   continue
!
! i boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)
!
!    2) Do first portion of the interior points.
!
!      Copy density values.
!
         do j=js,js
           do i=is-2,ie+2
             dcopy(i,j,ks) = d(i,j,ks)
           enddo ! i
         enddo ! j
!
       call pdv (is,ie-1,js,js,ks  ,ks &
                ,dcopy ,eod )
!
!       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
!     &                ,dlo,eod)
!c
!c    dlo   Mass density (copy)                      for transport step.
!c    eod   Specific energy density e/d (or (e+p)/d] for transport step.
!
!    3) Wait for communications to complete.
!
       if(nreq .ne. 0) &
          call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border.
!
       call pdv (ie  ,ie  ,js  ,js  ,ks  ,ks &
                ,dcopy ,eod )
!
!      Mark the boundary values out of date.
!
999   continue
      do i = 1,6
       bvstat(i,2) = 0      !  e
       if(leos .ne. 1) bvstat(i,8) = 0 ! T
      enddo
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 P D V _ D                 \\\\\\\\\!
!=======================================================================
!
