!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 A D V X 3                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine advx3 (dlo,den ,eod,edn ,mflx,s1,s2,s3,ero,ern ,abo,abn )
!
!    RAF, 2/19/97
!
!  PURPOSE: 
!    Controls the update of density, energy, and momenta
!    from the advection terms in the 3-direction.
!
!  INPUT:
!    dlo         Mass            density prior to update
!    eod         Specific energy density prior to update
!    ero         Specific radiation energy density prior to update
!
!  OUTPUT:
!    den         Mass            density    after update
!    edn         Specific energy density    after update
!    ern         Specific radiation energy density    after update
!
!  I/O:
!    s1,s2,s3    Momentum density components (get updated)
!
!  LOCAL:
!    mflx        Mass flux in the 3-direction at zone faces 
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    BVALD   , BVALE   , BVALER
!
!-----------------------------------------------------------------------
      use real_prec
      use param
      use config
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
      integer  :: k1,k2,i
      integer  :: kbeg, kend, krange, kblocks, kskip, ktlb, kblk, kpage
      integer  :: j, k
      real(rl) :: p3
!
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn)
      real(rl) :: s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn)
      real(rl) :: eod(in,jn,kn), edn(in,jn,kn)
      real(rl), optional :: ero(in,jn,kn), ern(in,jn,kn)
      real(rl), optional :: abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
!
      real(rl) :: atwid (ijkn)
      real(rl) :: mflux (ijkn,1)
      real(rl) :: dtwid (ijkn,1), dd  (ijkn,1)
      real(rl) :: etwid (ijkn,1), deod(ijkn,1)
!
      real(rl) :: atwid1 (ijkn), atwid2 (ijkn), atwid3 (ijkn)
      real(rl) :: atwidj1(ijkn), atwidj2(ijkn), atwidj3(ijkn)
!
      real(rl) :: sflx  (ijkn,1), dq   (ijkn,1)
!
!      Tunable data
!
      data p3 / 0.9 /   ! Fraction of interior points to do 3rd stage
!
! Set blocking factor for 3rd stage k-loops.  The number of iterations
! in a block should be small enough so that all the data fits on the
! number of pages that the TLB can hold.  Thus,
!
! iterations = <TLB entries * page size> / <2-D arrays * word size>
!
! In tranx3 and momx3, 11 different 2-D arrays of data are used.  
! Assume 16kB pages and 8B words.  The number of iterations per
! block is forced to be at least 5.
!
!#ifdef TLB
!       data ktlb / TLB /  ! Get TLB size in MB from cpp -DTLB=n
       data ktlb / 128 /  ! This is the default for compilation on 
!                           SGI O2K
!#else /* TLB */
!       data ktlb / 9999 /    ! Assume many TLB entries
!#endif /* TLB */
!#ifdef BLK_MIN
!       data kblk / BLK_MIN /  ! Get min block size from cpp -DBLK_MIN=n
!#else /* BLK_MIN */
       data kblk / 5 /        ! Default min block size is this many
!#endif /* BLK_MIN */
!#ifdef PAGE
!       data kpage / PAGE /    ! Get PAGE size in bytes from cpp -DPAGE=n
!#else /* PAGE */
       data kpage / 16384 /   ! Assume this many bytes per page
!#endif /* PAGE */
!
!-----------------------------------------------------------------------
!
! Divide up the work.  Since we must have
!
!   ks+4 < k1-1
!   k1   < k2-1
!   k2   < ke-3
!
!   ke - ks .ge. 12   --- this is the smallest allowable k range.
!
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k2   = ke - int( p3 * nx3z )
       k2   = min( k2, ke - 3 )
       k2   = max( k2, ks + 6 )
       k1   = ( k2 + ks ) / 2
       k1   = max( k1, ks + 5 )
!......................................................................
!
! i boundaries
!
!    1) Post sends and receives. 
!       By exchanging the i, j, and k boundary data in three separate 
!       steps, we ensure that the corner and edge data are correctly
!       updated.
!
       nreq = 0
       nsub = nsub + 1
       call bvalv3 (1,0,0,0,0,0,v3 )
!
      if(xiso) then
       if (nseq .eq. 1) then
!
! We need to make a copy of the density, since we skipped pdv.
!
         do 30 k=ks-2,ke+2
           do 20 j=js-2,je+2
             do 10 i=is-2,ie+2
               dlo(i,j,k) = den(i,j,k)
10           continue
20         continue
30       continue
       endif
!
      endif ! xiso
       call bvald  (1,0,0,0,0,0,dlo)
!
!    2) Do first portion of the interior points.
!
       if (lrad .ne. 0) then
         call tranx3 (is+1,ie,js+1,je,ks+3,k1,dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero ,ern  ,abo,abn ) 
       else
         call tranx3 (is+1,ie,js+1,je,ks+3,k1,dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is+2,ie,js+2,je,ks+4,k1-1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
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
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,1,0,0,0,dlo)
       call bvalv3 (0,0,1,0,0,0,v3 )
!
!    2) Do second portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx3 (is  ,is  ,js+1,je,ks+3,k1,dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero ,ern  ,abo,abn ) 
       else
         call tranx3 (is  ,is  ,js+1,je,ks+3,k1,dlo,den ,eod,edn  &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is  ,is+1,js+2,je,ks+4,k1-1,s1,s2,s3,mflx,  &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx3 (is  ,ie  ,js+1,je,k1+1,k2,dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero ,ern  ,abo,abn ) 
       else
         call tranx3 (is  ,ie  ,js+1,je,k1+1,k2,dlo,den ,eod,edn  &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is  ,ie  ,js+2,je,k1  ,k2-1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
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
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
!
! We need all the density slabs.
!
       call bvald  (0,0,0,0,3,3,dlo)
       if(nspec .gt. 1) call bvalabuns(0,0,0,0,3,3,abo)
!
!
! We need all slabs of eod.
!
       if(xiso .eqv. .false.) call bvale  (0,0,0,0,3,3,eod)
!
       if(lrad .ne. 0) call bvalert(0,0,0,0,3,3,ero)
!
! TRANX3 and MOMX3 together need all 3 velocities at ks-2 through ke+2.
!
       call bvalv1 (0,0,0,0,3,3,v1 )
       call bvalv2 (0,0,0,0,3,3,v2 )
       call bvalv3 (0,0,0,0,3,3,v3 )
!
!    2) Do last portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx3 (is  ,ie  ,js  ,js  ,ks+3,k2,dlo,den ,eod,edn  &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero ,ern  ,abo,abn ) 
       else
         call tranx3 (is  ,ie  ,js  ,js  ,ks+3,k2,dlo,den ,eod,edn  &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is  ,ie  ,js  ,js+1,ks+4,k2-1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
!
! Block the k loop to reduce TLB misses; the k ranges for the earlier
! stages above should be small enough already.
!
       kskip   = max( ktlb * kpage / (11*8*in*jn), kblk)
       krange = ke-2 - (k2+1) + 1
       kblocks = krange / kskip
       if ( mod(krange,kskip) .eq. 0) kblocks = max(kblocks-1,0)
!
       do 35 kbeg = k2+1, k2+1 + kblocks*kskip, kskip 
         kend = min( kbeg + kskip - 1, ke-2 )
         if (lrad .ne. 0) then
           call tranx3 (is  ,ie  ,js  ,je  ,kbeg,kend,dlo,den ,eod,edn &
                       ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                       ,ero ,ern  ,abo,abn )
         else 
           call tranx3 (is  ,ie  ,js  ,je  ,kbeg,kend,dlo,den ,eod,edn &
                       ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
         endif
         call momx3  (is  ,ie  ,js  ,je  ,kbeg-1,kend-1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
35     continue
!
! Mark d and e/d (e) boundary values out of date.
!
       do 40 i=1,6
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
40     continue
!
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
       if (lrad .ne. 0) then
         call tranx3 (is  ,ie  ,js  ,je  ,ks, ks+2, dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero,ern  ,abo,abn ) 
       else
         call tranx3 (is  ,ie  ,js  ,je  ,ks, ks+2, dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is  ,ie  ,js  ,je  ,ks, ks+3, s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx3 (is  ,ie  ,js  ,je  ,ke-1, ke, dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod &
                     ,ero,ern  ,abo,abn )
       else
         call tranx3 (is  ,ie  ,js  ,je  ,ke-1, ke, dlo,den ,eod,edn &
                     ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
       endif
       call momx3  (is  ,ie  ,js  ,je  ,ke-2, ke, s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3, &
                    sflx,dq)
!
       return
end subroutine advx3
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 A D V X 3                 \\\\\\\\\\
!
!=======================================================================
!
!
