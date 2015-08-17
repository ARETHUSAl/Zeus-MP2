!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 A D V X 2                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine advx2 (dlo,den ,eod,edn ,mflx,s1,s2,s3,ero,ern ,abo,abn )
!
!    RAF, 2/17/97
!
!  PURPOSE: 
!    Controls the update of density, energy, and momenta
!    from the advection terms in the 2-direction.
!
!  INPUT:
!    dlo         Mass            density prior to update
!    eod         Specific energy density prior to update
!
!  OUTPUT:
!    den         Mass            density    after update
!    edn         Specific energy density    after update
!
!  I/O:
!    s1,s2,s3    Momentum density components (get updated)
!
!  LOCAL:
!    mflx        Mass flux in the 2-direction at zone faces 
!
!  EXTERNALS:
!    BVALV1  , BVALV2  , BVALV3
!    BVALD   , BVALE
!
!-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use root
      use field
      use bndry
      use grid
      use mpiyes
      use mpipar
!
      implicit none
!
      integer  :: k1, k2, i, j, k
      real(rl) :: p2
!
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn)
      real(rl) :: s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn)
      real(rl) :: eod(in,jn,kn), edn(in,jn,kn)
      real(rl), optional :: ero(in,jn,kn), ern(in,jn,kn)
      real(rl), optional :: abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
!
      real(rl) :: atwid (ijkn)
      real(rl) :: atwid1 (ijkn), atwid2 (ijkn), atwid3 (ijkn)
      real(rl) :: atwidj (ijkn)
!
!      Tunable data
!
      data p2 / 0.9 /
!
!-----------------------------------------------------------------------
!
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k1   = nint( ( real( ke + ks ) - p2 * nx3z ) * haf )
       k1   = max( k1, ks + 3)
       k2   = ke + ks - k1
       k2   = min( k2, ke - 1 )
!
!-----------------------------------------------------------------------
      if(ldimen .eq. 2) go to 222
!-----------------------------------------------------------------------
!
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
       call bvalv2 (1,0,0,0,0,0,v2 )
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
         call tranx2 (is+1,ie,js+3,je-2,ks+1,k1,dlo,den &
                   ,eod,edn ,mflx,atwid,ero,ern  ,abo,abn)
       else 
         call tranx2 (is+1,ie,js+3,je-2,ks+1,k1,dlo,den &
                   ,eod,edn ,mflx,atwid)
       endif
       call momx2  (is+2,ie,js+4,je-3,ks+2,k1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
!
! We need all the density slabs.
!
       call bvald  (0,0,3,3,0,0,dlo)
       if(nspec .gt. 1) call bvalabuns(0,0,3,3,0,0,abo)
!
!
! We need all slabs of eod.
!
       if(xiso .eqv. .false.) call bvale  (0,0,3,3,0,0,eod)
!
       if(lrad .ne. 0) call bvalert(0,0,3,3,0,0,ero)
!
! TRANX2 and MOMX2 together need all 3 velocities at js-2 through je+2.
!
       call bvalv1 (0,0,3,3,0,0,v1 )
       call bvalv2 (0,0,3,3,0,0,v2 )
       call bvalv3 (0,0,3,3,0,0,v3 )
!
!    2) Do middle portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,is  ,js+3,je-2,ks+1,k1,dlo,den &
                     ,eod,edn  ,mflx,atwid,ero,ern  ,abo,abn)
       else
         call tranx2 (is  ,is  ,js+3,je-2,ks+1,k1,dlo,den &
                     ,eod,edn ,mflx,atwid)
       endif
       call momx2  (is  ,is+1,js+4,je-3,ks+2,k1,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,js+3,je-2,k1+1,k2,dlo,den &
                     ,eod,edn ,mflx,atwid ,ero,ern  ,abo,abn)
       else 
         call tranx2 (is  ,ie  ,js+3,je-2,k1+1,k2,dlo,den &
                     ,eod,edn ,mflx,atwid)
       endif
       call momx2  (is  ,ie  ,js+4,je-3,k1+1,k2,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! k boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,0,0,1,0,dlo)
       call bvalv2 (0,0,0,0,1,0,v2 )
!
!    2) Do last portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,js  ,js+2,ks+1,k2,dlo,den &
                     ,eod,edn  ,mflx,atwid,ero,ern  ,abo,abn)
       else
         call tranx2 (is  ,ie  ,js  ,js+2,ks+1,k2,dlo,den &
                     ,eod,edn ,mflx,atwid)
       endif
       call momx2  (is  ,ie  ,js  ,js+3,ks+2,k2,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,je-1,je  ,ks+1,k2,dlo,den &
                     ,eod,edn ,mflx,atwid,ero,ern   ,abo,abn )
       else
         call tranx2 (is  ,ie  ,je-1,je  ,ks+1,k2,dlo,den &
                     ,eod,edn ,mflx,atwid) 
       endif
       call momx2  (is  ,ie  ,je-2,je  ,ks+2,k2,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,js  ,je  ,k2+1,ke,dlo,den &
                     ,eod,edn  ,mflx,atwid,ero,ern  ,abo,abn )
       else
         call tranx2 (is  ,ie  ,js  ,je  ,k2+1,ke,dlo,den &
                     ,eod,edn  ,mflx,atwid)
       endif
       call momx2  (is  ,ie  ,js  ,je  ,k2+1,ke,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
! Mark d and e/d (e) boundary values out of date.
!
       do 40 i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
40     continue
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border zones.
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den &
                     ,eod,edn  ,mflx,atwid,ero,ern  ,abo,abn)
       else
         call tranx2 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den &
                     ,eod,edn ,mflx,atwid)
       endif
!      write(*,"('ADVX2 before: s2 =',1pd13.5)")s2(3,3,3)
       call momx2  (is  ,ie  ,js  ,je  ,ks, ks+1, s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!      write(*,"('ADVX2 after: s2 =',1pd13.5)")s2(3,3,3)
!
      go to 999
!
!=======================================================================
!     2D TRANSPORT
!=======================================================================
!
222   continue
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
       call bvalv2 (1,0,0,0,0,0,v2 )
!
      if(xiso) then
       if (nseq .eq. 1) then
!
! We need to make a copy of the density, since we skipped pdv.
!
           do j=js-2,je+2
             do i=is-2,ie+2
               dlo(i,j,ks) = den(i,j,ks)
             enddo
           enddo
       endif
!
      endif ! xiso
       call bvald  (1,0,0,0,0,0,dlo)
!
!    2) Do first portion of the interior points.
!
       if (lrad .ne. 0) then
         call tranx2 (is+1,ie,js+3,je-2,ks,ks,dlo,den &
                     ,eod,edn ,mflx,atwid,ero,ern  ,abo,abn )
       else
         call tranx2 (is+1,ie,js+3,je-2,ks,ks,dlo,den &
                     ,eod,edn ,mflx,atwid)
       endif
!
       call momx2  (is+2,ie,js+4,je-3,ks,ks,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! j boundaries
!
!    1) Post sends and receives.
!
       nreq = 0
       nsub = nsub + 1
!
! We need all the density slabs.
!
       call bvald  (0,0,3,3,0,0,dlo)
       if(nspec .gt. 1) call bvalabuns(0,0,3,3,0,0,abo)
!
!
! We need all slabs of eod.
!
       if(xiso .eqv. .false.) call bvale  (0,0,3,3,0,0,eod)
!
       if(lrad .ne. 0) call bvalert(0,0,3,3,0,0,ero)
!
! TRANX2 and MOMX2 together need all 3 velocities at js-2 through je+2.
!
       call bvalv1 (0,0,3,3,0,0,v1 )
       call bvalv2 (0,0,3,3,0,0,v2 )
       call bvalv3 (0,0,3,3,0,0,v3 )
!
!    2) Do middle portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,is  ,js+3,je-2,ks,ks,dlo,den &
                     ,eod,edn   ,mflx,atwid,ero,ern  ,abo,abn )
       else
         call tranx2 (is  ,is  ,js+3,je-2,ks,ks,dlo,den &
                     ,eod,edn  ,mflx,atwid)
       endif
!
       call momx2  (is  ,is+1,js+4,je-3,ks,ks,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border zones.
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,js  ,js+2,ks, ks, dlo,den &
                     ,eod,edn  ,mflx,atwid,ero,ern  ,abo,abn )
       else
         call tranx2 (is  ,ie  ,js  ,js+2,ks, ks, dlo,den &
                     ,eod,edn  ,mflx,atwid)
       endif
       call momx2  (is  ,ie  ,js  ,js+3,ks,ks,s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
       if (lrad .ne. 0) then
         call tranx2 (is  ,ie  ,je-1,je  ,ks,ks,dlo,den &
                     ,eod,edn   ,mflx,atwid,ero,ern  ,abo,abn )
       else
         call tranx2 (is  ,ie  ,je-1,je  ,ks,ks,dlo,den &
                      ,eod,edn  ,mflx,atwid)
       endif
!
       call momx2  (is  ,ie  ,je-2,je  ,ks, ks, s1,s2,s3,mflx, &
                    atwid1,atwid2,atwid3,atwidj)
!
! Mark d and e/d (e) boundary values out of date.
!
999    continue
!
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
       enddo
!
       return
end subroutine advx2
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 A D V X 2                 \\\\\\\\\!
!=======================================================================
!
!
