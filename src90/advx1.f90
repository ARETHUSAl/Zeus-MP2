!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 A D V X 1                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!    
!=======================================================================
!
subroutine advx1 (dlo,den ,eod,edn ,mflx,s1,s2,s3, &
                  ero,ern ,abo,abn )
!
!    RAF, 2/17/97
!
!  PURPOSE: 
!    Controls the update of density, energy, and momenta
!    from the advection terms in the 1-direction.
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
!    mflx        Mass flux in the 1-direction at zone faces 
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
      use mpiyes
      use mpipar
!
      implicit none
!
      integer  :: k1, k2, i, j, k
!
      real(rl) :: p1
!
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn)
      real(rl) :: s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn)
      real(rl) :: eod(in,jn,kn), edn(in,jn,kn)
      real(rl), optional :: ero(in,jn,kn), ern(in,jn,kn)
      real(rl), optional :: abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
!
      real(rl) :: atwid (ijkn),  mflux (ijkn)
      real(rl) :: dtwid (ijkn),  dd    (ijkn)
      real(rl) :: etwid (ijkn),  deod  (ijkn)
!
      real(rl) :: atwid1(ijkn)
      real(rl) :: vtwid ( ijkn ), sflx  ( ijkn )
      real(rl) :: dq ( ijkn )
!
!      Tunable data
!
      data p1 / 0.9 /
!
!-----------------------------------------------------------------------
!
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k1   = int( real( nx3z ) * p1 ) + ks   !  The lion's share...
       k1   = min( k1, ke - 2 )               !  but not too much!
       k2   = ( ke + k1 ) / 2                 !  Half the remainder.
      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111
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
!
! TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
!
       call bvalv1 (3,3,0,0,0,0,v1 )
       call bvalv2 (3,3,0,0,0,0,v2 )
       call bvalv3 (3,3,0,0,0,0,v3 )
!
! We need all slabs of er/d.
!
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
!
! We need all slabs of eod.
!
       call bvale  (3,3,0,0,0,0,eod)
!
      else ! xiso
!
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
!
! We need all the density slabs.
!
       call bvald  (3,3,0,0,0,0,dlo)
       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)
!
!    2) Do first portion of the interior points.
!
       if (lrad .ne. 0) then
         call tranx1 (is+3,ie-2,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod , edn, mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero , ern, abo , abn )
       else
         call tranx1 (is+3,ie-2,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod , edn, mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is+4,ie-3,js+2,je  ,ks+2,k1,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
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
       call bvald  (0,0,1,0,0,0,dlo)
       call bvalv1 (0,0,1,0,0,0,v1 )
!
!    2) Do second portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,is+2,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo ,abn)
       else
         call tranx1 (is  ,is+2,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,is+3,js+2,je  ,ks+2,k1,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx1 (ie-1,ie  ,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern  ,abo,abn)
       else  
         call tranx1 (ie-1,ie  ,js+1,je  ,ks+1,k1,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (ie-2,ie  ,js+2,je  ,ks+2,k1,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,ie  ,js+1,je  ,k1+1,k2,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod,&
                      ero ,ern ,abo ,abn )
       else
         call tranx1 (is  ,ie  ,js+1,je  ,k1+1,k2,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,ie  ,js+2,je  ,k1+1,k2,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
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
       call bvalv1 (0,0,0,0,1,0,v1 )
!       call bvalv3 (0,0,0,0,0,1,v3 )     !   Did this already in pdv.
!
!    2) Do last portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,ie  ,js  ,js  ,ks+1,k2,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo,abn )
       else 
         call tranx1 (is  ,ie  ,js  ,js  ,ks+1,k2,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,ie  ,js  ,js+1,ks+2,k2,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,ie  ,js  ,je  ,k2+1,ke  ,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo ,abn )
       else 
         call tranx1 (is  ,ie  ,js  ,je  ,k2+1,ke  ,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,ie  ,js  ,je  ,k2+1,ke  ,s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
!
! Mark d and e/d (e) boundary values out of date.
!
       do 40 i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0
40     continue
!
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border zones.
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo,abn )
       else
         call tranx1 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
!
       call momx1  (is  ,ie  ,js  ,je  ,ks, ks+1, s1,s2,s3,mflx, &
                    atwid1,vtwid,sflx,dq)
!
       go to 999
!=======================================================================
!      2D TRANSPORT
!=======================================================================
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
!
! TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
!
      call bvalv1 (3,3,0,0,0,0,v1 )
      call bvalv2 (3,3,0,0,0,0,v2 )
      call bvalv3 (3,3,0,0,0,0,v3 )
!
! We need all slabs of er/d.
!
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
!
! We need all slabs of eod.
!
       call bvale  (3,3,0,0,0,0,eod)
!
      else ! xiso
!
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
!
! We need all the density slabs.
!
       call bvald  (3,3,0,0,0,0,dlo)
       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)
!
!    2) Do first portion of the interior points.
! 
       if (lrad .ne. 0) then
         call tranx1 (is+3,ie-2,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo ,abn )
       else
         call tranx1 (is+3,ie-2,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is+4,ie-3,js+2,je  ,ks,ks,s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
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
       call bvald  (0,0,1,0,0,0,dlo)
       call bvalv1 (0,0,1,0,0,0,v1 )
!
!    2) Do second portion of the interior points, plus some on borders.
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,is+2,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo ,abn )
       else
         call tranx1 (is  ,is+2,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,is+3,js+2,je  ,ks,ks,s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx1 (ie-1,ie  ,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero ,ern ,abo,abn )
       else
         call tranx1 (ie-1,ie  ,js+1,je  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (ie-2,ie  ,js+2,je  ,ks,ks,s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!
! Finally, do the remaining border zones.
! 
       if (lrad .ne. 0) then
         call tranx1 (is  ,ie  ,js  ,js  ,ks, ks, dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod & 
                     ,ero ,ern ,abo,abn )
       else
         call tranx1 (is  ,ie  ,js  ,js  ,ks, ks, dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
!
       call momx1  (is  ,ie  ,js  ,js+1,ks, ks, s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
!
! Mark d and e/d (e) boundary values out of date.
!
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0
       enddo
!
      go to 999
!=======================================================================
!      1D TRANSPORT
!=======================================================================
111   continue
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
!
! TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
!
      call bvalv1 (3,3,0,0,0,0,v1 )
      call bvalv2 (3,3,0,0,0,0,v2 )
      call bvalv3 (3,3,0,0,0,0,v3 )
!
! We need all slabs of er/d.
!
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
!
! We need all slabs of eod.
!
       call bvale  (3,3,0,0,0,0,eod)
!
      else ! xiso
!
       if (nseq .eq. 1) then
!
! We need to make a copy of the density, since we skipped pdv.
!
           do j=js,js
             do i=is-2,ie+2
               dlo(i,j,ks) = den(i,j,ks)
             enddo
           enddo
       endif
!
      endif ! xiso
!
! We need all the density slabs.
!
       call bvald  (3,3,0,0,0,0,dlo)
       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)
!
!    2) Do first portion of the interior points.
!      
       if  (lrad .ne. 0) then
         call tranx1 (is+3,ie-2,js,js  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                      ,ero,ern ,abo,abn )
       else 
         call tranx1 (is+3,ie-2,js,js  ,ks,ks,dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is+4,ie-3,js,js  ,ks,ks,s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
      
!
!    3) Wait for communications to complete.
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
!......................................................................
!
! Finally, do the remaining border zones.
!
       if (lrad .ne. 0) then
         call tranx1 (is  ,is+2,js  ,js,ks, ks, dlo,den &
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                      ,ero,ern ,abo,abn )
       else 
         call tranx1 (is  ,is+2,js  ,js,ks, ks, dlo,den & 
                     ,eod ,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (is  ,is+3,js  ,js,ks, ks, s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
!
       if (lrad .ne. 0) then
         call tranx1 (ie-1,ie ,js  ,js,ks, ks, dlo,den &
                     ,eod,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod &
                     ,ero,ern ,abo,abn )
       else
         call tranx1 (ie-1,ie ,js  ,js,ks, ks, dlo,den &
                     ,eod,edn ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
       endif
       call momx1  (ie-2,ie  ,js  ,js,ks, ks, s1,s2,s3,mflx, &
                      atwid1,vtwid,sflx,dq)
!
! Mark d and e/d (e) boundary values out of date.
!
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(7,1) = 0
       enddo
!
999    return
end subroutine advx1
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 A D V X 1                 \\\\\\\\\!
!=======================================================================
!
!
