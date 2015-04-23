!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              T R A N S P R T              \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine transprt
!
!    mln:zeus3d.transprt <------------------ controls the transport step
!                                                          october, 1987
!
!    written by: Mike Norman
!    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
!    modified 2: February, 1990 by David Clarke; incorporated into
!                ZEUS3D
!    modified 3: Feb. 15, 1996 by Robert Fiedler; completely
!                rewritten for ZEUS-MP.
!    modified 4: Dec. 20, 1996 by Robert Fiedler; added radiation.
!    modified 5: Jan. 21, 1997 by Robert Fiedler; NO_TRANSPORT switch
!    modified 6: Dec. 30, 1999 by PSLi; added update of momenta.
!
!  PURPOSE: This subroutine transports the field variables through the
!  mesh in a directionally split manner.  In each succesive call to
!  TRANSPRT, the order of the directions is permuted (resulting in
!  XYZ...YXZ...YZX...ZYX...ZXY...XZY...XYZ...etc.).  This MAY be better
!  than leaving the order the same each time (XYZ...XYZ...etc), and
!  seems to be better than unsplit schemes (Hawley).  Momenta are
!  computed from velocities in "avisc" and then transported.  Velocities
!  are not updated until the end of the transport step.  
!
!  The magnetic field components are updated by CT which is a merger (as
!  implemented by Jim Stone) of the method of characteristics and a
!  variant of the constrained transport algorithm of Evans and Hawley.
!
!  Note that the order in which variables are transported is important
!  (especially d).  
!
!  LOCAL VARIABLES:
!
!  EXTERNALS:
!    CT
!    ADVX1   , ADVX2   , ADVX3
!
!-----------------------------------------------------------------------
!
      use config
      use param
      use root
      use grid
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
      integer i,j,k
!
!-----------------------------------------------------------------------
      if(xhydro .eqv. .false.) go to 666
!
      if(xmhd) then
!
!      Transport the three components of B using Constrained Transport.
!
       call ct
      endif ! xmhd
!
!      Momentum densities were computed from velocities in the
!      artificial viscosity substep (which must therefore not be
!      skipped, even if qcon = 0.)  Momentum density boundary
!      values are not needed.
!
!PS 
      do k=ks,ke
         do j=js,je
           do i=is,ie
             w3da(i,j,k) = v1(i,j,k) * 0.5 * (d(i-1,j  ,k  ) + d(i,j,k))
             w3db(i,j,k) = v2(i,j,k) * 0.5 * (d(i  ,j-1,k  ) + d(i,j,k)) &
                           * g2b(i)
             w3dc(i,j,k) = v3(i,j,k) * 0.5 * (d(i  ,j  ,k-1) + d(i,j,k)) &
                           * g31b(i) * g32b(j)
           enddo
         enddo
       enddo
!X       if(nhy .eq. 21) then
!X        write(*,"('TRANSPRT: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
!X       endif
!
!---------------- directional split in X1-X2-X3 fashion ----------------
!
       nseq = 0  ! in /root/
       if (ix1x2x3 .eq. 1) then
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,ero,ern
!     &                  ,abo,abn
!     &                  ,mflx,s1,s2,s3)
!        
         call advx1 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
         call advx2 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
         call advx3 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!
         ix1x2x3 = 2
         goto 10
!
!---------------- directional split in X2-X1-X3 fashion ----------------
!
       else if (ix1x2x3 .eq. 2) then
!
         call advx2 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
         call advx1 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
         call advx3 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!
         ix1x2x3 = 3
         goto 10
!
!---------------- directional split in X2-X3-X1 fashion ----------------
!
       else if (ix1x2x3 .eq. 3) then
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,mflx,s1,s2,s3)
!
         call advx2 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
         call advx3 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
         call advx1 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!
         ix1x2x3 = 4
         goto 10
!
!---------------- directional split in X3-X2-X1 fashion ----------------
!
       else if (ix1x2x3 .eq. 4) then
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,mflx,s1,s2,s3)
!
         call advx3 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!X         if(nhy .eq. 21)
!X     .    write(*,"('ADVX3: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
         call advx2 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
!X         if(nhy .eq. 21)
!X     .    write(*,"('ADVX2: w3db = ',12pd16.8)")w3db(4,4,3),w3db(4,4,4)
         call advx1 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!X         if(nhy .eq. 21)
!X     .    write(*,"('ADVX1: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
!
         ix1x2x3 = 5
         goto 10
!
!---------------- directional split in X3-X1-X2 fashion ----------------
!
       else if (ix1x2x3 .eq. 5) then
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,mflx,s1,s2,s3)
!
         call advx3 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
         call advx1 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
         call advx2 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!
         ix1x2x3 = 6
         goto 10
!
!---------------- directional split in X1-X3-X2 fashion ----------------
!
       else ! if (ix1x2x3 .eq. 6) then
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,mflx,s1,s2,s3)
!
         call advx1 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
         call advx3 (d   ,w3dd &
                    ,w3dg,w3de &
                    ,w3dh,er &
                    ,w4da,abun &
                    ,w3df,w3da,w3db,w3dc)
         call advx2 (w3dd,d &
                    ,w3de,w3dg &
                    ,er  ,w3dh &
                    ,abun,w4da &
                    ,w3df,w3da,w3db,w3dc)
!
         ix1x2x3 = 1
         goto 10
       endif
!
!-----------------------------------------------------------------------
!
10     continue
!
! Mark momentum density (velocity) boundary values out of date.  
! The d and e boundary values were maked out of date in advx*.
!
       do 20 i=1,6
         bvstat(i,3) = 0  !  v1
         bvstat(i,4) = 0  !  v2
         bvstat(i,5) = 0  !  v3
20     continue
!
! The velocities need to be computed from momentum densities.
! This will be done in nudt/newdt.
!
      go to 777
666   continue
      if(lrad .ne. 0) then
!
! Skipping the transport step.  Emulate it by copying er to w3dh and 
! marking the er boundary data out of date.
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             w3dh(i,j,k) = er(i,j,k)  !  Really er/d; fixed in newdt
           enddo ! i
         enddo ! j
       enddo ! k
!
       do i=1,6
         bvstat(i,6) = 0  !  er
       enddo
!
      endif ! lrad
!
777   return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////              T R A N S P R T              \\\\\\\\\\
!
!=======================================================================
!
