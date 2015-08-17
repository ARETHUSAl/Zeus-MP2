!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              T R A N S P R T              \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine transprt_1D
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
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer i,j,k,l
!
!-----------------------------------------------------------------------
      if(xhydro .eqv. .false.) go to 666
!
      if(xmhd) then
!
!      Transport the three components of B using Constrained Transport.
!
       call ct_1D
      endif ! xmhd
!
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
             w3db(i,j,k) = v2(i,j,k) * 0.5 * (d(i  ,j  ,k  ) + d(i,j,k)) &
                           * g2b(i)
             w3dc(i,j,k) = v3(i,j,k) * 0.5 * (d(i  ,j  ,k  ) + d(i,j,k)) &
                           * g31b(i) * g32b(j)
           enddo
         enddo
       enddo
!
!---------------- directional split in X1-X2-X3 fashion ----------------
!
       nseq = 0  ! in /root/
!
!       subroutine advx1 (dlo,den
!     &                  ,eod,edn
!     &                  ,ero,ern
!     &                  ,mflx,s1,s2,s3)
!        
         if (lrad .ne. 0) then
           call advx1 (w3dd,d   ,w3de,w3dg,w3df,w3da,w3db,w3dc &
                      ,er  ,w3dh,abun,w4da)
         else
           call advx1 (w3dd,d   ,w3de,w3dg,w3df,w3da,w3db,w3dc)
         endif
!
10     continue
!
! Mark momentum density (velocity) boundary values out of date.  
! The d and e boundary values were maked out of date in advx*.
!
       do 20 i=1,6
         bvstat(i,3) = 0  !  v1
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
      endif ! lrad
!
!      Also need an array copy for ABUN and W4DA...
!
      if(nspec .gt. 1) then
       do l =  1, nspec
        do k = ks, ke
         do j = js, je
          do i = is, ie
           w4da(i,j,k,l) = abun(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      endif
!
      do i=1,6
       bvstat(i,6) = 0  !  er
       bvstat(i,7) = 0  !  abun
      enddo
!
777   return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////              T R A N S P R T              \\\\\\\\\!
!=======================================================================
!
