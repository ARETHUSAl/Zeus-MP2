!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 B V A L D                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine bvald ( rl1, ru1, rl2, ru2, rl3, ru3, d )
!
!    dac:zeus3d.bvald <------------------------- density boundary values
!    from mln:zeus04.bval; jms:zeus2d.bvald               february, 1990
!
!    written by: David Clarke, February, 1990.
!    modified 1: RAF, 3/5/96; completely rewritten for ZEUS-MP
!    modified 2: RAF, 8/27/96; correction for periodic BCs iib, no MPI.
!
!  PURPOSE: This routine sets boundary values for the density.  The
!  active zones for "d" are "is" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation. The ranges for HD boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  However, edge and corner boundary values are required the MOC 
!  algorithm.
!
!  Boundary values are set for the first two zones beyond the boundary
!  to allow for third order interpolations.
!
!  Boundary values are set according to the the basic flow types:
!
!      nflo = 1  =>  reflecting [v(normal) = b(normal) = 0]
!           =-1  =>  reflecting (XYZ: same as 1; ZRP: same as 1 with
!                    inversion of 3-components at ijb; RTP: same as 1
!                    with inversion of 2- and 3-components at iib and
!                    inversion of 3-components at ijb and ojb.)
!           = 2  =>  flow out
!           = 3  =>  flow in
!           = 4  =>  periodic
!           = 5  =>  reflecting [v(normal) = 0, b(tangential) = 0]
!
!  If desired, every boundary zone may be given a different boundary
!  type.  These types are stored in the following six arrays:
!
!      niib(j,k) = nflo of inner i boundary on sweep j,k
!      noib(j,k) = nflo of outer i boundary on sweep j,k
!      nijb(k,i) = nflo of inner j boundary on sweep k,i
!      nojb(k,i) = nflo of outer j boundary on sweep k,i
!      nikb(i,j) = nflo of inner k boundary on sweep i,j
!      nokb(i,j) = nflo of outer k boundary on sweep i,j
!
!  In addition, there are "niib2", "niib3", and "niib23" to account for
!  the different centring of the variables.  Similar arrays are defined
!  for all other boundaries (see discussion in BNDYFLGS).
!
!  Note that there is no point in setting the boundaries if the grid
!  stretching routine (EXTEND) is used, until the solution actually
!  reaches the edge of the grid.
!
!  Flags rl1, ru1, rl2, ru2, rl3, and ru3 request boundary data:
!      = 0  => do nothing
!      = 1  => pass 1st         slab
!      = 3  => pass 1st and 2nd slab
!
!  Array bvstat in /bndryi/ records the status of boundary values:
!      = 0  => needs updating
!      = 1  => 1st         slab  is  up to date, but not the second.
!      = 3  => 1st and 2nd slabs are up to date.
!
!  The values of bvstat must be reset to 0 when the corresponding
!  function is updated in some external routine.
!
!  For the mass density only, if an "rl" = 3, we pass an extra layer to 
!  the "m" tile, allowing us to compute the mass flux at the point
!  x[123]a([ijk]s-1) in TRANX[123].
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      real(rl) :: d(in,jn,kn)
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,1)
       u1 = ru1 - bvstat(2,1)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(d (is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,1100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l1.ge.2) then
             nreq = nreq + 1
             call MPI_IRECV(diib(1,1,3)        ,  jn*kn &
             , MPI_FLOAT, n1m ,1175+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(1,1) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND(d (is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,1200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,1) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is+1,j,k)
               diib(j,k,3) = d(is+2,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is-1,j,k)
               diib(j,k,3) = d(is-2,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               d(is-1,j,k) = diib (j,k,1)
               d(is-2,j,k) = diib (j,k,2)
!              do nothing for diib(j,k,3)
             endif
             if (niib(j,k) .eq. 5) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is+1,j,k)
               diib(j,k,3) = d(is+2,j,k)
             endif
           enddo
         enddo
         bvstat(1,1) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV(d (ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,1200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,1) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(d (ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,1100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l1 .ge. 2) then
             do k=ks-2,ke+3
               do j=js-2,je+3
                 doib(j,k,3) = d(ie-2,j,k)
               enddo
             enddo
             nreq = nreq + 1
             call MPI_ISEND(doib (1,1,3), jn*kn &
             , MPI_FLOAT, n1p ,1175+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(1,1) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               d(ie+1,j,k) = d(ie,j,k)
               d(ie+2,j,k) = d(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               d(ie+1,j,k) = d(ie,j,k)
               d(ie+2,j,k) = d(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               d(ie+1,j,k) = doib (j,k,1)
               d(ie+2,j,k) = doib (j,k,2)
             endif
             if (noib(j,k) .eq. 5) then
               d(ie+1,j,k) = d(ie  ,j,k)
               d(ie+2,j,k) = d(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,1) = ru1
         endif
       endif
#endif /* MPI_USED */
#ifndef MPI_USED
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is+1,j,k)
               diib(j,k,3) = d(is+2,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is-1,j,k)
               diib(j,k,3) = d(is-2,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               d(is-1,j,k) = diib (j,k,1)
               d(is-2,j,k) = diib (j,k,2)
!              do nothing for diib(j,k,3)
             endif
             if (niib(j,k) .eq. 4) then
               d(is-1,j,k) = d(ie  ,j,k)
               d(is-2,j,k) = d(ie-1,j,k)
               diib(j,k,3) = d(ie-2,j,k)
             endif
             if (niib(j,k) .eq. 5) then
               d(is-1,j,k) = d(is  ,j,k)
               d(is-2,j,k) = d(is+1,j,k)
               diib(j,k,3) = d(is+2,j,k)
             endif
           enddo
         enddo
         bvstat(1,1) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               d(ie+1,j,k) = d(ie,j,k)
               d(ie+2,j,k) = d(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               d(ie+1,j,k) = d(ie,j,k)
               d(ie+2,j,k) = d(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               d(ie+1,j,k) = doib (j,k,1)
               d(ie+2,j,k) = doib (j,k,2)
             endif
             if (noib(j,k) .eq. 4) then
               d(ie+1,j,k) = d(is  ,j,k)
               d(ie+2,j,k) = d(is+1,j,k)
             endif
             if (noib(j,k) .eq. 5) then
               d(ie+1,j,k) = d(ie  ,j,k)
               d(ie+2,j,k) = d(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,1) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,1)
       u2 = ru2 - bvstat(4,1)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(d (   1,js-ll+j-1,   1), 1, j_slice, n2m &
                           ,1300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l2 .ge. 2) then
             nreq = nreq + 1
             call MPI_IRECV(dijb(1,1,3)        , in*kn &
             , MPI_FLOAT, n2m ,1375+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(3,1) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND(d (   1,js+uu+j-1,   1), 1, j_slice, n2m &
                           ,1400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,1) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js+1,k)
               dijb(i,k,3) = d(i,js+2,k)
             endif
             if (nijb(i,k) .eq. 2) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js-1,k)
               dijb(i,k,3) = d(i,js-2,k)
             endif
             if (nijb(i,k) .eq. 3) then
               d(i,js-1,k) = dijb (i,k,1)
               d(i,js-2,k) = dijb (i,k,2)
!              do nothing for dijb(i,k,3)
             endif
             if (nijb(i,k) .eq. 5) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js+1,k)
               dijb(i,k,3) = d(i,js+2,k)
             endif
           enddo
         enddo
         bvstat(3,1) = rl2
         endif
       endif
!
!      Outer j boundary.
!
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV(d (   1,je+j+uu,   1), 1, j_slice, n2p &
                           ,1400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,1) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(d (   1,je+j-ll,   1), 1, j_slice, n2p &
                           ,1300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l2 .ge. 2) then
             do k=ks-2,ke+3
               do i=is-2,ie+3
                 dojb(i,k,3) = d(i,je-2,k)
               enddo
             enddo
             nreq = nreq + 1
             call MPI_ISEND(dojb(1,1,3)          , in*kn &
             , MPI_FLOAT, n2p,1375+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(3,1) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               d(i,je+1,k) = dojb (i,k,1)
               d(i,je+2,k) = dojb (i,k,2)
             endif
             if (nojb(i,k) .eq. 5) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,1) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js+1,k)
               dijb(i,k,3) = d(i,js+2,k)
             endif
             if (nijb(i,k) .eq. 2) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js-1,k)
               dijb(i,k,3) = d(i,js-2,k)
             endif
             if (nijb(i,k) .eq. 3) then
               d(i,js-1,k) = dijb (i,k,1)
               d(i,js-2,k) = dijb (i,k,2)
!              do nothing for dijb(i,k,3)
             endif
             if (nijb(i,k) .eq. 4) then
               d(i,js-1,k) = d(i,je  ,k)
               d(i,js-2,k) = d(i,je-1,k)
               dijb(i,k,3) = d(i,je-2,k)
             endif
             if (nijb(i,k) .eq. 5) then
               d(i,js-1,k) = d(i,js  ,k)
               d(i,js-2,k) = d(i,js+1,k)
               dijb(i,k,3) = d(i,js+2,k)
             endif
           enddo
         enddo
         bvstat(3,1) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               d(i,je+1,k) = dojb (i,k,1)
               d(i,je+2,k) = dojb (i,k,2)
             endif
             if (nojb(i,k) .eq. 4) then
               d(i,je+1,k) = d(i,js  ,k)
               d(i,je+2,k) = d(i,js+1,k)
             endif
             if (nojb(i,k) .eq. 5) then
               d(i,je+1,k) = d(i,je  ,k)
               d(i,je+2,k) = d(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,1) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,1)
       u3 = ru3 - bvstat(6,1)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(d (   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,1500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l3 .ge. 2) then
           nreq = nreq + 1
             call MPI_IRECV(dikb(1,1,3)        , in*jn &
             , MPI_FLOAT, n3m ,1575+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(5,1) = rl3
         endif
         if (u3 .gt. 0) then
         do k=1,us
             nreq = nreq + 1
             call MPI_ISEND(d (   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,1600+25*k+nsub, comm3d, req(nreq), ierr)
         enddo
         bvstat(6,1) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks+1)
               dikb(i,j,3) = d(i,j,ks+2)
             endif
             if (nikb(i,j) .eq. 2) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks-1)
               dikb(i,j,3) = d(i,j,ks-2)
             endif
             if (nikb(i,j) .eq. 3) then
               d(i,j,ks-1) = dikb (i,j,1)
               d(i,j,ks-2) = dikb (i,j,2)
!              do nothing for dikb(i,j,3)
             endif
             if (nikb(i,j) .eq. 5) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks+1)
               dikb(i,j,3) = d(i,j,ks+2)
             endif
           enddo
         enddo
         bvstat(5,1) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV(d (   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,1600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,1) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(d (   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,1500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           if (l3 .ge. 2) then
             do j=js-2,je+3
               do i=is-2,ie+3
                 dokb(i,j,3) = d(i,j,ke-2)
               enddo
             enddo
             nreq = nreq + 1
             call MPI_ISEND(dokb(1,1,3)          , in*jn &
             , MPI_FLOAT, n3p,1575+nsub, comm3d &
             , req(nreq), ierr)
           endif
           bvstat(5,1) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               d(i,j,ke+1) = dokb (i,j,1)
               d(i,j,ke+2) = dokb (i,j,2)
             endif
             if (nokb(i,j) .eq. 5) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,1) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks+1)
               dikb(i,j,3) = d(i,j,ks+2)
             endif
             if (nikb(i,j) .eq. 2) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks-1)
               dikb(i,j,3) = d(i,j,ks-2)
             endif
             if (nikb(i,j) .eq. 3) then
               d(i,j,ks-1) = dikb (i,j,1)
               d(i,j,ks-2) = dikb (i,j,2)
!              do nothing for dikb(i,j,3)
             endif
             if (nikb(i,j) .eq. 4) then
               d(i,j,ks-1) = d(i,j,ke  )
               d(i,j,ks-2) = d(i,j,ke-1)
               dikb(i,j,3) = d(i,j,ke-2)
             endif
             if (nikb(i,j) .eq. 5) then
               d(i,j,ks-1) = d(i,j,ks  )
               d(i,j,ks-2) = d(i,j,ks+1)
               dikb(i,j,3) = d(i,j,ks+2)
             endif
           enddo
         enddo
         bvstat(5,1) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               d(i,j,ke+1) = dokb (i,j,1)
               d(i,j,ke+2) = dokb (i,j,2)
             endif
             if (nokb(i,j) .eq. 4) then
               d(i,j,ke+1) = d(i,j,ks  )
               d(i,j,ke+2) = d(i,j,ks+1)
             endif
             if (nokb(i,j) .eq. 5) then
               d(i,j,ke+1) = d(i,j,ke  )
               d(i,j,ke+2) = d(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,1) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 B V A L D                 \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 B V A L E                 \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvale ( rl1, ru1, rl2, ru2, rl3, ru3, e )
!
!    dac:zeus3d.bvale <--------- internal energy density boundary values
!    from mln:zeus04.bval; jms:zeus2d.bvale               february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!    modified 2: RAF, 8/27/96; correction for periodic BCs iib, no MPI.
!
!  PURPOSE: This routine sets boundary values for the internal energy.
!  The active zones for "e" are "is" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation.  No edge or corner boundary values are required.
!  Thus, the ranges for the boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!  Flag iupper activates enables sends    in the "m" direction and
!                                receives in the "p" direction when
!  nonzero. 
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      real(rl) :: e(in,jn,kn)
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,2)
       u1 = ru1 - bvstat(2,2)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV( e(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,2100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,2) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND( e(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                         ,2200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,2) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is+1,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is-1,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               if (l1.gt.1 ) then   !  PASS e/d (NOT TOTAL_ENERGY)
                 e(is-1,j,k) = eiib (j,k,1) / diib (j,k,1)
                 e(is-2,j,k) = eiib (j,k,2) / diib (j,k,2)
               else                 !  PASS e
                 e(is-1,j,k) = eiib (j,k,1)
                 e(is-2,j,k) = eiib (j,k,2)
               endif
             endif
             if (niib(j,k) .eq. 5) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,2) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV( e(ie+i+uu,   1,   1), 1, i_slice, n1p &
                         ,2200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,2) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND( e(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,2100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,2) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               e(ie+1,j,k) = e(ie,j,k)
               e(ie+2,j,k) = e(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               e(ie+1,j,k) = e(ie,j,k)
               e(ie+2,j,k) = e(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               if (l1.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(ie+1,j,k) = eoib (j,k,1) / doib (j,k,1)
                 e(ie+2,j,k) = eoib (j,k,2) / doib (j,k,2)
               else                 !  Pass e
                 e(ie+1,j,k) = eoib (j,k,1)
                 e(ie+2,j,k) = eoib (j,k,2)
               endif
             endif
             if (noib(j,k) .eq. 5) then
               e(ie+1,j,k) = e(ie  ,j,k)
               e(ie+2,j,k) = e(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,2) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is+1,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is-1,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               if (l1.gt.1 ) then   !  PASS e/d (NOT TOTAL_ENERGY)
                 e(is-1,j,k) = eiib (j,k,1) / diib (j,k,1)
                 e(is-2,j,k) = eiib (j,k,2) / diib (j,k,2)
               else                 !  PASS e
                 e(is-1,j,k) = eiib (j,k,1)
                 e(is-2,j,k) = eiib (j,k,2)
               endif
             endif
             if (niib(j,k) .eq. 4) then
               e(is-1,j,k) = e(ie  ,j,k)
               e(is-2,j,k) = e(ie-1,j,k)
             endif
             if (niib(j,k) .eq. 5) then
               e(is-1,j,k) = e(is  ,j,k)
               e(is-2,j,k) = e(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,2) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               e(ie+1,j,k) = e(ie,j,k)
               e(ie+2,j,k) = e(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               e(ie+1,j,k) = e(ie,j,k)
               e(ie+2,j,k) = e(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               if (l1.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(ie+1,j,k) = eoib (j,k,1) / doib (j,k,1)
                 e(ie+2,j,k) = eoib (j,k,2) / doib (j,k,2)
               else                 !  Pass e
                 e(ie+1,j,k) = eoib (j,k,1)
                 e(ie+2,j,k) = eoib (j,k,2)
               endif
             endif
             if (noib(j,k) .eq. 4) then
               e(ie+1,j,k) = e(is  ,j,k)
               e(ie+2,j,k) = e(is+1,j,k)
             endif
             if (noib(j,k) .eq. 5) then
               e(ie+1,j,k) = e(ie  ,j,k)
               e(ie+2,j,k) = e(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,2) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,2)
       u2 = ru2 - bvstat(4,2)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV( e(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                         ,2300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,2) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND( e(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                         ,2400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,2) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js+1,k)
             endif
             if (nijb(i,k) .eq. 2) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 3) then
               if (l2.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,js-1,k) = eijb (i,k,1) / dijb (i,k,1)
                 e(i,js-2,k) = eijb (i,k,2) / dijb (i,k,2)
               else                 !  Pass e
                 e(i,js-1,k) = eijb (i,k,1)
                 e(i,js-2,k) = eijb (i,k,2)
               endif
             endif
             if (nijb(i,k) .eq. 5) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,2) = rl2
         endif
       endif
!
!      Outer j boundary.
!
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV( e(   1,je+j+uu,   1), 1, j_slice, n2p &
                         ,2400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,2) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND( e(   1,je+j-ll,   1), 1, j_slice, n2p &
                         ,2300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,2) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               if (l2.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,je+1,k) = eojb (i,k,1) / dojb (i,k,1)
                 e(i,je+2,k) = eojb (i,k,2) / dojb (i,k,2)
               else                 !  Pass e
                 e(i,je+1,k) = eojb (i,k,1)
                 e(i,je+2,k) = eojb (i,k,2)
               endif
             endif
             if (nojb(i,k) .eq. 5) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,2) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js+1,k)
             endif
             if (nijb(i,k) .eq. 2) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 3) then
               if (l2.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,js-1,k) = eijb (i,k,1) / dijb (i,k,1)
                 e(i,js-2,k) = eijb (i,k,2) / dijb (i,k,2)
               else                 !  Pass e
                 e(i,js-1,k) = eijb (i,k,1)
                 e(i,js-2,k) = eijb (i,k,2)
               endif
             endif
             if (nijb(i,k) .eq. 4) then
               e(i,js-1,k) = e(i,je  ,k)
               e(i,js-2,k) = e(i,je-1,k)
             endif
             if (nijb(i,k) .eq. 5) then
               e(i,js-1,k) = e(i,js  ,k)
               e(i,js-2,k) = e(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,2) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               if (l2.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,je+1,k) = eojb (i,k,1) / dojb (i,k,1)
                 e(i,je+2,k) = eojb (i,k,2) / dojb (i,k,2)
               else                 !  Pass e
                 e(i,je+1,k) = eojb (i,k,1)
                 e(i,je+2,k) = eojb (i,k,2)
               endif
             endif
             if (nojb(i,k) .eq. 4) then
               e(i,je+1,k) = e(i,js  ,k)
               e(i,je+2,k) = e(i,js+1,k)
             endif
             if (nojb(i,k) .eq. 5) then
               e(i,je+1,k) = e(i,je  ,k)
               e(i,je+2,k) = e(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,2) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,2)
       u3 = ru3 - bvstat(6,2)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV( e(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,2500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,2) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND( e(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,2600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,2) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks+1)
             endif
             if (nikb(i,j) .eq. 2) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 3) then
               if (l3.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,j,ks-1) = eikb (i,j,1) / dikb (i,j,1)
                 e(i,j,ks-2) = eikb (i,j,2) / dikb (i,j,2)
               else                 !  Pass e
                 e(i,j,ks-1) = eikb (i,j,1)
                 e(i,j,ks-2) = eikb (i,j,2)
               endif
             endif
             if (nikb(i,j) .eq. 5) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,2) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV( e(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,2600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,2) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND( e(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,2500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,2) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               if (l3.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,j,ke+1) = eokb (i,j,1) / dokb (i,j,1)
                 e(i,j,ke+2) = eokb (i,j,2) / dokb (i,j,2)
               else                 !  Pass e
                 e(i,j,ke+1) = eokb (i,j,1)
                 e(i,j,ke+2) = eokb (i,j,2)
               endif
             endif
             if (nokb(i,j) .eq. 5) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,2) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks+1)
             endif
             if (nikb(i,j) .eq. 2) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 3) then
               if (l3.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,j,ks-1) = eikb (i,j,1) / dikb (i,j,1)
                 e(i,j,ks-2) = eikb (i,j,2) / dikb (i,j,2)
               else                 !  Pass e
                 e(i,j,ks-1) = eikb (i,j,1)
                 e(i,j,ks-2) = eikb (i,j,2)
               endif
             endif
             if (nikb(i,j) .eq. 4) then
               e(i,j,ks-1) = e(i,j,ke  )
               e(i,j,ks-2) = e(i,j,ke-1)
             endif
             if (nikb(i,j) .eq. 5) then
               e(i,j,ks-1) = e(i,j,ks  )
               e(i,j,ks-2) = e(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,2) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               if (l3.gt.1 ) then   !  Pass e/d (not for TOTAL_ENERGY)
                 e(i,j,ke+1) = eokb (i,j,1) / dokb (i,j,1)
                 e(i,j,ke+2) = eokb (i,j,2) / dokb (i,j,2)
               else                 !  Pass e
                 e(i,j,ke+1) = eokb (i,j,1)
                 e(i,j,ke+2) = eokb (i,j,2)
               endif
             endif
             if (nokb(i,j) .eq. 4) then
               e(i,j,ke+1) = e(i,j,ks  )
               e(i,j,ke+2) = e(i,j,ks+1)
             endif
             if (nokb(i,j) .eq. 5) then
               e(i,j,ke+1) = e(i,j,ke  )
               e(i,j,ke+2) = e(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,2) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L E                  \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 B V A L E R S             \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalers ( rl1, ru1, rl2, ru2, rl3, ru3, er )
!
!    Written : RAF, 9/23/96 for ZEUS-MP
!    Last modified: 2/17/97
!
!  PURPOSE: This routine sets boundary values for the raditation
!  energy density.
!
!  This version is to be called from a source step routine.  If er is
!  specified on a boundary, er**b is used for boundary values.
!
!  The active zones for "t" are "is" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation.  No edge or corner boundary values are required.
!  Thus, the ranges for the boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!  Flag iupper activates enables sends    in the "m" direction and
!                                receives in the "p" direction when
!  nonzero. 
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      real(rl) :: er(in,jn,kn)
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,6)
       u1 = ru1 - bvstat(2,6)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(2).eq.0 .or. niis(2).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,6100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,6) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                         ,6200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,6) = ru1
         endif
       else
         if (l1 .gt. 0) then
#ifdef MARSHAK
         call source
#endif
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(liib(j,k)) .eq. 1) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
             if (liib(j,k) .eq. 2) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is-1,j,k)
             endif
             if (liib(j,k) .eq. 3) then
               er(is-1,j,k) = eriib (j,k,1)
               er(is-2,j,k) = eriib (j,k,2)
             endif
             if (liib(j,k) .eq. 5) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,6) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(2).eq.0 .or. nois(2).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(ie+i+uu,   1,   1), 1, i_slice, n1p &
                         ,6200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,6) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,6100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,6) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(loib(j,k)) .eq. 1) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
             if (loib(j,k) .eq. 2) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie+1,j,k)
             endif
             if (loib(j,k) .eq. 3) then
               er(ie+1,j,k) = eroib (j,k,1)
               er(ie+2,j,k) = eroib (j,k,2)
             endif
             if (loib(j,k) .eq. 5) then
               er(ie+1,j,k) = er(ie  ,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,6) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(liib(j,k)) .eq. 1) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
             if (liib(j,k) .eq. 2) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is-1,j,k)
             endif
             if (liib(j,k) .eq. 3) then
               er(is-1,j,k) = eriib (j,k,1)
               er(is-2,j,k) = eriib (j,k,2)
             endif
             if (liib(j,k) .eq. 4) then
               er(is-1,j,k) = er(ie  ,j,k)
               er(is-2,j,k) = er(ie-1,j,k)
             endif
             if (liib(j,k) .eq. 5) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,6) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(loib(j,k)) .eq. 1) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
             if (loib(j,k) .eq. 2) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie+1,j,k)
             endif
             if (loib(j,k) .eq. 3) then
               er(ie+1,j,k) = eroib (j,k,1)
               er(ie+2,j,k) = eroib (j,k,2)
             endif
             if (loib(j,k) .eq. 4) then
               er(ie+1,j,k) = er(is  ,j,k)
               er(ie+2,j,k) = er(is+1,j,k)
             endif
             if (loib(j,k) .eq. 5) then
               er(ie+1,j,k) = er(ie  ,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,6) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,6)
       u2 = ru2 - bvstat(4,6)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(2).eq.0 .or. nijs(2).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                         ,6300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,6) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                         ,6400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,6) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lijb(i,k)) .eq. 1) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
             if (lijb(i,k) .eq. 2) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js-1,k)
             endif
             if (lijb(i,k) .eq. 3) then
               er(i,js-1,k) = erijb (i,k,1)
               er(i,js-2,k) = erijb (i,k,2)
             endif
             if (lijb(i,k) .eq. 5) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,6) = rl2
         endif
       endif
!
!      Outer j boundary.
!
       if (nojs(2).eq.0 .or. nojs(2).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(   1,je+j+uu,   1), 1, j_slice, n2p &
                         ,6400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,6) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(   1,je+j-ll,   1), 1, j_slice, n2p &
                         ,6300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,6) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lojb(i,k)) .eq. 1) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
             if (lojb(i,k) .eq. 2) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je+1,k)
             endif
             if (lojb(i,k) .eq. 3) then
               er(i,je+1,k) = erojb (i,k,1)
               er(i,je+2,k) = erojb (i,k,2)
             endif
             if (lojb(i,k) .eq. 5) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,6) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lijb(i,k)) .eq. 1) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
             if (lijb(i,k) .eq. 2) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js-1,k)
             endif
             if (lijb(i,k) .eq. 3) then
               er(i,js-1,k) = erijb (i,k,1)
               er(i,js-2,k) = erijb (i,k,2)
             endif
             if (lijb(i,k) .eq. 4) then
               er(i,js-1,k) = er(i,je  ,k)
               er(i,js-2,k) = er(i,je-1,k)
             endif
             if (lijb(i,k) .eq. 5) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,6) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lojb(i,k)) .eq. 1) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
             if (lojb(i,k) .eq. 2) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je+1,k)
             endif
             if (lojb(i,k) .eq. 3) then
               er(i,je+1,k) = erojb (i,k,1)
               er(i,je+2,k) = erojb (i,k,2)
             endif
             if (lojb(i,k) .eq. 4) then
               er(i,je+1,k) = er(i,js  ,k)
               er(i,je+2,k) = er(i,js+1,k)
             endif
             if (lojb(i,k) .eq. 5) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,6) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,6)
       u3 = ru3 - bvstat(6,6)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(2).eq.0 .or. niks(2).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,6500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,6) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,6600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,6) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(likb(i,j)) .eq. 1) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
             if (likb(i,j) .eq. 2) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks-1)
             endif
             if (likb(i,j) .eq. 3) then
               er(i,j,ks-1) = erikb (i,j,1)
               er(i,j,ks-2) = erikb (i,j,2)
             endif
             if (likb(i,j) .eq. 5) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,6) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(2).eq.0 .or. noks(2).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,6600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,6) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,6500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,6) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lokb(i,j)) .eq. 1) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
             if (lokb(i,j) .eq. 2) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke+1)
             endif
             if (lokb(i,j) .eq. 3) then
               er(i,j,ke+1) = erokb (i,j,1)
               er(i,j,ke+2) = erokb (i,j,2)
             endif
             if (lokb(i,j) .eq. 5) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,6) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(likb(i,j)) .eq. 1) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
             if (likb(i,j) .eq. 2) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks-1)
             endif
             if (likb(i,j) .eq. 3) then
               er(i,j,ks-1) = erikb (i,j,1)
               er(i,j,ks-2) = erikb (i,j,2)
             endif
             if (likb(i,j) .eq. 4) then
               er(i,j,ks-1) = er(i,j,ke  )
               er(i,j,ks-2) = er(i,j,ke-1)
             endif
             if (likb(i,j) .eq. 5) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,6) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lokb(i,j)) .eq. 1) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
             if (lokb(i,j) .eq. 2) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke+1)
             endif
             if (lokb(i,j) .eq. 3) then
               er(i,j,ke+1) = erokb (i,j,1)
               er(i,j,ke+2) = erokb (i,j,2)
             endif
             if (lokb(i,j) .eq. 4) then
               er(i,j,ke+1) = er(i,j,ks  )
               er(i,j,ke+2) = er(i,j,ks+1)
             endif
             if (lokb(i,j) .eq. 5) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,6) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L E R S              \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 B V A L E R T             \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalert ( rl1, ru1, rl2, ru2, rl3, ru3, er )
!
!    Written : RAF, 9/23/96 for ZEUS-MP
!    Last modified: 2/17/97
!
!  PURPOSE: This routine sets boundary values for the raditation
!  energy density.
!
!  This version is to be called from a transport step routine.  If er is
!  specified on a boundary, er**b/d**b is used for boundary values.
!
!  The active zones for "t" are "is" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation.  No edge or corner boundary values are required.
!  Thus, the ranges for the boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!  Flag iupper activates enables sends    in the "m" direction and
!                                receives in the "p" direction when
!  nonzero. 
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      real(rl) :: er(in,jn,kn)
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,6)
       u1 = ru1 - bvstat(2,6)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(2).eq.0 .or. niis(2).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,6100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,6) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                         ,6200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,6) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(liib(j,k)) .eq. 1) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
             if (liib(j,k) .eq. 2) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is-1,j,k)
             endif
             if (liib(j,k) .eq. 3) then
               er(is-1,j,k) = eriib (j,k,1) / diib (j,k,1)
               er(is-2,j,k) = eriib (j,k,2) / diib (j,k,2)
             endif
             if (liib(j,k) .eq. 5) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,6) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(2).eq.0 .or. nois(2).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(ie+i+uu,   1,   1), 1, i_slice, n1p &
                         ,6200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,6) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,6100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,6) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(loib(j,k)) .eq. 1) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
             if (loib(j,k) .eq. 2) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie+1,j,k)
             endif
             if (loib(j,k) .eq. 3) then
               er(ie+1,j,k) = eroib (j,k,1) / doib (j,k,1)
               er(ie+2,j,k) = eroib (j,k,2) / doib (j,k,2)
             endif
             if (loib(j,k) .eq. 5) then
               er(ie+1,j,k) = er(ie  ,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,6) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(liib(j,k)) .eq. 1) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
             if (liib(j,k) .eq. 2) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is-1,j,k)
             endif
             if (liib(j,k) .eq. 3) then
               er(is-1,j,k) = eriib (j,k,1) / diib (j,k,1)
               er(is-2,j,k) = eriib (j,k,2) / diib (j,k,2)
             endif
             if (liib(j,k) .eq. 4) then
               er(is-1,j,k) = er(ie  ,j,k)
               er(is-2,j,k) = er(ie-1,j,k)
             endif
             if (liib(j,k) .eq. 5) then
               er(is-1,j,k) = er(is  ,j,k)
               er(is-2,j,k) = er(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,6) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(loib(j,k)) .eq. 1) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
             if (loib(j,k) .eq. 2) then
               er(ie+1,j,k) = er(ie,j,k)
               er(ie+2,j,k) = er(ie+1,j,k)
             endif
             if (loib(j,k) .eq. 3) then
               er(ie+1,j,k) = eroib (j,k,1) / doib (j,k,1)
               er(ie+2,j,k) = eroib (j,k,2) / doib (j,k,2)
             endif
             if (loib(j,k) .eq. 4) then
               er(ie+1,j,k) = er(is  ,j,k)
               er(ie+2,j,k) = er(is+1,j,k)
             endif
             if (loib(j,k) .eq. 5) then
               er(ie+1,j,k) = er(ie  ,j,k)
               er(ie+2,j,k) = er(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,6) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,6)
       u2 = ru2 - bvstat(4,6)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(2).eq.0 .or. nijs(2).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                         ,6300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,6) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                         ,6400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,6) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lijb(i,k)) .eq. 1) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
             if (lijb(i,k) .eq. 2) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js-1,k)
             endif
             if (lijb(i,k) .eq. 3) then
               er(i,js-1,k) = erijb (i,k,1) / dijb (i,k,1)
               er(i,js-2,k) = erijb (i,k,2) / dijb (i,k,2)
             endif
             if (lijb(i,k) .eq. 5) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,6) = rl2
         endif
       endif
!
!      Outer j boundary.
!
       if (nojs(2).eq.0 .or. nojs(2).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(   1,je+j+uu,   1), 1, j_slice, n2p &
                         ,6400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,6) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(   1,je+j-ll,   1), 1, j_slice, n2p &
                         ,6300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,6) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lojb(i,k)) .eq. 1) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
             if (lojb(i,k) .eq. 2) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je+1,k)
             endif
             if (lojb(i,k) .eq. 3) then
               er(i,je+1,k) = erojb (i,k,1) / dojb (i,k,1)
               er(i,je+2,k) = erojb (i,k,2) / dojb (i,k,2)
             endif
             if (lojb(i,k) .eq. 5) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,6) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lijb(i,k)) .eq. 1) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
             if (lijb(i,k) .eq. 2) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js-1,k)
             endif
             if (lijb(i,k) .eq. 3) then
               er(i,js-1,k) = erijb (i,k,1) / dijb (i,k,1)
               er(i,js-2,k) = erijb (i,k,2) / dijb (i,k,2)
             endif
             if (lijb(i,k) .eq. 4) then
               er(i,js-1,k) = er(i,je  ,k)
               er(i,js-2,k) = er(i,je-1,k)
             endif
             if (lijb(i,k) .eq. 5) then
               er(i,js-1,k) = er(i,js  ,k)
               er(i,js-2,k) = er(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,6) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lojb(i,k)) .eq. 1) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
             if (lojb(i,k) .eq. 2) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je+1,k)
             endif
             if (lojb(i,k) .eq. 3) then
               er(i,je+1,k) = erojb (i,k,1) / dojb (i,k,1)
               er(i,je+2,k) = erojb (i,k,2) / dojb (i,k,2)
             endif
             if (lojb(i,k) .eq. 4) then
               er(i,je+1,k) = er(i,js  ,k)
               er(i,je+2,k) = er(i,js+1,k)
             endif
             if (lojb(i,k) .eq. 5) then
               er(i,je+1,k) = er(i,je  ,k)
               er(i,je+2,k) = er(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,6) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,6)
       u3 = ru3 - bvstat(6,6)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(2).eq.0 .or. niks(2).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV( er(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,6500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,6) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND( er(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,6600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,6) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(likb(i,j)) .eq. 1) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
             if (likb(i,j) .eq. 2) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks-1)
             endif
             if (likb(i,j) .eq. 3) then
               er(i,j,ks-1) = erikb (i,j,1) / dikb (i,j,1)
               er(i,j,ks-2) = erikb (i,j,2) / dikb (i,j,2)
             endif
             if (likb(i,j) .eq. 5) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,6) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(2).eq.0 .or. noks(2).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV( er(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,6600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,6) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND( er(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,6500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,6) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lokb(i,j)) .eq. 1) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
             if (lokb(i,j) .eq. 2) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke+1)
             endif
             if (lokb(i,j) .eq. 3) then
               er(i,j,ke+1) = erokb (i,j,1) / dokb (i,j,1)
               er(i,j,ke+2) = erokb (i,j,2) / dokb (i,j,2)
             endif
             if (lokb(i,j) .eq. 5) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,6) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(likb(i,j)) .eq. 1) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
             if (likb(i,j) .eq. 2) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks-1)
             endif
             if (likb(i,j) .eq. 3) then
               er(i,j,ks-1) = erikb (i,j,1) / dikb (i,j,1)
               er(i,j,ks-2) = erikb (i,j,2) / dikb (i,j,2)
             endif
             if (likb(i,j) .eq. 4) then
               er(i,j,ks-1) = er(i,j,ke  )
               er(i,j,ks-2) = er(i,j,ke-1)
             endif
             if (likb(i,j) .eq. 5) then
               er(i,j,ks-1) = er(i,j,ks  )
               er(i,j,ks-2) = er(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,6) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(lokb(i,j)) .eq. 1) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
             if (lokb(i,j) .eq. 2) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke+1)
             endif
             if (lokb(i,j) .eq. 3) then
               er(i,j,ke+1) = erokb (i,j,1) / dokb (i,j,1)
               er(i,j,ke+2) = erokb (i,j,2) / dokb (i,j,2)
             endif
             if (lokb(i,j) .eq. 4) then
               er(i,j,ke+1) = er(i,j,ks  )
               er(i,j,ke+2) = er(i,j,ks+1)
             endif
             if (lokb(i,j) .eq. 5) then
               er(i,j,ke+1) = er(i,j,ke  )
               er(i,j,ke+2) = er(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,6) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L E R                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L V 1                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalv1 ( rl1, ru1, rl2, ru2, rl3, ru3, v1 )
!
!    dac:zeus3d.bvalv1 <----------- 1-direction velocity boundary values
!    from mln:zeus04.bflo; jms:zeus2d.bvalv1              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 1-velocity.  The
!  active zones for "v1" are "is+1" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation.  No edge or corner boundary values are required.
!  Thus, the ranges for the boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is+1, ie                      k = ks  , ke
!    k-boundaries:   i = is+1, ie     j = js  , je
!
!  Note that for periodic or tile-tile boundaries, "is" is also active.
!
!  The flow out boundary uses a switch to ensure fluid can only flow OUT
!  of the i boundary (boundary value set to 0 if it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3 activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array v1 is input so that velocity values at old time levels
!  as well as momentum components can be passed.
!
!  NOTE: Need to know whether to pass velocity or momentum density
!        for inflow boundaries.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
      real(rl) :: v1(in,jn,kn)
      real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,3)
       u1 = ru1 - bvstat(2,3)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v1(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                         ,3100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,3) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND(v1(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,3200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,3) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               v1(is  ,j,k) =       vg1(is)
               v1(is-1,j,k) = 2.0 * vg1(is) - v1(is+1,j,k)
               v1(is-2,j,k) = 2.0 * vg1(is) - v1(is+2,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               q1           = sign ( haf, v1(is+1,j,k) - vg1(is) )
               v1(is  ,j,k) = v1(is+1,j,k) * ( 0.5 - q1 )
               v1(is-1,j,k) = v1(is  ,j,k)
               v1(is-2,j,k) = v1(is  ,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               v1(is  ,j,k) = v1iib (j,k,1)
               v1(is-1,j,k) = v1iib (j,k,2)
               v1(is-2,j,k) = 2.0 * v1iib (j,k,2) - v1iib (j,k,1)
             endif
             if (niib(j,k) .eq. 5) then
               v1(is  ,j,k) =       vg1(is)
               v1(is-1,j,k) = 2.0 * vg1(is) - v1(is+1,j,k)
               v1(is-2,j,k) = 2.0 * vg1(is) - v1(is+2,j,k)
             endif
           enddo
         enddo
         bvstat(1,3) = rl1
         endif
       endif
!
!      Outer i boundary.
!
! 1-face-centered quantities need to be evolved on only one end
! (we have chosen to evolve them on the inner boundary), hence
! no change in subscripts compared to zone-centered quantities.
!
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV(v1(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,3200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,3) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v1(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,3100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,3) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               v1(ie+1,j,k) =       vg1(ie+1)
               v1(ie+2,j,k) = 2.0 * vg1(ie+1) - v1(ie  ,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               q1           = sign ( haf, v1(ie,j,k) - vg1(ie+1) )
               v1(ie+1,j,k) = v1(ie,j,k) * ( 0.5 + q1 )
               v1(ie+2,j,k) = v1(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               v1(ie+1,j,k) = v1oib (j,k,1)
               v1(ie+2,j,k) = v1oib (j,k,2)
             endif
             if (noib(j,k) .eq. 5) then
               v1(ie+1,j,k) =       vg1(ie+1)
               v1(ie+2,j,k) = 2.0 * vg1(ie+1) - v1(ie,j,k)
             endif
           enddo
         enddo
         bvstat(2,3) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               v1(is  ,j,k) =       vg1(is)
               v1(is-1,j,k) = 2.0 * vg1(is) - v1(is+1,j,k)
               v1(is-2,j,k) = 2.0 * vg1(is) - v1(is+2,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               q1           = sign ( haf, v1(is+1,j,k) - vg1(is) )
               v1(is  ,j,k) = v1(is+1,j,k) * ( 0.5 - q1 )
               v1(is-1,j,k) = v1(is  ,j,k)
               v1(is-2,j,k) = v1(is  ,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               v1(is  ,j,k) = v1iib (j,k,1)
               v1(is-1,j,k) = v1iib (j,k,2)
               v1(is-2,j,k) = 2.0 * v1iib (j,k,2) - v1iib (j,k,1)
             endif
             if (niib(j,k) .eq. 4) then
               v1(is-1,j,k) = v1(ie  ,j,k)
               v1(is-2,j,k) = v1(ie-1,j,k)
             endif
             if (niib(j,k) .eq. 5) then
               v1(is  ,j,k) =       vg1(is)
               v1(is-1,j,k) = 2.0 * vg1(is) - v1(is+1,j,k)
               v1(is-2,j,k) = 2.0 * vg1(is) - v1(is+2,j,k)
             endif
           enddo
         enddo
         bvstat(1,3) = rl1
         endif
!
!      Outer i boundary.
!
! 1-face-centered quantities need to be evolved on only one end
! (we have chosen to evolve them on the inner boundary), hence
! no change in subscripts compared to zone-centered quantities.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               v1(ie+1,j,k) =       vg1(ie+1)
               v1(ie+2,j,k) = 2.0 * vg1(ie+1) - v1(ie  ,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               q1           = sign ( haf, v1(ie,j,k) - vg1(ie+1) )
               v1(ie+1,j,k) = v1(ie,j,k) * ( 0.5 + q1 )
               v1(ie+2,j,k) = v1(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               v1(ie+1,j,k) = v1oib (j,k,1)
               v1(ie+2,j,k) = v1oib (j,k,2)
             endif
             if (noib(j,k) .eq. 4) then
               v1(ie+1,j,k) = v1(is  ,j,k)
               v1(ie+2,j,k) = v1(is+1,j,k)
             endif
             if (noib(j,k) .eq. 5) then
               v1(ie+1,j,k) =       vg1(ie+1)
               v1(ie+2,j,k) = 2.0 * vg1(ie+1) - v1(ie,j,k)
             endif
           enddo
         enddo
         bvstat(2,3) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,3)
       u2 = ru2 - bvstat(4,3)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v1(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                           ,3300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,3) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND(v1(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                           ,3400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,3) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb1(i,k)) .eq. 1) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js+1,k)
             endif
             if (nijb1(i,k) .eq. 2) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js-1,k)
             endif
             if (nijb1(i,k) .eq. 3) then
               v1(i,js-1,k) = v1ijb (i,k,1)
               v1(i,js-2,k) = v1ijb (i,k,2)
             endif
             if (nijb1(i,k) .eq. 5) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,3) = rl2
         endif
       endif
!
!      Outer j boundary.
!
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV(v1(   1,je+j+uu,   1), 1, j_slice, n2p &
                           ,3400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,3) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v1(   1,je+j-ll,   1), 1, j_slice, n2p &
                           ,3300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,3) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb1(i,k)) .eq. 1) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je-1,k)
             endif
             if (nojb1(i,k) .eq. 2) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je+1,k)
             endif
             if (nojb1(i,k) .eq. 3) then
               v1(i,je+1,k) = v1ojb (i,k,1)
               v1(i,je+2,k) = v1ojb (i,k,2)
             endif
             if (nojb1(i,k) .eq. 5) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,3) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb1(i,k)) .eq. 1) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js+1,k)
             endif
             if (nijb1(i,k) .eq. 2) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js-1,k)
             endif
             if (nijb1(i,k) .eq. 3) then
               v1(i,js-1,k) = v1ijb (i,k,1)
               v1(i,js-2,k) = v1ijb (i,k,2)
             endif
             if (nijb1(i,k) .eq. 4) then
               v1(i,js-1,k) = v1(i,je  ,k)
               v1(i,js-2,k) = v1(i,je-1,k)
             endif
             if (nijb1(i,k) .eq. 5) then
               v1(i,js-1,k) = v1(i,js  ,k)
               v1(i,js-2,k) = v1(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,3) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb1(i,k)) .eq. 1) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je-1,k)
             endif
             if (nojb1(i,k) .eq. 2) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je+1,k)
             endif
             if (nojb1(i,k) .eq. 3) then
               v1(i,je+1,k) = v1ojb (i,k,1)
               v1(i,je+2,k) = v1ojb (i,k,2)
             endif
             if (nojb1(i,k) .eq. 4) then
               v1(i,je+1,k) = v1(i,js  ,k)
               v1(i,je+2,k) = v1(i,js+1,k)
             endif
             if (nojb1(i,k) .eq. 5) then
               v1(i,je+1,k) = v1(i,je  ,k)
               v1(i,je+2,k) = v1(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,3) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,3)
       u3 = ru3 - bvstat(6,3)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v1(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,3500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,3) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND(v1(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,3600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,3) = ru3
         endif
       else
       if (l3 .gt. 0) then
       do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb1(i,j)) .eq. 1) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks+1)
             endif
             if (nikb1(i,j) .eq. 2) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks-1)
             endif
             if (nikb1(i,j) .eq. 3) then
               v1(i,j,ks-1) = v1ikb (i,j,1)
               v1(i,j,ks-2) = v1ikb (i,j,2)
             endif
             if (nikb1(i,j) .eq. 5) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,3) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV(v1(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,3600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,3) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v1(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,3500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,3) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb1(i,j)) .eq. 1) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke-1)
             endif
             if (nokb1(i,j) .eq. 2) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke+1)
             endif
             if (nokb1(i,j) .eq. 3) then
               v1(i,j,ke+1) = v1okb (i,j,1)
               v1(i,j,ke+2) = v1okb (i,j,2)
             endif
             if (nokb1(i,j) .eq. 5) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,3) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
       if (l3 .gt. 0) then
       do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb1(i,j)) .eq. 1) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks+1)
             endif
             if (nikb1(i,j) .eq. 2) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks-1)
             endif
             if (nikb1(i,j) .eq. 3) then
               v1(i,j,ks-1) = v1ikb (i,j,1)
               v1(i,j,ks-2) = v1ikb (i,j,2)
             endif
             if (nikb1(i,j) .eq. 4) then
               v1(i,j,ks-1) = v1(i,j,ke  )
               v1(i,j,ks-2) = v1(i,j,ke-1)
             endif
             if (nikb1(i,j) .eq. 5) then
               v1(i,j,ks-1) = v1(i,j,ks  )
               v1(i,j,ks-2) = v1(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,3) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb1(i,j)) .eq. 1) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke-1)
             endif
             if (nokb1(i,j) .eq. 2) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke+1)
             endif
             if (nokb1(i,j) .eq. 3) then
               v1(i,j,ke+1) = v1okb (i,j,1)
               v1(i,j,ke+2) = v1okb (i,j,2)
             endif
             if (nokb1(i,j) .eq. 4) then
               v1(i,j,ke+1) = v1(i,j,ks  )
               v1(i,j,ke+2) = v1(i,j,ks+1)
             endif
             if (nokb1(i,j) .eq. 5) then
               v1(i,j,ke+1) = v1(i,j,ke  )
               v1(i,j,ke+2) = v1(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,3) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L V 1                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L V 2                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalv2 ( rl1, ru1, rl2, ru2, rl3, ru3, v2 )
!
!    dac:zeus3d.bvalv2 <----------- 2-direction velocity boundary values
!    from mln:zeus04.bflo; jms:zeus2d.bvalv2              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 2-velocity.  The
!  active zones for "v2" are "is" to "ie" in the 1-direction, "js+1" to
!  "je" in the 2-direction, and "ks" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for third order
!  interpolation.  No edge or corner boundary values are required.
!  Thus, the ranges for the boundary value calculations are:
!
!    i-boundaries:                    j = js+1, je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js+1, je
!
!  Note that for periodic or tile-tile boundaries, "js" is also active.
!
!  The flow out boundary uses a switch to ensure fluid can only flow OUT
!  of the j boundary (boundary value set to 0 if it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array v2 is input so that velocity values at old time levels
!  and momenta can be passed.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
      real(rl) :: v2(in,jn,kn)
      real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,4)
       u1 = ru1 - bvstat(2,4)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v2(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,4100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,4) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND(v2(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,4200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,4) = ru1
         endif
       else
!
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if (niib2(j,k) .eq. 1) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
             endif
             if (niib2(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v2(is-1,j,k) =-v2(is  ,j,k)
               v2(is-2,j,k) =-v2(is+1,j,k)
              else
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
              endif
             endif
             if (niib2(j,k) .eq. 2) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is-1,j,k)
             endif
             if (niib2(j,k) .eq. 3) then
               v2(is-1,j,k) = v2iib (j,k,1)
               v2(is-2,j,k) = v2iib (j,k,2)
             endif
             if (niib2(j,k) .eq. 5) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,4) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV(v2(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,4200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,4) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v2(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,4100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,4) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib2(j,k)) .eq. 1) then
               v2(ie+1,j,k) = v2(ie,j,k)
               v2(ie+2,j,k) = v2(ie-1,j,k)
             endif
             if (noib2(j,k) .eq. 2) then
!#if PROBLEM == advect  --   I wish cpp could do this!
!               v2(ie+1,j,k) = 2.0 * v2(ie  ,j,k) - v2(ie-1,j,k)
!               v2(ie+2,j,k) = 2.0 * v2(ie+1,j,k) - v2(ie  ,j,k)
!#else
               v2(ie+1,j,k) = v2(ie,j,k)
               v2(ie+2,j,k) = v2(ie+1,j,k)
!#endif
             endif
             if (noib2(j,k) .eq. 3) then
               v2(ie+1,j,k) = v2oib (j,k,1)
               v2(ie+2,j,k) = v2oib (j,k,2)
             endif
             if (noib2(j,k) .eq. 5) then
               v2(ie+1,j,k) = v2(ie  ,j,k)
               v2(ie+2,j,k) = v2(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,4) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if (niib2(j,k) .eq. 1) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
             endif
             if (niib2(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v2(is-1,j,k) =-v2(is  ,j,k)
               v2(is-2,j,k) =-v2(is+1,j,k)
              else
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
              endif
             endif
             if (niib2(j,k) .eq. 2) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is-1,j,k)
             endif
             if (niib2(j,k) .eq. 3) then
               v2(is-1,j,k) = v2iib (j,k,1)
               v2(is-2,j,k) = v2iib (j,k,2)
             endif
             if (niib2(j,k) .eq. 4) then
               v2(is-1,j,k) = v2(ie  ,j,k)
               v2(is-2,j,k) = v2(ie-1,j,k)
             endif
             if (niib2(j,k) .eq. 5) then
               v2(is-1,j,k) = v2(is  ,j,k)
               v2(is-2,j,k) = v2(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,4) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib2(j,k)) .eq. 1) then
               v2(ie+1,j,k) = v2(ie,j,k)
               v2(ie+2,j,k) = v2(ie-1,j,k)
             endif
             if (noib2(j,k) .eq. 2) then
!#if PROBLEM == advect  --   I wish cpp could do this!
!               v2(ie+1,j,k) = 2.0 * v2(ie  ,j,k) - v2(ie-1,j,k)
!               v2(ie+2,j,k) = 2.0 * v2(ie+1,j,k) - v2(ie  ,j,k)
!#else
               v2(ie+1,j,k) = v2(ie,j,k)
               v2(ie+2,j,k) = v2(ie+1,j,k)
!#endif
             endif
             if (noib2(j,k) .eq. 3) then
               v2(ie+1,j,k) = v2oib (j,k,1)
               v2(ie+2,j,k) = v2oib (j,k,2)
             endif
             if (noib2(j,k) .eq. 4) then
               v2(ie+1,j,k) = v2(is  ,j,k)
               v2(ie+2,j,k) = v2(is+1,j,k)
             endif
             if (noib2(j,k) .eq. 5) then
               v2(ie+1,j,k) = v2(ie  ,j,k)
               v2(ie+2,j,k) = v2(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,4) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,4)
       u2 = ru2 - bvstat(4,4)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v2(   1,js-ll+j-1,   1), 1 , j_slice, n2m &
                           ,4300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,4) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND(v2(   1,js+uu+j-1,   1), 1 , j_slice, n2m &
                           ,4400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,4) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               v2(i,js  ,k) =       vg2(js)
               v2(i,js-1,k) = 2.0 * vg2(js) - v2(i,js+1,k)
               v2(i,js-2,k) = 2.0 * vg2(js) - v2(i,js+2,k)
             endif
             if (nijb(i,k) .eq. 2) then
               q1           = sign ( haf, v2(i,js+1,k) - vg2(js) )
               v2(i,js  ,k) = v2(i,js+1,k) * ( 0.5 - q1 )
               v2(i,js-1,k) = v2(i,js  ,k)
               v2(i,js-2,k) = v2(i,js  ,k)
             endif
             if (nijb(i,k) .eq. 3) then
               v2(i,js  ,k) = v2ijb (i,k,1)
               v2(i,js-1,k) = v2ijb (i,k,2)
               v2(i,js-2,k) = 2.0 * v2ijb (i,k,2) - v2ijb (i,k,1)
             endif
             if (nijb(i,k) .eq. 5) then
               v2(i,js  ,k) =       vg2(js)
               v2(i,js-1,k) = 2.0 * vg2(js) - v2(i,js+1,k)
               v2(i,js-2,k) = 2.0 * vg2(js) - v2(i,js+2,k)
             endif
           enddo
         enddo
         bvstat(3,4) = rl2
         endif
       endif
!
!      Outer j boundary.
!
      if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV(v2(   1,je+j+uu,   1), 1 , j_slice, n2p &
                           ,4400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,4) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v2(   1,je+j-ll,   1), 1 , j_slice, n2p &
                           ,4300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,4) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               v2(i,je+1,k) =       vg2(je+1)
               v2(i,je+2,k) = 2.0 * vg2(je+1) - v2(i,je  ,k)
             endif
             if (nojb(i,k) .eq. 2) then
               q1           = sign ( haf, v2(i,je,k) - vg2(je+1) )
               v2(i,je+1,k) = v2(i,je  ,k) * ( 0.5 + q1 )
               v2(i,je+2,k) = v2(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               v2(i,je+1,k) = v2ojb (i,k,1)
               v2(i,je+2,k) = v2ojb (i,k,2)
             endif
             if (nojb(i,k) .eq. 5) then
               v2(i,je+1,k) =       vg2(je+1)
               v2(i,je+2,k) = 2.0 * vg2(je+1) - v2(i,je  ,k)
             endif
           enddo
         enddo
         bvstat(4,4) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               v2(i,js  ,k) =       vg2(js)
               v2(i,js-1,k) = 2.0 * vg2(js) - v2(i,js+1,k)
               v2(i,js-2,k) = 2.0 * vg2(js) - v2(i,js+2,k)
             endif
             if (nijb(i,k) .eq. 2) then
               q1           = sign ( haf, v2(i,js+1,k) - vg2(js) )
               v2(i,js  ,k) = v2(i,js+1,k) * ( 0.5 - q1 )
               v2(i,js-1,k) = v2(i,js  ,k)
               v2(i,js-2,k) = v2(i,js  ,k)
             endif
             if (nijb(i,k) .eq. 3) then
               v2(i,js  ,k) = v2ijb (i,k,1)
               v2(i,js-1,k) = v2ijb (i,k,2)
               v2(i,js-2,k) = 2.0 * v2ijb (i,k,2) - v2ijb (i,k,1)
             endif
             if (nijb(i,k) .eq. 4) then
               v2(i,js-1,k) = v2(i,je  ,k)
               v2(i,js-2,k) = v2(i,je-1,k)
             endif
             if (nijb(i,k) .eq. 5) then
               v2(i,js  ,k) =       vg2(js)
               v2(i,js-1,k) = 2.0 * vg2(js) - v2(i,js+1,k)
               v2(i,js-2,k) = 2.0 * vg2(js) - v2(i,js+2,k)
             endif
           enddo
         enddo
         bvstat(3,4) = rl2
         endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               v2(i,je+1,k) =       vg2(je+1)
               v2(i,je+2,k) = 2.0 * vg2(je+1) - v2(i,je  ,k)
             endif
             if (nojb(i,k) .eq. 2) then
               q1           = sign ( haf, v2(i,je,k) - vg2(je+1) )
               v2(i,je+1,k) = v2(i,je  ,k) * ( 0.5 + q1 )
               v2(i,je+2,k) = v2(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               v2(i,je+1,k) = v2ojb (i,k,1)
               v2(i,je+2,k) = v2ojb (i,k,2)
             endif
             if (nojb(i,k) .eq. 4) then
               v2(i,je+1,k) = v2(i,js  ,k)
               v2(i,je+2,k) = v2(i,js+1,k)
             endif
             if (nojb(i,k) .eq. 5) then
               v2(i,je+1,k) =       vg2(je+1)
               v2(i,je+2,k) = 2.0 * vg2(je+1) - v2(i,je  ,k)
             endif
           enddo
         enddo
         bvstat(4,4) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,4)
       u3 = ru3 - bvstat(6,4)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v2(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,4500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,4) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND(v2(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,4600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,4) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb2(i,j)) .eq. 1) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks+1)
             endif
             if (nikb2(i,j) .eq. 2) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks-1)
             endif
             if (nikb2(i,j) .eq. 3) then
               v2(i,j,ks-1) = v2ikb (i,j,1)
               v2(i,j,ks-2) = v2ikb (i,j,2)
             endif
             if (nikb2(i,j) .eq. 5) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,4) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV(v2(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,4600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,4) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v2(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,4500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,4) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb2(i,j)) .eq. 1) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke-1)
             endif
             if (nokb2(i,j) .eq. 2) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke+1)
             endif
             if (nokb2(i,j) .eq. 3) then
               v2(i,j,ke+1) = v2okb (i,j,1)
               v2(i,j,ke+2) = v2okb (i,j,2)
             endif
             if (nokb2(i,j) .eq. 5) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,4) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb2(i,j)) .eq. 1) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks+1)
             endif
             if (nikb2(i,j) .eq. 2) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks-1)
             endif
             if (nikb2(i,j) .eq. 3) then
               v2(i,j,ks-1) = v2ikb (i,j,1)
               v2(i,j,ks-2) = v2ikb (i,j,2)
             endif
             if (nikb2(i,j) .eq. 4) then
               v2(i,j,ks-1) = v2(i,j,ke  )
               v2(i,j,ks-2) = v2(i,j,ke-1)
             endif
             if (nikb2(i,j) .eq. 5) then
               v2(i,j,ks-1) = v2(i,j,ks  )
               v2(i,j,ks-2) = v2(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,4) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb2(i,j)) .eq. 1) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke-1)
             endif
             if (nokb2(i,j) .eq. 2) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke+1)
             endif
             if (nokb2(i,j) .eq. 3) then
               v2(i,j,ke+1) = v2okb (i,j,1)
               v2(i,j,ke+2) = v2okb (i,j,2)
             endif
             if (nokb2(i,j) .eq. 4) then
               v2(i,j,ke+1) = v2(i,j,ks  )
               v2(i,j,ke+2) = v2(i,j,ks+1)
             endif
             if (nokb2(i,j) .eq. 5) then
               v2(i,j,ke+1) = v2(i,j,ke  )
               v2(i,j,ke+2) = v2(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,4) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L V 2                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L V 3                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalv3 ( rl1, ru1, rl2, ru2, rl3, ru3, v3 )
!
!    dac:zeus3d.bvalv3 <----------- 3-direction velocity boundary values
!    from mln:zeus04.bflo; jms:zeus2d.bvalv3              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 3-velocity.  The
!  active zones for "v3" are "is" to "ie" in the 1-direction, "js" to
!  "je" in the 2-direction, and "ks+1" to "ke" in the 3-direction.  Two
!  layers of boundary values at each face are needed for interpolation.
!  No edge or corner boundary values are required.  Thus, the ranges for
!  the boundary value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks+1, ke
!    j-boundaries:   i = is  , ie                      k = ks+1, ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  Note that for periodic or tile-tile boundaries, "ks" is also active.
!
!  The flow out boundary uses a switch to ensure fluid can only flow OUT
!  of the k boundary (boundary value set to 0 if it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array v3 is input so that velocity values at old time levels
!  can be passed.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3
      integer  :: rl1,rl2,rl3,ru1,ru2,ru3
      integer  :: ls,ll,lu,us,ul,uu
!
      real(rl) :: v3(in,jn,kn)
      real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,5)
       u1 = ru1 - bvstat(2,5)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v3(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,5100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,5) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND(v3(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,5200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,5) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if (niib3(j,k) .eq. 1) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
             endif
             if (niib3(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3(is-1,j,k) =-v3(is  ,j,k)
               v3(is-2,j,k) =-v3(is+1,j,k)
              else
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
              endif
             endif
             if (niib3(j,k) .eq. 2) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is-1,j,k)
             endif
             if (niib3(j,k) .eq. 3) then
               v3(is-1,j,k) = v3iib (j,k,1)
               v3(is-2,j,k) = v3iib (j,k,2)
             endif
             if (niib3(j,k) .eq. 5) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,5) = rl1
         endif
       endif
!
!      Outer i boundary.
!
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV(v3(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,5200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,5) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v3(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,5100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,5) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib3(j,k)) .eq. 1) then
               v3(ie+1,j,k) = v3(ie,j,k)
               v3(ie+2,j,k) = v3(ie-1,j,k)
             endif
             if (noib3(j,k) .eq. 2) then
!#if PROBLEM == advect  --  I wish!
!               v3(ie+1,j,k) = 2.0 * v3(ie  ,j,k) - v3(ie-1,j,k)
!               v3(ie+2,j,k) = 2.0 * v3(ie+1,j,k) - v3(ie  ,j,k)
!#else
               v3(ie+1,j,k) = v3(ie,j,k)
               v3(ie+2,j,k) = v3(ie+1,j,k)
!#endif
             endif
             if (noib3(j,k) .eq. 3) then
               v3(ie+1,j,k) = v3oib (j,k,1)
               v3(ie+2,j,k) = v3oib (j,k,2)
             endif
             if (noib3(j,k) .eq. 5) then
               v3(ie+1,j,k) = v3(ie  ,j,k)
               v3(ie+2,j,k) = v3(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,5) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if (niib3(j,k) .eq. 1) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
             endif
             if (niib3(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3(is-1,j,k) =-v3(is  ,j,k)
               v3(is-2,j,k) =-v3(is+1,j,k)
              else
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
              endif
             endif
             if (niib3(j,k) .eq. 2) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is-1,j,k)
             endif
             if (niib3(j,k) .eq. 3) then
               v3(is-1,j,k) = v3iib (j,k,1)
               v3(is-2,j,k) = v3iib (j,k,2)
             endif
             if (niib3(j,k) .eq. 4) then
               v3(is-1,j,k) = v3(ie  ,j,k)
               v3(is-2,j,k) = v3(ie-1,j,k)
             endif
             if (niib3(j,k) .eq. 5) then
               v3(is-1,j,k) = v3(is  ,j,k)
               v3(is-2,j,k) = v3(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,5) = rl1
         endif
!
!      Outer i boundary.
!
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do j=js-1,je+1
             if ( abs(noib3(j,k)) .eq. 1) then
               v3(ie+1,j,k) = v3(ie,j,k)
               v3(ie+2,j,k) = v3(ie-1,j,k)
             endif
             if (noib3(j,k) .eq. 2) then
!#if PROBLEM == advect  --  I wish!
!               v3(ie+1,j,k) = 2.0 * v3(ie  ,j,k) - v3(ie-1,j,k)
!               v3(ie+2,j,k) = 2.0 * v3(ie+1,j,k) - v3(ie  ,j,k)
!#else
               v3(ie+1,j,k) = v3(ie,j,k)
               v3(ie+2,j,k) = v3(ie+1,j,k)
!#endif
             endif
             if (noib3(j,k) .eq. 3) then
               v3(ie+1,j,k) = v3oib (j,k,1)
               v3(ie+2,j,k) = v3oib (j,k,2)
             endif
             if (noib3(j,k) .eq. 4) then
               v3(ie+1,j,k) = v3(is  ,j,k)
               v3(ie+2,j,k) = v3(is+1,j,k)
             endif
             if (noib3(j,k) .eq. 5) then
               v3(ie+1,j,k) = v3(ie  ,j,k)
               v3(ie+2,j,k) = v3(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,5) = ru1
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,5)
       u2 = ru2 - bvstat(4,5)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v3(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                           ,5300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,5) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND(v3(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                           ,5400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,5) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if (nijb3(i,k) .eq. 1) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
             endif
             if (nijb3(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
              else ! lgeom
               v3(i,js-1,k) =-v3(i,js  ,k)
               v3(i,js-2,k) =-v3(i,js+1,k)
              endif ! lgeom
             endif
             if (nijb3(i,k) .eq. 2) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js-1,k)
             endif
             if (nijb3(i,k) .eq. 3) then
               v3(i,js-1,k) = v3ijb (i,k,1)
               v3(i,js-2,k) = v3ijb (i,k,2)
             endif
             if (nijb3(i,k) .eq. 5) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
             endif
         enddo
       enddo
       bvstat(3,5) = rl2
       endif
       endif
!
!      Outer j boundary.
!
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV(v3(   1,je+j+uu,   1), 1, j_slice, n2p &
                           ,5400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,5) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v3(   1,je+j-ll,   1), 1, j_slice, n2p &
                           ,5300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,5) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if (nojb3(i,k) .eq. 1) then
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
             endif
             if (nojb3(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3(i,je+1,k) =-v3(i,je  ,k)
               v3(i,je+2,k) =-v3(i,je-1,k)
              else
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
              endif
             endif
             if (nojb3(i,k) .eq. 2) then
!#if alias PROBLEM.eq.advect  -- I wish!
               v3(i,je+1,k) = 2.0 * v3(i,je  ,k) - v3(i,je-1,k)
               v3(i,je+2,k) = 2.0 * v3(i,je+1,k) - v3(i,je  ,k)
!#else
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je+1,k)
!#endif
             endif
             if (nojb3(i,k) .eq. 3) then
               v3(i,je+1,k) = v3ojb (i,k,1)
               v3(i,je+2,k) = v3ojb (i,k,2)
             endif
             if (nojb3(i,k) .eq. 5) then
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,5) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if (nijb3(i,k) .eq. 1) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
             endif
             if (nijb3(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
              else ! lgeom
               v3(i,js-1,k) =-v3(i,js  ,k)
               v3(i,js-2,k) =-v3(i,js+1,k)
              endif ! lgeom
             endif
             if (nijb3(i,k) .eq. 2) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js-1,k)
             endif
             if (nijb3(i,k) .eq. 3) then
               v3(i,js-1,k) = v3ijb (i,k,1)
               v3(i,js-2,k) = v3ijb (i,k,2)
             endif
             if (nijb3(i,k) .eq. 4) then
               v3(i,js-1,k) = v3(i,je  ,k)
               v3(i,js-2,k) = v3(i,je-1,k)
             endif
             if (nijb3(i,k) .eq. 5) then
               v3(i,js-1,k) = v3(i,js  ,k)
               v3(i,js-2,k) = v3(i,js+1,k)
             endif
         enddo
       enddo
       bvstat(3,5) = rl2
       endif
!
!      Outer j boundary.
!
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
!dir$ ivdep
           do i=is-1,ie+1
             if (nojb3(i,k) .eq. 1) then
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
             endif
             if (nojb3(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3(i,je+1,k) =-v3(i,je  ,k)
               v3(i,je+2,k) =-v3(i,je-1,k)
              else
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
              endif
             endif
             if (nojb3(i,k) .eq. 2) then
!#if alias PROBLEM.eq.advect  -- I wish!
               v3(i,je+1,k) = 2.0 * v3(i,je  ,k) - v3(i,je-1,k)
               v3(i,je+2,k) = 2.0 * v3(i,je+1,k) - v3(i,je  ,k)
!#else
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je+1,k)
!#endif
             endif
             if (nojb3(i,k) .eq. 3) then
               v3(i,je+1,k) = v3ojb (i,k,1)
               v3(i,je+2,k) = v3ojb (i,k,2)
             endif
             if (nojb3(i,k) .eq. 4) then
               v3(i,je+1,k) = v3(i,js  ,k)
               v3(i,je+2,k) = v3(i,js+1,k)
             endif
             if (nojb3(i,k) .eq. 5) then
               v3(i,je+1,k) = v3(i,je  ,k)
               v3(i,je+2,k) = v3(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,5) = ru2
         endif
#endif /* NO MPI */
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,5)
       u3 = ru3 - bvstat(6,5)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(v3(   1,   1,ks-ll+k-1), 1 , k_slice, n3m &
                           ,5500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,5) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND(v3(   1,   1,ks+uu+k-1), 1 , k_slice, n3m &
                           ,5600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,5) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               v3(i,j,ks  ) =       vg3(ks)
               v3(i,j,ks-1) = 2.0 * vg3(ks) - v3(i,j,ks+1)
               v3(i,j,ks-2) = 2.0 * vg3(ks) - v3(i,j,ks+2)
             endif
             if (nikb(i,j) .eq. 2) then
               q1           = sign ( haf, v3(i,j,ks+1) - vg3(ks) )
               v3(i,j,ks  ) = v3(i,j,ks+1) * ( 0.5 - q1 )
               v3(i,j,ks-1) = v3(i,j,ks  )
               v3(i,j,ks-2) = v3(i,j,ks  )
             endif
             if (nikb(i,j) .eq. 3) then
               v3(i,j,ks  ) = v3ikb (i,j,1)
               v3(i,j,ks-1) = v3ikb (i,j,2)
               v3(i,j,ks-2) = 2.0 * v3ikb (i,j,2) - v3ikb (i,j,1)
             endif
             if (nikb(i,j) .eq. 5) then
               v3(i,j,ks  ) =       vg3(ks)
               v3(i,j,ks-1) = 2.0 * vg3(ks) - v3(i,j,ks+1)
               v3(i,j,ks-2) = 2.0 * vg3(ks) - v3(i,j,ks+2)
             endif
           enddo
         enddo
         bvstat(5,5) = rl3
         endif
       endif
!
!      Outer k boundary.
!
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV(v3(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,5600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,5) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(v3(   1,   1,ke+k-ll), 1 , k_slice, n3p &
                           ,5500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,5) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               v3(i,j,ke+1) =       vg3(ke+1)
               v3(i,j,ke+2) = 2.0 * vg3(ke+1) - v3(i,j,ke  )
             endif
             if (nokb(i,j) .eq. 2) then
               q1           = sign ( haf, v3(i,j,ke) - vg3(ke+1) )
               v3(i,j,ke+1) = v3(i,j,ke  ) * ( 0.5 + q1 )
               v3(i,j,ke+2) = v3(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               v3(i,j,ke+1) = v3okb (i,j,1)
               v3(i,j,ke+2) = v3okb (i,j,2)
             endif
             if (nokb(i,j) .eq. 5) then
               v3(i,j,ke+1) =       vg3(ke+1)
               v3(i,j,ke+2) = 2.0 * vg3(ke+1) - v3(i,j,ke)
             endif
           enddo
         enddo
         bvstat(6,5) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               v3(i,j,ks  ) =       vg3(ks)
               v3(i,j,ks-1) = 2.0 * vg3(ks) - v3(i,j,ks+1)
               v3(i,j,ks-2) = 2.0 * vg3(ks) - v3(i,j,ks+2)
             endif
             if (nikb(i,j) .eq. 2) then
               q1           = sign ( haf, v3(i,j,ks+1) - vg3(ks) )
               v3(i,j,ks  ) = v3(i,j,ks+1) * ( 0.5 - q1 )
               v3(i,j,ks-1) = v3(i,j,ks  )
               v3(i,j,ks-2) = v3(i,j,ks  )
             endif
             if (nikb(i,j) .eq. 3) then
               v3(i,j,ks  ) = v3ikb (i,j,1)
               v3(i,j,ks-1) = v3ikb (i,j,2)
               v3(i,j,ks-2) = 2.0 * v3ikb (i,j,2) - v3ikb (i,j,1)
             endif
             if (nikb(i,j) .eq. 4) then
               v3(i,j,ks-1) = v3(i,j,ke  )
               v3(i,j,ks-2) = v3(i,j,ke-1)
             endif
             if (nikb(i,j) .eq. 5) then
               v3(i,j,ks  ) =       vg3(ks)
               v3(i,j,ks-1) = 2.0 * vg3(ks) - v3(i,j,ks+1)
               v3(i,j,ks-2) = 2.0 * vg3(ks) - v3(i,j,ks+2)
             endif
           enddo
         enddo
         bvstat(5,5) = rl3
         endif
!
!      Outer k boundary.
!
         if (u3 .gt. 0) then
         do j=js-1,je+1
!dir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               v3(i,j,ke+1) =       vg3(ke+1)
               v3(i,j,ke+2) = 2.0 * vg3(ke+1) - v3(i,j,ke  )
             endif
             if (nokb(i,j) .eq. 2) then
               q1           = sign ( haf, v3(i,j,ke) - vg3(ke+1) )
               v3(i,j,ke+1) = v3(i,j,ke  ) * ( 0.5 + q1 )
               v3(i,j,ke+2) = v3(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               v3(i,j,ke+1) = v3okb (i,j,1)
               v3(i,j,ke+2) = v3okb (i,j,2)
             endif
             if (nokb(i,j) .eq. 4) then
               v3(i,j,ke+1) = v3(i,j,ks  )
               v3(i,j,ke+2) = v3(i,j,ks+1)
             endif
             if (nokb(i,j) .eq. 5) then
               v3(i,j,ke+1) =       vg3(ke+1)
               v3(i,j,ke+2) = 2.0 * vg3(ke+1) - v3(i,j,ke)
             endif
           enddo
         enddo
         bvstat(6,5) = ru3
         endif
#endif /* NO MPI */
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L V 3                \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L B 1                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalb1 ( rl1, ru1, rl2, ru2, rl3, ru3, b1 )
!
!    dac:zeus3d.bvalb1 <----------- 1-direction 
!                    magnetic field boundary values 
!    from mln:zeus04.bflo; jms:zeus2d.bvalb1              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 1-component of 
!  the magnetic filed.  The active zones for "b1" are "is+1" to "ie" 
!  in the 1-direction, "js" toc  "je" in the 2-direction, and "ks"  
!  to "ke" in the 3-direction.  Two layers of boundary values at each 
!  face are needed for third order interpolation.  No edge or corner 
!  boundary values are required. Thus, the ranges for the boundary 
!  value calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks  , ke
!    j-boundaries:   i = is+1, ie                      k = ks  , ke
!    k-boundaries:   i = is+1, ie     j = js  , je
!
!  Note that for periodic or tile-tile boundaries, "is" is also active.
!
!  The flow out boundary uses a switch to ensure fluid can only flow OUT
!  of the i boundary (boundary value set to 0 if it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3 activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array b1 is input so that velocity values at old time levels
!  as well as momentum components can be passed.
!
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
      integer :: i,j,k,l1,l2,l3,u1,u2,u3
      integer :: rl1,rl2,rl3,ru1,ru2,ru3
#ifdef MPI_USED &
      integer :: ls,ll,lu,us,ul,uu
#endif
      real(rl) :: b1(in,jn,kn)
      real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,3)
       u1 = ru1 - bvstat(2,3)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB1: ',i2,' to receive ',i2,' l i-slices from ' &
             ,i2)") myid, ls, n1m
#endif /* DEBUG */
           do 410 i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b1(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                         ,3100+25*i+nsub, comm3d, req(nreq), ierr)
410        continue
           bvstat(1,3) = rl1
         endif
         if (u1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB1: ',i2,' to send ',i2,' u i-slices to ' &
             ,i2)") myid, us, n1m
#endif /* DEBUG */
           do 420 i=1,us
             nreq = nreq + 1
             call MPI_ISEND(b1(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,3200+25*i+nsub, comm3d, req(nreq), ierr)
420        continue
           bvstat(2,3) = ru1
         endif
       else
#endif
         if (l1 .gt. 0) then
         do 430 k=ks-1,ke+1
!dir$ ivdep
           do 425 j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               b1(is  ,j,k) =       vg1(is)
               b1(is-1,j,k) = 2.0 * vg1(is) - b1(is+1,j,k)
               b1(is-2,j,k) = 2.0 * vg1(is) - b1(is+2,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               q1           = sign ( haf, b1(is+1,j,k) - vg1(is) )
               b1(is  ,j,k) = b1(is+1,j,k) * ( 0.5 - q1 )
               b1(is-1,j,k) = b1(is  ,j,k)
               b1(is-2,j,k) = b1(is  ,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               b1(is  ,j,k) = v1iib (j,k,1)
               b1(is-1,j,k) = v1iib (j,k,2)
               b1(is-2,j,k) = 2.0 * b1iib (j,k,2) - b1iib (j,k,1)
             endif
#ifndef MPI_USED
             if (niib(j,k) .eq. 4) then
               b1(is-1,j,k) = b1(ie  ,j,k)
               b1(is-2,j,k) = b1(ie-1,j,k)
             endif
#endif
             if (niib(j,k) .eq. 5) then
               b1(is  ,j,k) =       vg1(is)
               b1(is-1,j,k) = 2.0 * vg1(is) - b1(is+1,j,k)
               b1(is-2,j,k) = 2.0 * vg1(is) - b1(is+2,j,k)
             endif
425        continue
430      continue
         bvstat(1,3) = rl1
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer i boundary.
!
! 1-face-centered quantities need to be evolved on only one end
! (we have chosen to evolve them on the inner boundary), hence
! no change in subscripts compared to zone-centered quantities.
!
#ifdef MPI_USED
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB1: ',i2,' to receive ',i2,' u i-slices from ' &
             ,i2)") myid, us, n1p
#endif /* DEBUG */
           do 440 i=1,us
             nreq = nreq + 1
             call MPI_IRECV(b1(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,3200+25*i+nsub, comm3d, req(nreq), ierr)
440        continue
           bvstat(2,3) = ru1
         endif
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB1: ',i2,' to send ',i2,' l i-slices to ' &
             ,i2)") myid, ls, n1p
#endif /* DEBUG */
           do 450 i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b1(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,3100+25*i+nsub, comm3d, req(nreq), ierr)
450        continue
           bvstat(1,3) = rl1
         endif
       else
#endif
         if (u1 .gt. 0) then
         do 460 k=ks-1,ke+1
!dir$ ivdep
           do 455 j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               b1(ie+1,j,k) =       vg1(ie+1)
               b1(ie+2,j,k) = 2.0 * vg1(ie+1) - b1(ie  ,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               q1           = sign ( haf, b1(ie,j,k) - vg1(ie+1) )
               b1(ie+1,j,k) = b1(ie,j,k) * ( 0.5 + q1 )
               b1(ie+2,j,k) = b1(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               b1(ie+1,j,k) = b1oib (j,k,1)
               b1(ie+2,j,k) = b1oib (j,k,2)
             endif
#ifndef MPI_USED
             if (noib(j,k) .eq. 4) then
               b1(ie+1,j,k) = b1(is  ,j,k)
               b1(ie+2,j,k) = b1(is+1,j,k)
             endif
#endif
             if (noib(j,k) .eq. 5) then
               b1(ie+1,j,k) =       vg1(ie+1)
               b1(ie+2,j,k) = 2.0 * vg1(ie+1) - b1(ie,j,k)
             endif
455        continue
460      continue
         bvstat(2,3) = ru1
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,3)
       u2 = ru2 - bvstat(4,3)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do 470 j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b1(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                           ,3300+25*j+nsub, comm3d, req(nreq), ierr)
470        continue
           bvstat(3,3) = rl2
         endif
         if (u2 .gt. 0) then
           do 480 j=1,us
             nreq = nreq + 1
             call MPI_ISEND(b1(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                           ,3400+25*j+nsub, comm3d, req(nreq), ierr)
480        continue
           bvstat(4,3) = ru2
         endif
       else
#endif
         if (l2 .gt. 0) then
         do 490 k=ks-1,ke+1
!dir$ ivdep
           do 485 i=is-1,ie+1
             if ( abs(nijb1(i,k)) .eq. 1) then
               b1(i,js-1,k) = b1(i,js  ,k)
               b1(i,js-2,k) = b1(i,js+1,k)
             endif
             if (nijb1(i,k) .eq. 2) then
               b1(i,js-1,k) = b1(i,js  ,k)
               b1(i,js-2,k) = b1(i,js-1,k)
             endif
             if (nijb1(i,k) .eq. 3) then
               b1(i,js-1,k) = b1ijb (i,k,1)
               b1(i,js-2,k) = b1ijb (i,k,2)
             endif
#ifndef MPI_USED
             if (nijb1(i,k) .eq. 4) then
               b1(i,js-1,k) = b1(i,je  ,k)
               b1(i,js-2,k) = b1(i,je-1,k)
             endif
#endif
             if (nijb1(i,k) .eq. 5) then
               b1(i,js-1,k) = b1(i,js  ,k)
               b1(i,js-2,k) = b1(i,js+1,k)
             endif
485        continue
490      continue
         bvstat(3,3) = rl2
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer j boundary.
!
#ifdef MPI_USED
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do 500 j=1,us
             nreq = nreq + 1
             call MPI_IRECV(b1(   1,je+j+uu,   1), 1, j_slice, n2p &
                           ,3400+25*j+nsub, comm3d, req(nreq), ierr)
500        continue
           bvstat(4,3) = ru2
         endif
         if (l2 .gt. 0) then
           do 510 j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b1(   1,je+j-ll,   1), 1, j_slice, n2p &
                           ,3300+25*j+nsub, comm3d, req(nreq), ierr)
510        continue
           bvstat(3,3) = rl2
         endif
       else
#endif
         if (u2 .gt. 0) then
         do 520 k=ks-1,ke+1
!dir$ ivdep
           do 515 i=is-1,ie+1
             if ( abs(nojb1(i,k)) .eq. 1) then
               b1(i,je+1,k) = b1(i,je  ,k)
               b1(i,je+2,k) = b1(i,je-1,k)
             endif
             if (nojb1(i,k) .eq. 2) then
               b1(i,je+1,k) = b1(i,je  ,k)
               b1(i,je+2,k) = b1(i,je+1,k)
             endif
             if (nojb1(i,k) .eq. 3) then
               b1(i,je+1,k) = b1ojb (i,k,1)
               b1(i,je+2,k) = b1ojb (i,k,2)
             endif
#ifndef MPI_USED
             if (nojb1(i,k) .eq. 4) then
               b1(i,je+1,k) = b1(i,js  ,k)
               b1(i,je+2,k) = b1(i,js+1,k)
             endif
#endif
             if (nojb1(i,k) .eq. 5) then
               b1(i,je+1,k) = b1(i,je  ,k)
               b1(i,je+2,k) = b1(i,je-1,k)
             endif
515        continue
520      continue
         bvstat(4,3) = ru2
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,3)
       u3 = ru3 - bvstat(6,3)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do 530 k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b1(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,3500+25*k+nsub, comm3d, req(nreq), ierr)
530        continue
           bvstat(5,3) = rl3
         endif
         if (u3 .gt. 0) then
           do 540 k=1,us
             nreq = nreq + 1
             call MPI_ISEND(b1(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,3600+25*k+nsub, comm3d, req(nreq), ierr)
540        continue
           bvstat(6,3) = ru3
         endif
       else
#endif
       if (l3 .gt. 0) then
       do 550 j=js-1,je+1
!dir$ ivdep
           do 545 i=is-1,ie+1
             if ( abs(nikb1(i,j)) .eq. 1) then
               b1(i,j,ks-1) = b1(i,j,ks  )
               b1(i,j,ks-2) = b1(i,j,ks+1)
             endif
             if (nikb1(i,j) .eq. 2) then
               b1(i,j,ks-1) = b1(i,j,ks  )
               b1(i,j,ks-2) = b1(i,j,ks-1)
             endif
             if (nikb1(i,j) .eq. 3) then
               b1(i,j,ks-1) = b1ikb (i,j,1)
               b1(i,j,ks-2) = b1ikb (i,j,2)
             endif
#ifndef MPI_USED
             if (nikb1(i,j) .eq. 4) then
               b1(i,j,ks-1) = b1(i,j,ke  )
               b1(i,j,ks-2) = b1(i,j,ke-1)
             endif
#endif
             if (nikb1(i,j) .eq. 5) then
               b1(i,j,ks-1) = b1(i,j,ks  )
               b1(i,j,ks-2) = b1(i,j,ks+1)
             endif
545        continue
550      continue
         bvstat(5,3) = rl3
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer k boundary.
!
#ifdef MPI_USED
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do 560 k=1,us
             nreq = nreq + 1
             call MPI_IRECV(b1(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,3600+25*k+nsub, comm3d, req(nreq), ierr)
560        continue
           bvstat(6,3) = ru3
         endif
         if (l3 .gt. 0) then
           do 570 k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b1(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,3500+25*k+nsub, comm3d, req(nreq), ierr)
570        continue
           bvstat(5,3) = rl3
         endif
       else
#endif
         if (u3 .gt. 0) then
         do 580 j=js-1,je+1
!dir$ ivdep
           do 575 i=is-1,ie+1
             if ( abs(nokb1(i,j)) .eq. 1) then
               b1(i,j,ke+1) = b1(i,j,ke  )
               b1(i,j,ke+2) = b1(i,j,ke-1)
             endif
             if (nokb1(i,j) .eq. 2) then
               b1(i,j,ke+1) = b1(i,j,ke  )
               b1(i,j,ke+2) = b1(i,j,ke+1)
             endif
             if (nokb1(i,j) .eq. 3) then
               b1(i,j,ke+1) = b1okb (i,j,1)
               b1(i,j,ke+2) = b1okb (i,j,2)
             endif
#ifndef MPI_USED
             if (nokb1(i,j) .eq. 4) then
               b1(i,j,ke+1) = b1(i,j,ks  )
               b1(i,j,ke+2) = b1(i,j,ks+1)
             endif
#endif
             if (nokb1(i,j) .eq. 5) then
               b1(i,j,ke+1) = b1(i,j,ke  )
               b1(i,j,ke+2) = b1(i,j,ke-1)
             endif
575        continue
580      continue
         bvstat(6,3) = ru3
         endif
#ifdef MPI_USED
       endif
#endif
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L V 1                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L B 2                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalb2 ( rl1, ru1, rl2, ru2, rl3, ru3, b2 )
!
!    dac:zeus3d.bvalb2 <----------- 2-direction velocity boundary values
!    from mln:zeus04.bflo; jms:zeus2d.bvalb2              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 2-component of the 
!  magnetic field. The active zones for "v2" are "is" to "ie" in the 
!  1-direction, "js+1" to "je" in the 2-direction, and "ks" to "ke" 
!  in the 3-direction. Two layers of boundary values at each face are 
!  needed for third order interpolation.  No edge or corner boundary 
!  values are required. Thus, the ranges for the boundary value 
!  calculations are:
!
!    i-boundaries:                    j = js+1, je     k = ks  , ke
!    j-boundaries:   i = is  , ie                      k = ks  , ke
!    k-boundaries:   i = is  , ie     j = js+1, je
!
!  Note that for periodic or tile-tile boundaries, "js" is also active.
!
!  The flow out boundary uses a switch to ensure magnetic field can 
!  can only go OUT of the j boundary (boundary value set to 0 if 
!  it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array b2 is input so that velocity values at old time levels
!  and momenta can be passed.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
       integer :: i,j,k,l1,l2,l3,u1,u2,u3
       integer :: rl1,rl2,rl3,ru1,ru2,ru3
#ifdef MPI_USED &
       integer :: ls,ll,lu,us,ul,uu
#endif
       real(rl) :: b2(in,jn,kn)
       real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,4)
       u1 = ru1 - bvstat(2,4)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB2: ',i2,' to receive ',i2,' l i-slices from ' &
             ,i2)") myid, ls, n1m
#endif /* DEBUG */
           do 610 i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b2(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,4100+25*i+nsub, comm3d, req(nreq), ierr)
610        continue
           bvstat(1,4) = rl1
         endif
         if (u1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB2: ',i2,' to send ',i2,' u i-slices to ' &
               ,i2)") myid, us, n1m
#endif /* DEBUG */
           do 620 i=1,us
             nreq = nreq + 1
             call MPI_ISEND(b2(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,4200+25*i+nsub, comm3d, req(nreq), ierr)
620        continue
           bvstat(2,4) = ru1
         endif
       else
#endif
!
         if (l1 .gt. 0) then
         do 630 k=ks-1,ke+1
!dir$ ivdep
           do 625 j=js-1,je+1
             if (niib2(j,k) .eq. 1) then
               b2(is-1,j,k) = b2(is  ,j,k)
               b2(is-2,j,k) = b2(is+1,j,k)
             endif
             if (niib2(j,k) .eq.-1) then
#ifdef RTP
               b2(is-1,j,k) =-b2(is  ,j,k)
               b2(is-2,j,k) =-b2(is+1,j,k)
#else
               b2(is-1,j,k) = b2(is  ,j,k)
               b2(is-2,j,k) = b2(is+1,j,k)
#endif
             endif
             if (niib2(j,k) .eq. 2) then
               b2(is-1,j,k) = b2(is  ,j,k)
               b2(is-2,j,k) = b2(is-1,j,k)
             endif
             if (niib2(j,k) .eq. 3) then
               b2(is-1,j,k) = b2iib (j,k,1)
               b2(is-2,j,k) = b2iib (j,k,2)
             endif
#ifndef MPI_USED
             if (niib2(j,k) .eq. 4) then
               b2(is-1,j,k) = b2(ie  ,j,k)
               b2(is-2,j,k) = b2(ie-1,j,k)
             endif
#endif
             if (niib2(j,k) .eq. 5) then
               b2(is-1,j,k) = b2(is  ,j,k)
               b2(is-2,j,k) = b2(is+1,j,k)
             endif
625        continue
630      continue
         bvstat(1,4) = rl1
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer i boundary.
!
#ifdef MPI_USED
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
#ifdef DEBUG
        write(*,"('BVALB2: ',i2,' to receive ',i2,' u i-slices from ' &
             ,i2)") myid, us, n1p
#endif /* DEBUG */
           do 640 i=1,us
             nreq = nreq + 1
             call MPI_IRECV(b2(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,4200+25*i+nsub, comm3d, req(nreq), ierr)
640        continue
           bvstat(2,4) = ru1
         endif
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB2: ',i2,' to send ',i2,' l i-slices to ' &
             ,i2)") myid, ls, n1p
#endif /* DEBUG */
           do 650 i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b2(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,4100+25*i+nsub, comm3d, req(nreq), ierr)
650        continue
           bvstat(1,4) = rl1
         endif
       else
#endif
         if (u1 .gt. 0) then
         do 660 k=ks-1,ke+1
!dir$ ivdep
           do 655 j=js-1,je+1
             if ( abs(noib2(j,k)) .eq. 1) then
               b2(ie+1,j,k) = b2(ie,j,k)
               b2(ie+2,j,k) = b2(ie-1,j,k)
             endif
             if (noib2(j,k) .eq. 2) then
!#if PROBLEM == advect  --   I wish cpp could do this!
!               b2(ie+1,j,k) = 2.0 * b2(ie  ,j,k) - b2(ie-1,j,k)
!               b2(ie+2,j,k) = 2.0 * b2(ie+1,j,k) - b2(ie  ,j,k)
!#else
               b2(ie+1,j,k) = b2(ie,j,k)
               b2(ie+2,j,k) = b2(ie+1,j,k)
!#endif
             endif
             if (noib2(j,k) .eq. 3) then
               b2(ie+1,j,k) = b2oib (j,k,1)
               b2(ie+2,j,k) = b2oib (j,k,2)
             endif
#ifndef MPI_USED
             if (noib2(j,k) .eq. 4) then
               b2(ie+1,j,k) = b2(is  ,j,k)
               b2(ie+2,j,k) = b2(is+1,j,k)
             endif
#endif
             if (noib2(j,k) .eq. 5) then
               b2(ie+1,j,k) = b2(ie  ,j,k)
               b2(ie+2,j,k) = b2(ie-1,j,k)
             endif
655        continue
660      continue
         bvstat(2,4) = ru1
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,4)
       u2 = ru2 - bvstat(4,4)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do 670 j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b2(   1,js-ll+j-1,   1), 1 , j_slice, n2m &
                           ,4300+25*j+nsub, comm3d, req(nreq), ierr)
670        continue
           bvstat(3,4) = rl2
         endif
         if (u2 .gt. 0) then
           do 680 j=1,us
             nreq = nreq + 1
             call MPI_ISEND(b2(   1,js+uu+j-1,   1), 1 , j_slice, n2m &
                           ,4400+25*j+nsub, comm3d, req(nreq), ierr)
680        continue
           bvstat(4,4) = ru2
         endif
       else
#endif
         if (l2 .gt. 0) then
         do 690 k=ks-1,ke+1
!dir$ ivdep
           do 685 i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               b2(i,js  ,k) =       vg2(js)
               b2(i,js-1,k) = 2.0 * vg2(js) - b2(i,js+1,k)
               b2(i,js-2,k) = 2.0 * vg2(js) - b2(i,js+2,k)
             endif
             if (nijb(i,k) .eq. 2) then
               q1           = sign ( haf, b2(i,js+1,k) - vg2(js) )
               b2(i,js  ,k) = b2(i,js+1,k) * ( 0.5 - q1 )
               b2(i,js-1,k) = b2(i,js  ,k)
               b2(i,js-2,k) = b2(i,js  ,k)
             endif
             if (nijb(i,k) .eq. 3) then
               b2(i,js  ,k) = b2ijb (i,k,1)
               b2(i,js-1,k) = b2ijb (i,k,2)
               b2(i,js-2,k) = 2.0 * v2ijb (i,k,2) - v2ijb (i,k,1)
             endif
#ifndef MPI_USED
             if (nijb(i,k) .eq. 4) then
               b2(i,js-1,k) = b2(i,je  ,k)
               b2(i,js-2,k) = b2(i,je-1,k)
             endif
#endif
             if (nijb(i,k) .eq. 5) then
               b2(i,js  ,k) =       vg2(js)
               b2(i,js-1,k) = 2.0 * vg2(js) - b2(i,js+1,k)
               b2(i,js-2,k) = 2.0 * vg2(js) - b2(i,js+2,k)
             endif
685        continue
690      continue
         bvstat(3,4) = rl2
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer j boundary.
!
#ifdef MPI_USED
      if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do 700 j=1,us
             nreq = nreq + 1
             call MPI_IRECV(b2(   1,je+j+uu,   1), 1 , j_slice, n2p &
                           ,4400+25*j+nsub, comm3d, req(nreq), ierr)
700        continue
           bvstat(4,4) = ru2
         endif
         if (l2 .gt. 0) then
           do 710 j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b2(   1,je+j-ll,   1), 1 , j_slice, n2p &
                           ,4300+25*j+nsub, comm3d, req(nreq), ierr)
710        continue
           bvstat(3,4) = rl2
         endif
       else
#endif
         if (u2 .gt. 0) then
         do 720 k=ks-1,ke+1
!dir$ ivdep
           do 715 i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               b2(i,je+1,k) =       vg2(je+1)
               b2(i,je+2,k) = 2.0 * vg2(je+1) - b2(i,je  ,k)
             endif
             if (nojb(i,k) .eq. 2) then
               q1           = sign ( haf, b2(i,je,k) - vg2(je+1) )
               b2(i,je+1,k) = b2(i,je  ,k) * ( 0.5 + q1 )
               b2(i,je+2,k) = b2(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               b2(i,je+1,k) = b2ojb (i,k,1)
               b2(i,je+2,k) = b2ojb (i,k,2)
             endif
#ifndef MPI_USED
             if (nojb(i,k) .eq. 4) then
               b2(i,je+1,k) = b2(i,js  ,k)
               b2(i,je+2,k) = b2(i,js+1,k)
             endif
#endif
             if (nojb(i,k) .eq. 5) then
               b2(i,je+1,k) =       vg2(je+1)
               b2(i,je+2,k) = 2.0 * vg2(je+1) - b2(i,je  ,k)
             endif
715        continue
720      continue
         bvstat(4,4) = ru2
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,4)
       u3 = ru3 - bvstat(6,4)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do 730 k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b2(   1,   1,ks-ll+k-1), 1, k_slice, n3m &
                           ,4500+25*k+nsub, comm3d, req(nreq), ierr)
730        continue
           bvstat(5,4) = rl3
         endif
         if (u3 .gt. 0) then
           do 740 k=1,us
             nreq = nreq + 1
             call MPI_ISEND(b2(   1,   1,ks+uu+k-1), 1, k_slice, n3m &
                           ,4600+25*k+nsub, comm3d, req(nreq), ierr)
740        continue
           bvstat(6,4) = ru3
         endif
       else
#endif
         if (l3 .gt. 0) then
         do 750 j=js-1,je+1
!dir$ ivdep
           do 745 i=is-1,ie+1
             if ( abs(nikb2(i,j)) .eq. 1) then
               b2(i,j,ks-1) = b2(i,j,ks  )
               b2(i,j,ks-2) = b2(i,j,ks+1)
             endif
             if (nikb2(i,j) .eq. 2) then
               b2(i,j,ks-1) = b2(i,j,ks  )
               b2(i,j,ks-2) = b2(i,j,ks-1)
             endif
             if (nikb2(i,j) .eq. 3) then
               b2(i,j,ks-1) = b2ikb (i,j,1)
               b2(i,j,ks-2) = b2ikb (i,j,2)
             endif
#ifndef MPI_USED
             if (nikb2(i,j) .eq. 4) then
               b2(i,j,ks-1) = b2(i,j,ke  )
               b2(i,j,ks-2) = b2(i,j,ke-1)
             endif
#endif
             if (nikb2(i,j) .eq. 5) then
               b2(i,j,ks-1) = b2(i,j,ks  )
               b2(i,j,ks-2) = b2(i,j,ks+1)
             endif
745        continue
750      continue
         bvstat(5,4) = rl3
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer k boundary.
!
#ifdef MPI_USED
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do 760 k=1,us
             nreq = nreq + 1
             call MPI_IRECV(b2(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,4600+25*k+nsub, comm3d, req(nreq), ierr)
760        continue
           bvstat(6,4) = ru3
         endif
         if (l3 .gt. 0) then
           do 770 k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b2(   1,   1,ke+k-ll), 1, k_slice, n3p &
                           ,4500+25*k+nsub, comm3d, req(nreq), ierr)
770        continue
           bvstat(5,4) = rl3
         endif
       else
#endif
         if (u3 .gt. 0) then
         do 780 j=js-1,je+1
!dir$ ivdep
           do 775 i=is-1,ie+1
             if ( abs(nokb2(i,j)) .eq. 1) then
               b2(i,j,ke+1) = b2(i,j,ke  )
               b2(i,j,ke+2) = b2(i,j,ke-1)
             endif
             if (nokb2(i,j) .eq. 2) then
               b2(i,j,ke+1) = b2(i,j,ke  )
               b2(i,j,ke+2) = b2(i,j,ke+1)
             endif
             if (nokb2(i,j) .eq. 3) then
               b2(i,j,ke+1) = b2okb (i,j,1)
               b2(i,j,ke+2) = b2okb (i,j,2)
             endif
#ifndef MPI_USED
             if (nokb2(i,j) .eq. 4) then
               b2(i,j,ke+1) = b2(i,j,ks  )
               b2(i,j,ke+2) = b2(i,j,ks+1)
             endif
#endif
             if (nokb2(i,j) .eq. 5) then
               b2(i,j,ke+1) = b2(i,j,ke  )
               b2(i,j,ke+2) = b2(i,j,ke-1)
             endif
775        continue
780      continue
         bvstat(6,4) = ru3
         endif
#ifdef MPI_USED
       endif
#endif
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L B 2                \\\\\\\\\\
!
!=======================================================================
!
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                B V A L B 3                \\\\\\\\\\
!
!=======================================================================
!
       subroutine bvalb3 ( rl1, ru1, rl2, ru2, rl3, ru3, b3 )
!
!    dac:zeus3d.bvalv3 <----------- 3-direction velocity boundary values
!    from mln:zeus04.bflo; jms:zeus2d.bvalv3              february, 1990
!
!    written by: David Clarke, February 1990
!    modified 1: RAF, 3/5/96 for ZEUS-MP
!
!  PURPOSE: This routine sets boundary values for the 3-component of the 
!  magnetic field. The active zones for "v3" are "is" to "ie" in the 
!  1-direction, "js" to "je" in the 2-direction, and "ks+1" to "ke" 
!  in the 3-direction.  Two layers of boundary values at each face 
!  are needed for interpolation. No edge or corner boundary values 
!  are required.  Thus, the ranges for the boundary value 
!  calculations are:
!
!    i-boundaries:                    j = js  , je     k = ks+1, ke
!    j-boundaries:   i = is  , ie                      k = ks+1, ke
!    k-boundaries:   i = is  , ie     j = js  , je
!
!  Note that for periodic or tile-tile boundaries, "ks" is also active.
!
!  The flow out boundary uses a switch to ensure fluid can only flow OUT
!  of the k boundary (boundary value set to 0 if it tries to flow in).
!
!  See comments in BVALD.
!
!  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
!  Their values give the number of layers to pass.
!
!  Array b3 is input so that velocity values at old time levels
!  can be passed.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
       integer :: i,j,k,l1,l2,l3,u1,u2,u3
       integer :: rl1,rl2,rl3,ru1,ru2,ru3
#ifdef MPI_USED &
       integer :: ls,ll,lu,us,ul,uu
#endif
       real(rl) :: b3(in,jn,kn)
       real(rl) :: q1
!
!-----------------------------------------------------------------------
!------------------------  I - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l1 = rl1 - bvstat(1,5)
       u1 = ru1 - bvstat(2,5)
!
!      Inner i boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
!
! Post a receive for a slab of data from the interior of the 
! neighboring tile to fill my ghost zones.  Initiate a send 
! to pass a slab of my interior data for my neighbor's ghost zones.
! 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB3: ',i2,' to receive ',i2,' l i-slices from ' &
             ,i2)") myid, ls, n1m
#endif /* DEBUG */
           do 810 i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b3(is-ll+i-1,   1,   1), 1, i_slice, n1m &
                           ,5100+25*i+nsub, comm3d, req(nreq), ierr)
810        continue
           bvstat(1,5) = rl1
         endif
         if (u1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB3: ',i2,' to send ',i2,' u i-slices to ' &
             ,i2)") myid, us, n1m
#endif /* DEBUG */
           do 820 i=1,us
             nreq = nreq + 1
             call MPI_ISEND(b3(is+uu+i-1,   1,   1), 1, i_slice, n1m &
                           ,5200+25*i+nsub, comm3d, req(nreq), ierr)
820        continue
           bvstat(2,5) = ru1
         endif
       else
#endif
         if (l1 .gt. 0) then
         do 830 k=ks-1,ke+1
!dir$ ivdep
           do 825 j=js-1,je+1
             if (niib3(j,k) .eq. 1) then
               b3(is-1,j,k) = b3(is  ,j,k)
               b3(is-2,j,k) = b3(is+1,j,k)
             endif
             if (niib3(j,k) .eq.-1) then
#ifdef RTP
               b3(is-1,j,k) =-b3(is  ,j,k)
               b3(is-2,j,k) =-b3(is+1,j,k)
#else
               b3(is-1,j,k) = b3(is  ,j,k)
               b3(is-2,j,k) = b3(is+1,j,k)
#endif
             endif
             if (niib3(j,k) .eq. 2) then
               b3(is-1,j,k) = b3(is  ,j,k)
               b3(is-2,j,k) = b3(is-1,j,k)
             endif
             if (niib3(j,k) .eq. 3) then
               b3(is-1,j,k) = b3iib (j,k,1)
               b3(is-2,j,k) = b3iib (j,k,2)
             endif
#ifndef MPI_USED
             if (niib3(j,k) .eq. 4) then
               b3(is-1,j,k) = b3(ie  ,j,k)
               b3(is-2,j,k) = b3(ie-1,j,k)
             endif
#endif
             if (niib3(j,k) .eq. 5) then
               b3(is-1,j,k) = b3(is  ,j,k)
               b3(is-2,j,k) = b3(is+1,j,k)
             endif
825        continue
830      continue
         bvstat(1,5) = rl1
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer i boundary.
!
#ifdef MPI_USED
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB3: ',i2,' to receive ',i2,' u i-slices from ' &
             ,i2)") myid, us, n1p
#endif /* DEBUG */
           do 840 i=1,us
             nreq = nreq + 1
             call MPI_IRECV(b3(ie+i+uu,   1,   1), 1, i_slice, n1p &
                           ,5200+25*i+nsub, comm3d, req(nreq), ierr)
840        continue
           bvstat(2,5) = ru1
         endif
         if (l1 .gt. 0) then
#ifdef DEBUG
           write(*,"('BVALB3: ',i2,' to send ',i2,' l i-slices to ' &
             ,i2)") myid, ls, n1p
#endif /* DEBUG */
           do 850 i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b3(ie+i-ll,   1,   1), 1, i_slice, n1p &
                           ,5100+25*i+nsub, comm3d, req(nreq), ierr)
850        continue
           bvstat(1,5) = rl1
         endif
       else
#endif
         if (u1 .gt. 0) then
         do 860 k=ks-1,ke+1
!dir$ ivdep
           do 855 j=js-1,je+1
             if ( abs(noib3(j,k)) .eq. 1) then
               b3(ie+1,j,k) = b3(ie,j,k)
               b3(ie+2,j,k) = b3(ie-1,j,k)
             endif
             if (noib3(j,k) .eq. 2) then
!#if PROBLEM == advect  --  I wish!
!               b3(ie+1,j,k) = 2.0 * b3(ie  ,j,k) - b3(ie-1,j,k)
!               b3(ie+2,j,k) = 2.0 * b3(ie+1,j,k) - b3(ie  ,j,k)
!#else
               b3(ie+1,j,k) = b3(ie,j,k)
               b3(ie+2,j,k) = b3(ie+1,j,k)
!#endif
             endif
             if (noib3(j,k) .eq. 3) then
               b3(ie+1,j,k) = b3oib (j,k,1)
               b3(ie+2,j,k) = b3oib (j,k,2)
             endif
#ifndef MPI_USED
             if (noib3(j,k) .eq. 4) then
               b3(ie+1,j,k) = b3(is  ,j,k)
               b3(ie+2,j,k) = b3(is+1,j,k)
             endif
#endif
             if (noib3(j,k) .eq. 5) then
               b3(ie+1,j,k) = b3(ie  ,j,k)
               b3(ie+2,j,k) = b3(ie-1,j,k)
             endif
855        continue
860      continue
         bvstat(2,5) = ru1
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  J - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l2 = rl2 - bvstat(3,5)
       u2 = ru2 - bvstat(4,5)
!
!      Inner j boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
!
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do 870 j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b3(   1,js-ll+j-1,   1), 1, j_slice, n2m &
                           ,5300+25*j+nsub, comm3d, req(nreq), ierr)
870        continue
           bvstat(3,5) = rl2
         endif
         if (u2 .gt. 0) then
           do 880 j=1,us
             nreq = nreq + 1
             call MPI_ISEND(b3(   1,js+uu+j-1,   1), 1, j_slice, n2m &
                           ,5400+25*j+nsub, comm3d, req(nreq), ierr)
880        continue
           bvstat(4,5) = ru2
         endif
       else
#endif
         if (l2 .gt. 0) then
         do 890 k=ks-1,ke+1
!dir$ ivdep
           do 885 i=is-1,ie+1
             if (nijb3(i,k) .eq. 1) then
               b3(i,js-1,k) = b3(i,js  ,k)
               b3(i,js-2,k) = b3(i,js+1,k)
             endif
             if (nijb3(i,k) .eq.-1) then
#ifdef XYZ
               b3(i,js-1,k) = b3(i,js  ,k)
               b3(i,js-2,k) = b3(i,js+1,k)
#else
               b3(i,js-1,k) =-b3(i,js  ,k)
               b3(i,js-2,k) =-b3(i,js+1,k)
#endif
             endif
             if (nijb3(i,k) .eq. 2) then
               b3(i,js-1,k) = b3(i,js  ,k)
               b3(i,js-2,k) = b3(i,js-1,k)
             endif
             if (nijb3(i,k) .eq. 3) then
               b3(i,js-1,k) = b3ijb (i,k,1)
               b3(i,js-2,k) = b3ijb (i,k,2)
             endif
#ifndef MPI_USED
             if (nijb3(i,k) .eq. 4) then
               b3(i,js-1,k) = b3(i,je  ,k)
               b3(i,js-2,k) = b3(i,je-1,k)
             endif
#endif
             if (nijb3(i,k) .eq. 5) then
               b3(i,js-1,k) = b3(i,js  ,k)
               b3(i,js-2,k) = b3(i,js+1,k)
             endif
885      continue
890    continue
       bvstat(3,5) = rl2
       endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer j boundary.
!
#ifdef MPI_USED
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do 900 j=1,us
             nreq = nreq + 1
             call MPI_IRECV(b3(   1,je+j+uu,   1), 1, j_slice, n2p &
                           ,5400+25*j+nsub, comm3d, req(nreq), ierr)
900        continue
           bvstat(4,5) = ru2
         endif
         if (l2 .gt. 0) then
           do 910 j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b3(   1,je+j-ll,   1), 1, j_slice, n2p &
                           ,5300+25*j+nsub, comm3d, req(nreq), ierr)
910        continue
           bvstat(3,5) = rl2
         endif
       else
#endif
         if (u2 .gt. 0) then
         do 920 k=ks-1,ke+1
!dir$ ivdep
           do 915 i=is-1,ie+1
             if (nojb3(i,k) .eq. 1) then
               b3(i,je+1,k) = b3(i,je  ,k)
               b3(i,je+2,k) = b3(i,je-1,k)
             endif
             if (nojb3(i,k) .eq.-1) then
#ifdef RTP
               b3(i,je+1,k) =-b3(i,je  ,k)
               b3(i,je+2,k) =-b3(i,je-1,k)
#else
               b3(i,je+1,k) = b3(i,je  ,k)
               b3(i,je+2,k) = b3(i,je-1,k)
#endif
             endif
             if (nojb3(i,k) .eq. 2) then
!#if alias PROBLEM.eq.advect  -- I wish!
               b3(i,je+1,k) = 2.0 * b3(i,je  ,k) - b3(i,je-1,k)
               b3(i,je+2,k) = 2.0 * b3(i,je+1,k) - b3(i,je  ,k)
!#else
               b3(i,je+1,k) = b3(i,je  ,k)
               b3(i,je+2,k) = b3(i,je+1,k)
!#endif
             endif
             if (nojb3(i,k) .eq. 3) then
               b3(i,je+1,k) = b3ojb (i,k,1)
               b3(i,je+2,k) = b3ojb (i,k,2)
             endif
#ifndef MPI_USED
             if (nojb3(i,k) .eq. 4) then
               b3(i,je+1,k) = b3(i,js  ,k)
               b3(i,je+2,k) = b3(i,js+1,k)
             endif
#endif
             if (nojb3(i,k) .eq. 5) then
               b3(i,je+1,k) = b3(i,je  ,k)
               b3(i,je+2,k) = b3(i,je-1,k)
             endif
915        continue
920      continue
         bvstat(4,5) = ru2
         endif
#ifdef MPI_USED
       endif
#endif
!
!-----------------------------------------------------------------------
!------------------------  K - B O U N D A R Y  ------------------------
!-----------------------------------------------------------------------
!
       l3 = rl3 - bvstat(5,5)
       u3 = ru3 - bvstat(6,5)
!
!      Inner k boundary.
!
#ifdef MPI_USED
!
! Count slabs, compute positioning indices.
!
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
!
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do 930 k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(b3(   1,   1,ks-ll+k-1), 1 , k_slice, n3m &
                           ,5500+25*k+nsub, comm3d, req(nreq), ierr)
930        continue
           bvstat(5,5) = rl3
         endif
         if (u3 .gt. 0) then
           do 940 k=1,us
             nreq = nreq + 1
             call MPI_ISEND(b3(   1,   1,ks+uu+k-1), 1 , k_slice, n3m &
                           ,5600+25*k+nsub, comm3d, req(nreq), ierr)
940        continue
           bvstat(6,5) = ru3
         endif
       else
#endif
         if (l3 .gt. 0) then
         do 950 j=js-1,je+1
!dir$ ivdep
           do 945 i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               b3(i,j,ks  ) =       vg3(ks)
               b3(i,j,ks-1) = 2.0 * vg3(ks) - b3(i,j,ks+1)
               b3(i,j,ks-2) = 2.0 * vg3(ks) - b3(i,j,ks+2)
             endif
             if (nikb(i,j) .eq. 2) then
               q1           = sign ( haf, b3(i,j,ks+1) - vg3(ks) )
               b3(i,j,ks  ) = b3(i,j,ks+1) * ( 0.5 - q1 )
               b3(i,j,ks-1) = b3(i,j,ks  )
               b3(i,j,ks-2) = b3(i,j,ks  )
             endif
             if (nikb(i,j) .eq. 3) then
               b3(i,j,ks  ) = b3ikb (i,j,1)
               b3(i,j,ks-1) = b3ikb (i,j,2)
               b3(i,j,ks-2) = 2.0 * v3ikb (i,j,2) - v3ikb (i,j,1)
             endif
#ifndef MPI_USED
             if (nikb(i,j) .eq. 4) then
               b3(i,j,ks-1) = b3(i,j,ke  )
               b3(i,j,ks-2) = b3(i,j,ke-1)
             endif
#endif
             if (nikb(i,j) .eq. 5) then
               b3(i,j,ks  ) =       vg3(ks)
               b3(i,j,ks-1) = 2.0 * vg3(ks) - b3(i,j,ks+1)
               b3(i,j,ks-2) = 2.0 * vg3(ks) - b3(i,j,ks+2)
             endif
945        continue
950      continue
         bvstat(5,5) = rl3
         endif
#ifdef MPI_USED
       endif
#endif
!
!      Outer k boundary.
!
#ifdef MPI_USED
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do 960 k=1,us
             nreq = nreq + 1
             call MPI_IRECV(b3(   1,   1,ke+k+uu), 1, k_slice, n3p &
                           ,5600+25*k+nsub, comm3d, req(nreq), ierr)
960        continue
           bvstat(6,5) = ru3
         endif
         if (l3 .gt. 0) then
           do 970 k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(b3(   1,   1,ke+k-ll), 1 , k_slice, n3p &
                           ,5500+25*k+nsub, comm3d, req(nreq), ierr)
970        continue
           bvstat(5,5) = rl3
         endif
       else
#endif
         if (u3 .gt. 0) then
         do 980 j=js-1,je+1
!dir$ ivdep
           do 975 i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               b3(i,j,ke+1) =       vg3(ke+1)
               b3(i,j,ke+2) = 2.0 * vg3(ke+1) - b3(i,j,ke  )
             endif
             if (nokb(i,j) .eq. 2) then
               q1           = sign ( haf, b3(i,j,ke) - vg3(ke+1) )
               b3(i,j,ke+1) = b3(i,j,ke  ) * ( 0.5 + q1 )
               b3(i,j,ke+2) = b3(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               b3(i,j,ke+1) = b3okb (i,j,1)
               b3(i,j,ke+2) = b3okb (i,j,2)
             endif
#ifndef MPI_USED
             if (nokb(i,j) .eq. 4) then
               b3(i,j,ke+1) = b3(i,j,ks  )
               b3(i,j,ke+2) = b3(i,j,ks+1)
             endif
#endif
             if (nokb(i,j) .eq. 5) then
               b3(i,j,ke+1) =       vg3(ke+1)
               b3(i,j,ke+2) = 2.0 * vg3(ke+1) - b3(i,j,ke)
             endif
975        continue
980      continue
         bvstat(6,5) = ru3
         endif
#ifdef MPI_USED
       endif
#endif
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                B V A L B 3                \\\\\\\\\\
!
!=======================================================================
!
!
!
