!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               O F F D I A G               \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
!
      subroutine offdiag(ddp1,ddp2,ddp3)
!
!     evaluates all off-diagonal elements in the jacobian for the
!     solution to the implicit radiation energy equation.  all
!     terms in the sub- and super-diagonals depend only on quantities
!     which are constant during the n-r iteration, and thus only
!     need computing at the start.
!
!     written by: John C. Hayes; March 1997
!     Modified by: John Hayes, 2003; rewritten in F90 for Version 2
!
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use radiation
      use opac
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
!     loop indices and loop bounds
!
      integer  :: i, j, k
!
!     matrix arrays
!
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
!
!     scalar work variables
!
      real(rl) :: d_f1_d_em1, d_f1_d_ep1, &
                  d_f2_d_em1, d_f2_d_ep1, &
                  d_f3_d_em1, d_f3_d_ep1
!
      real(rl) :: d_deldotf_d_eip1, d_deldotf_d_eim1, &
                  d_deldotf_d_ejp1, d_deldotf_d_ejm1, &
                  d_deldotf_d_ekp1, d_deldotf_d_ekm1
!
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         ddp1(1,1,i,j,k) = 0.D0
         ddp2(1,1,i,j,k) = 0.D0
         ddp3(1,1,i,j,k) = 0.D0
        enddo
       enddo
      enddo
!
!     terms for the first super- and sub-diagonals
!
      do k = ks, ke
       do j = js, je
        do i = is, ie
         ddp1(1,1,i,j,k)  = -radth*dt*dvl2a(j)*dvl3a(k)* &
                            (g2a(i+1)*g31a(i+1))**2 * &
                            dr1(i+1,j  ,k  )*dvl1bi(i+1)
         if(ldimen .gt. 1) then
          ddp2(1,1,i,j,k)  = -radth*dt*dvl1a(i)*dvl3a(k)* &
                             (g2bi(i  )*g32a(j+1))**2 * &
                              dr2(i  ,j+1,k  )*dvl2bi(j+1)
          if(ldimen .eq. 3) then
           ddp3(1,1,i,j,k)  = -radth*dt*dvl1a(i)*dvl2a(j)* &
                              (g31bi(i  )*g32bi(j  ))**2 * &
                               dr3(i  ,j  ,k+1)*dvl3bi(k+1)
          endif
         endif
        enddo
       enddo
      enddo
!
!======================================================================
!     compute off-diagonal elements at physical domain boundaries,
!     if necessary
!======================================================================
!
!     Outer boundary of 1-coordinate, so zero out first super-diagonal
!     over the 2-3 outer coordinate face.  Other diagonals are 
!     unaffected.
!
      if(coords(1) .eq. ntiles(1)-1) then
       do k = ks, ke
        do j = js, je
         ddp1(1,1,ie,j,k) = 0.D0
        enddo
       enddo
      endif
      if(ldimen .gt. 1) then
       if(coords(2) .eq. ntiles(2)-1) then
        do k = ks, ke
         do i = is, ie
          ddp2(1,1,i,je,k) = 0.D0
         enddo
        enddo
       endif
       if(ldimen .eq. 3) then
        if(coords(3) .eq. ntiles(3)-1) then
         do j = js, je
          do i = is, ie
           ddp3(1,1,i,j,ke) = 0.D0
          enddo
         enddo
        endif
       endif ! ldimen = 3
      endif ! ldimen > 1
!
!======================================================================
!    Fill ghost zones when periodic B.C.'s are used and MPI is not
!    used.
!======================================================================
!
#ifndef MPI_USED
      do k = ks, ke
       do j = js, je
        if(liib(j,k) .eq. 4) then
         ddp1(1,1,is-1,j,k) = ddp1(1,1,ie  ,j,k)
         ddp1(1,1,is-2,j,k) = ddp1(1,1,ie-1,j,k)
         ddp1(1,1,ie+1,j,k) = ddp1(1,1,is  ,j,k)
         ddp1(1,1,ie+2,j,k) = ddp1(1,1,is+1,j,k)
        endif
       enddo
      enddo
      if(ldimen .gt. 1) then
       do k = ks, ke
        do i = is, ie
         if(lijb(i,k) .eq. 4) then
          ddp2(1,1,i,js-1,k) = ddp2(1,1,i,je  ,k)
          ddp2(1,1,i,js-2,k) = ddp2(1,1,i,je-1,k)
          ddp2(1,1,i,je+1,k) = ddp2(1,1,i,js  ,k)
          ddp2(1,1,i,je+2,k) = ddp2(1,1,i,js+1,k)
         endif
        enddo
       enddo
       if(ldimen .gt. 2) then
        do j = js, je
         do i = is, ie
          if(likb(i,j) .eq. 4) then
           ddp3(1,1,i,j,ks-1) = ddp3(1,1,i,j,ke  )
           ddp3(1,1,i,j,ks-2) = ddp3(1,1,i,j,ke-1)
           ddp3(1,1,i,j,ke+1) = ddp3(1,1,i,j,ks  )
           ddp3(1,1,i,j,ke+2) = ddp3(1,1,i,j,ks+1)
          endif
         enddo
        enddo
       endif ! ldimen > 2
      endif ! ldimen > 1
#endif /* NO MPI */
!
      return
      end
