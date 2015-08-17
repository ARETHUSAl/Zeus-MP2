!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 N E W X 3                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
      subroutine newx3
!
!  PURPOSE: Computes "new" x3 grid variables (grid variables at advanced
!  timestep) to be used in TRANSPRT.  Grid values are calculated for
!  k=ks-2 to ke+2, except for dvl3a (k=ks,ke+2) and dvl3b (k=ks+1,ke+2).
!
!  EXTERNALS: [none]
!
!  LOCALS:
!   vol3an,vol3bn = volume factors used to compute dvl3*n
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use root
      use scratch
!
      implicit NONE
!
      integer  :: k
      real(rl) :: vol3an(kn) ,  vol3bn(kn), qa,qb,qc,qd
!
!=======================================================================
!
      x3an(ks-2) = x3a(ks-2) + vg3(ks-2)*dt
      do 10 k=ks-1,ke+2
         x3an(k  ) = x3a (k) + vg3(k)*dt
        dx3an(k-1) = x3an(k) - x3an(k-1)
10    continue
      dx3an(ke+2) = (dx3an(ke+1)/dx3an(ke)) * dx3an(ke+1)
!
      dx3bn(ks-2) = dx3an(ks-2)
       x3bn(ks-2) =  x3an(ks-1) - 0.5*dx3an(ks-2)
      do 20 k=ks-1,ke+2
         x3bn(k) = x3an(k) + 0.5*dx3an(k)
        dx3bn(k) = x3bn(k) - x3bn(k-1)
20    continue
!
!  New volume factors
!
      vol3an(ks-2) = x3an(ks-2)
      do 40 k=ks-2,ke+1
        vol3an(k+1) = x3an(k+1)
        dvl3an(k  ) = vol3an(k+1) - vol3an(k)
40    continue
!
      vol3bn(ks-2) = x3bn(ks-2)
      do 50 k=ks-2,ke+1
        vol3bn(k+1) = x3bn(k+1)
        dvl3bn(k+1) = vol3bn(k+1) - vol3bn(k)
50    continue
      do k = ks-2, ke+1
       dvl3ani(k) = 1.0D0/dvl3an(k)
       dvl3bni(k) = 1.0D0/dvl3bn(k)
      enddo
!
      return
      end
