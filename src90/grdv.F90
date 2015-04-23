!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 G R D V                   \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine grdv (ibeg,iend,jbeg,jend,kbeg,kend, &
                       erold,gvcf,deldotv)
!
!======================================================================
!     Worker routine for computing the tensor inner product of grad(v)
!     with P, the radiation energy stress tensor, and with f, the
!     Eddington tensor.  In our application, the only difference 
!     P and f is a factor of er(i,j,k), the central radiation energy
!     density.  This routine is basically a modified copy of Robert
!     Fiedler's NR_PF routine, with unneeded portions deleted.
!
!======================================================================
!
!     AUTHORS
!       Robert Fiedler and John Hayes
!
!     LAST MODIFIED
!       5/26/2003
!
!     INPUT
!       ibeg, iend   The range of i-values to cover.
!       jbeg, jend   The range of j-values to cover.
!       kbeg, kend   The range of k-values to cover.
!
!     OUTPUT
!
!       gvcf     grad(v) : f
!       gvcp     grad(v) : P
!======================================================================
!
      use real_prec
      use config
      use param
      use cons
      use field
      use grid
      use radiation
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use root
      use opac
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, jbeg, kbeg, iend, jend, kend
!
      real(rl) :: et11, et22, et33, et12, et13, et23, ef, grad_er, &
                  derdx1, derdx2, derdx3, ros_mfp, flx_lim, big_R, temp
!
      real(rl) :: dv11(in), dv12(in), dv13(in), dv21(in), dv22(in), &
                  dv23(in), dv31(in), dv32(in), dv33(in)
!
      real(rl) :: gvcf   (in,jn,kn), gvcp(in,jn,kn), erold(in,jn,kn), &
                  deldotv(in,jn,kn)
!
      real(rl) :: onethird
!
!.......................................................................
!
!     The following computation of the inner product of 
!     grad(v) and the Eddington tensor should be performed 
!     only once per time step, just before the Newton-Raphson
!     iterations begin.  
!.......................................................................
!
      if(.not. xhydro) then
       do k = kbeg, kend
        do j = jbeg, jend
         do i = ibeg, iend
          gvcf   (i,j,k) = 0.0D0
          deldotv(i,j,k) = 0.0D0
         enddo
        enddo
       enddo
       return
      endif ! no hydro
      onethird = 1.D0 / 3.D0
!
      do k = kbeg, kend
       do j = jbeg, jend
!
!      Compute grad(v) and save the 9 elements in 1-D arrays.  
!
!      The following terms are not averaged over i when computing the 
!      inner product.
!
        do i = ibeg, iend
         dv11(i) = (v1(i+1,j  ,k  ) - v1(i  ,j  ,k  )) * dx1ai(i)
         dv22(i) = &
                   0.5 * (v1(i+1,j  ,k  ) + v1(i  ,j  ,k  )) &
                 * g2bi(i) * dg2bd1(i) 
         dv33(i) = &
                   0.5 * (v1(i+1,j  ,k  ) + v1(i  ,j  ,k  )) &
                 * g31bi(i) * dg31bd1(i) 
        enddo ! i
        if(ldimen .gt. 1) then
         do i = ibeg, iend
          dv22(i) = dv22(i) &
                  + (v2(i  ,j+1,k  ) - v2(i  ,j  ,k  )) * dx2ai(j) &
                  * g2bi(i)
          dv33(i) = dv33(i) &
                  + 0.5 * (v2(i  ,j+1,k  ) + v2(i  ,j  ,k  )) &
                  * g32bi(j) * g2bi(i) * dg32bd2(j) &
                  + (v3(i  ,j  ,k+1) - v3(i  ,j  ,k  )) * dx3ai(k) &
                  * g31bi(i) * g32bi(j)
         enddo ! i
        endif ! ldimen > 1
        do i = ibeg, iend
         if(ldimen .gt. 1) then
          dv23(i) = (dx2bi(j  ) * (v3(i,j  ,k  ) - v3(i,j-1,k  )) &
                  +  dx2bi(j+1) * (v3(i,j+1,k  ) - v3(i,j  ,k  )) &
                  +  dx2bi(j  ) * (v3(i,j  ,k+1) - v3(i,j-1,k+1)) &
                  +  dx2bi(j+1) * (v3(i,j+1,k+1) - v3(i,j  ,k+1))) &
                  * g2bi(i) * 0.25 
         else
          dv23(i) = 0.0D0
         endif ! ldimen > 1
         if(ldimen .gt. 1) then
          dv32(i) = ((dx3bi(k  ) * (v2(i,j  ,k  ) - v2(i,j  ,k-1)) &
                  +   dx3bi(k+1) * (v2(i,j  ,k+1) - v2(i,j  ,k  ))) &
                  *  g32ai(j  ) &
                  +  (dx3bi(k  ) * (v2(i,j+1,k  ) - v2(i,j+1,k-1)) &
                  +   dx3bi(k+1) * (v2(i,j+1,k+1) - v2(i,j+1,k  ))) &
                  *  g32ai(j+1)) &
                  * g31bi(i) * 0.25 &
                  - (0.5 * (v3(i,j  ,k  ) + v3(i,j-1,k  )) &
                  * g32ai(j  ) * dg32ad2(j  ) &
                  +  0.5 * (v3(i,j+1,k  ) + v3(i,j  ,k  )) &
                  * g32ai(j+1) * dg32ad2(j+1) &
                  +  0.5 * (v3(i,j  ,k+1) + v3(i,j-1,k+1)) &
                  * g32ai(j  ) * dg32ad2(j  ) &
                  +  0.5 * (v3(i,j+1,k+1) + v3(i,j  ,k+1)) &
                  * g32ai(j+1) * dg32ad2(j+1)) &
                  * g2bi(i) * 0.25
         else
          dv32(i) = 0.0D0
         endif
        enddo ! i
!
!.......................................................................
!   These terms are averaged over i when computing the inner product.
!.......................................................................
!
        do i = ibeg, iend+1
         if(ldimen .gt. 1) then
          dv12(i) = 0.5 * dx1bi(i) * (v2(i  ,j  ,k) - v2(i-1,j  ,k) &
                                   +  v2(i  ,j+1,k) - v2(i-1,j+1,k))
          dv13(i) = 0.5 * dx1bi(i) * (v3(i  ,j,k  ) - v3(i-1,j,k  ) &
                                   +  v3(i  ,j,k+1) - v3(i-1,j,k+1))
         else
          dv12(i) = 0.0D0
          dv13(i) = 0.0D0
         endif
         dv21(i) = 0.5 * (dx2bi(j  ) * (v1(i,j  ,k) - v1(i,j-1,k)) &
                       +  dx2bi(j+1) * (v1(i,j+1,k) - v1(i,j  ,k))) &
                 * g2ai(i)
         dv31(i) = 0.5 * (dx3bi(k  ) * (v1(i,j,k  ) - v1(i,j,k-1)) &
                       +  dx3bi(k+1) * (v1(i,j,k+1) - v1(i,j,k  ))) &
                 * g31ai(i) * g32bi(j)
        enddo ! i
        if(ldimen .gt. 1) then
         do i = ibeg, iend+1
          dv21(i) = dv21(i) &
                  - 0.25 * (v2(i  ,j  ,k) + v2(i-1,j  ,k) &
                         +  v2(i  ,j+1,k) + v2(i-1,j+1,k)) &
                  * g2ai(i) * dg2ad1(i) 
          dv31(i) = dv31(i) &
                  - 0.25 * (v3(i  ,j,k  ) + v3(i-1,j,k  ) &
                         +  v3(i  ,j,k+1) + v3(i-1,j,k+1)) &
                  * g31ai(i) * dg31ad1(i)
         enddo ! i
        endif ! ldimen > 1
!
! Compute the Eddington tensor, a function of the Eddington factor ef.
!
        do i = ibeg, iend
!
!   Get the zone-centered gradient of the radiation energy density er.
!
         derdx1 = (erold(i+1,j  ,k  ) - erold(i-1,j  ,k  )) &
                / ( dx1b(i  ) + dx1b(i+1) )
         if(ldimen .eq. 1) then
          derdx2 = 0.0
         else
          derdx2 = (erold(i  ,j+1,k  ) - erold(i  ,j-1,k  ))*g2bi(i) &
                 / ( dx2b(j  ) + dx2b(j+1) )
         endif ! ldimen = 1
         if(ldimen .ne. 3) then
          derdx3 = 0.0
         else
          derdx3 = (erold(i  ,j  ,k+1) - erold(i  ,j  ,k-1)) &
                 * g31bi(i) * g32bi(i) &
                  / ( dx3b(k  ) + dx3b(k+1) )
         endif ! ldimen /= 3
         grad_er = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
!
!   Get the Rosseland mean free path at zone centers.
!
         ros_mfp = 1.0 / (kapr(i,j,k) + sig(i,j,k))
!
!   Get the dimensionless quantity R that the flux limiter depends on.
!
         big_R  = grad_er * ros_mfp &
               / max(erold(i,j,k), tiny)
!
!   Get the flux limiter flx_lim
#ifdef MINERBO
!
!     Minerbo flux limiter
!
         if (big_R .le. 1.5) then
          flx_lim = 2.0/(3.0 + sqrt(9.0 + 12.0*big_R**2) )
         else
          flx_lim = 1.0/(1.0 + big_R + sqrt(1.0 + 2.0*big_R) )
         endif
#else /* MINERBO */
!
!     LP flux limiter
!
         flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
#endif /* MINERBO */
!
         ef = flx_lim + flx_lim**2 * big_R**2
!
!   Renormalize grad er.
!
         derdx1 = derdx1 / max ( grad_er, tiny )
         derdx2 = derdx2 / max ( grad_er, tiny )
         derdx3 = derdx3 / max ( grad_er, tiny )
!
!   Get the zone-centered Eddington tensor components.
!   Write them carefully, so that ef - 1/3 is exactly 0
!   for pure diffusion.
!
         et11 = 0.5 * (1.0 - ef) &
              + 1.5 * (ef - onethird) * derdx1**2
         et22 = 0.5 * (1.0 - ef) &
              + 1.5 * (ef - onethird) * derdx2**2
         et33 = 0.5 * (1.0 - ef) &
              + 1.5 * (ef - onethird) * derdx3**2
         et12 = 1.5 * (ef - onethird) * derdx1 * derdx2
         et13 = 1.5 * (ef - onethird) * derdx1 * derdx3
         et23 = 1.5 * (ef - onethird) * derdx2 * derdx3
!
! Take the inner product of the Eddington tensor with the gradient of 
! the velocity vector.  This is the grad(v):P term except for the
! factor of er (at the new time level).
!
         gvcf(i,j,k) = dv11(i) * et11 + dv22(i) * et22 &
                     + dv33(i) * et33 &
                     + (dv12(i  ) + dv12(i+1) + dv21(i  ) &
                     +  dv21(i+1)) * 0.5 * et12 &
                     + (dv13(i  ) + dv13(i+1) + dv31(i  ) &
                     +  dv31(i+1)) * 0.5 * et13 &
                     + (dv23(i) + dv32(i)) * et23
!
! velocity divergence term
!
         deldotv(i,j,k) = dvl1ai(i) * &
                         ( g2a(i+1)*g31a(i+1)*v1(i+1,j  ,k  ) - &
                           g2a(i  )*g31a(i  )*v1(i  ,j  ,k  )    ) &
                        + g2bi(i) * dvl2ai(j) * &
                         (          g32a(j+1)*v2(i  ,j+1,k  ) - &
                                    g32a(j  )*v2(i  ,j  ,k  )    ) &
                        + g31bi(i) * g32bi(j) * dvl3ai(k) * &
                         (                    v3(i  ,j  ,k+1) - &
                                              v3(i  ,j  ,k  )    )
        enddo ! i
       enddo ! j
      enddo ! k
!
      return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                  G R D V                  \\\\\\\\\\
!
!=======================================================================
