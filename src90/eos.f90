!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                    E O S                  \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine eos(ibeg,iend,jbeg,jend,kbeg,kend)
!
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
      use eos_par
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, &
                  kbeg, kend
!
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         if(leos .eq. 1) then
          if(.not. xiso) then
           if(.not. xtotnrg) then
            p(i,j,k) = gamm1*e(i,j,k)
           else
            p(i,j,k) = gamm1*( &
                       e(i,j,k) - 0.125D0*d(i,j,k)*( &
                        (v1(i  ,j  ,k  )+v1(i+1,j  ,k  ))**2 &
                       +(v2(i  ,j  ,k  )+v2(i  ,j+1,k  ))**2 &
                       +(v3(i  ,j  ,k  )+v3(i  ,j  ,k+1))**2 ) &
                              )
           endif
          else
           p(i,j,k) = ciso**2*d(i,j,k)
          endif
         endif ! IDEAL GAS
        enddo ! i
       enddo ! j
      enddo ! k
!
      return
      end
