subroutine gasdev(idum,vv)
! This function returns a normally distributed deviate with
! zero mean and unit variance, using ran1(IDUM) as the
! source of uniform deviates
       implicit none 
       integer, intent(inout) :: idum
       integer                :: iset,n,k,gotoc
       real(8)                :: v1,v2,r,fac,gset,gasdev1
       real(8)                :: temp1, temp2,vv
       integer :: ix1, ix2
       external ran1
       !write(*,*) 'want','random numbers'
       iset = 0
       v1 = 0.
       v2 = 0.
       r = 0.
       fac = 0.
       gset = 0.
       gasdev1 = 0.
       ix1=0.
       ix2=0.
       gotoc = 0
       if (iset.eq.0) then
1         call ran1(idum,temp1)  
          call ran1(idum,temp2)  
          v1 = 2.*temp1-1. 
          v2 = 2.*temp2-1.
!          write(*,*) 'v1',v1
!          write(*,*) 'v2',v2
!          write(*,*) 'idum',idum
          r = v1**2.+v2**2.
          if (r.ge.1.) then 
             !gotoc = gotoc + 1
             !if (gotoc.gt.n) then
             !   write(*,*) 'error in gasdev'
             !end if
             go to 1
          end if
          fac = dsqrt(-2.*dlog(r)/r)
          gset = v1*fac
          gasdev1 = v2*fac
          iset = 1
       else
          gasdev1 = gset
          iset = 0
       end if
       
       vv=gasdev1
!       write(*,*) 'random normal number',k,vv(k)
       return
end subroutine gasdev
