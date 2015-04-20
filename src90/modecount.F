      subroutine modecount(modeco,dd,nbox,norm,idx)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!! count modes
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use real_prec
      use param
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      implicit NONE
      real(rl)               :: sigma , modeco , box , window
      real(rl)               :: norm, idx, kmag
      real(rl), dimension(3) :: k, help
      integer                :: dd, J, iz, iy, ix, i, scale, nbox
      integer                :: N
      external sigma
    
      N = ijkn
      modeco = 0.
      box = (float(nbox)/float(N))/2.

      do 30 iz=N/2, N/2+dd-1
        do 20 iy=N/2-dd+1, N/2+dd-1         
          do 10 ix=N/2-dd+1, N/2+dd-1 

            if (iz.eq.N/2) then    ! lies in xy-plane, make sure only half-space
              if (iy.gt.N/2) go to 1
              if (iy.eq.N/2.and.ix.gt.N/2) go to 1                 
              go to 2
            endif 

1           continue

            k(1) = two*pi*float(ix-N/2)          
            k(2) = two*pi*float(iy-N/2)
            k(3) = two*pi*float(iz-N/2)

            kmag = two*pi*sqrt(float(ix-N/2)**2.+          ! L - box size =1
     &              float(iy-N/2)**2.+float(iz-N/2)**2.)


            do 5 i=1, 3
              if (k(i).ne.0.) then
                help(i) = (sin(k(i)*box)/(k(i)*box))
              else
                help(i) = 1.
              endif
5           continue

            window = help(1)*help(2)*help(3)

            modeco = modeco + ((sigma(kmag,idx)/norm)*kmag*window)**2.

2           continue

10        continue ! ix
20      continue ! iy
30    continue ! iz

      modeco = sqrt(2.*modeco)

      return
      end
