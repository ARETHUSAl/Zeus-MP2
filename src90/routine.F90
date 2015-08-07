!***********************************************************************
    subroutine random_isotropic_shell(f,jf,ampl0)
!
!   random_isotropic_shell
!
!   25-jun-15/axel: coded
!
      use Sub, only: cross, dot, dot2
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      character (len=1) :: dummy
      complex :: ii=(0.,1.)
      integer :: jf,i,nvect,ivect
      real, dimension (3) :: ee,kk,exk
      real, dimension (nx) :: kdotx
      real :: exk2,phik,ampl,ampl0
!
!  read header
!
      open(10,file='kvect.dat')
      read(10,*) nvect
      read(10,*) dummy
!
!  loop through all vectors
!
      f(:,:,:,jf:jf+2)=0.
      do ivect=1,nvect
        read(10,*) kk,ee,phik,ampl
        call cross(ee,kk,exk)
        call dot2(exk,exk2)
        exk=exk/sqrt(exk2)
!
        do n=n1,n2
          do m=m1,m2
            kdotx=kk(1)*x(l1:l2)+kk(2)*y(m)+kk(3)*z(n)
            do i=1,3
              f(l1:l2,m,n,jf+i-1)=f(l1:l2,m,n,jf+i-1) &
                +ampl0*ampl*exk(i)*real(exp(ii*(kdotx+phik)))
            enddo
          enddo
        enddo
      enddo
      close(10)

    endsubroutine random_isotropic_shell
