!
! JReppin 01/07/15
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////             B E N C H M A R K             \\\\\\\\\\
!=======================================================================
!
subroutine benchmark
!
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use subs, only: cross, dot, dot2
      implicit none
      integer :: i,j,k
      integer :: ip, jp, kp
!
      character (len=100) :: dummy
      complex :: ii=(0.,1.)
      integer :: jf,nvect,ivect
      real(rl), dimension(:,:,:), allocatable :: a1, a2, a3
      real(rl), dimension (3) :: ee,kk,exk
      real(rl), dimension (in) :: kdotx
      real(rl) :: exk2,phik,ampl,ampl0
      namelist / pgen  / ampl0
!
! set default values
!
      ampl0=1.0e-4
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
        buf_in(1) = ampl0
      endif
!      call MPI_BCAST( ibuf_in, 8, MPI_INTEGER &
!                     , 0, comm3d, ierr )
      call MPI_BCAST(buf_in, 1, MPI_FLOAT &
                     , 0, comm3d, ierr )
      if (myid .ne. 0) then
        ampl0 = buf_in(1)
      endif
!
!
! allocate arrays 
!
      allocate(a1(in,jn,kn))
      allocate(a2(in,jn,kn))
      allocate(a3(in,jn,kn))
      a1(:,:,:) = 0.0
      a2(:,:,:) = 0.0
      a3(:,:,:) = 0.0
      v1(:,:,:) = 0.0
      v2(:,:,:) = 0.0
      v3(:,:,:) = 0.0
      b1(:,:,:) = 0.0
      b2(:,:,:) = 0.0
      b3(:,:,:) = 0.0
      d (:,:,:) = 1.0
!
!  read header
!
      open(10,file='kvect.dat', status='old')
      read(10,*) nvect
      read(10,*) dummy
!
!  loop through all vectors
!
      do ivect=1,nvect
        read(10,*) kk,ee,ampl,phik
        call cross(ee,kk,exk)
        call dot2(exk,exk2)
        exk=exk/sqrt(exk2)
!
        do k=1,kn
          do j=1,jn
            do i=1,in
              kdotx(i) =kk(1)*x1b(i)+kk(2)*x2b(j)+kk(3)*x3b(k)
              a1(i,j,k) = a1(i,j,k) &
                          +ampl0*ampl*exk(1)*real(exp(ii*(kdotx(i)+phik)))
              a2(i,j,k) = a1(i,j,k) &
                          +ampl0*ampl*exk(2)*real(exp(ii*(kdotx(i)+phik)))
              a3(i,j,k) = a1(i,j,k) &
                          +ampl0*ampl*exk(3)*real(exp(ii*(kdotx(i)+phik)))
            enddo
          enddo
        enddo
      enddo
      close(10)
! ------------------------------------------------------------------------
! ---- compute b out of curl A -------------------------------------------
! ------------------------------------------------------------------------
      do k=ks-2,ke+2
        kp = k + 1
        do j=js-2,je+2
          jp = j + 1
          do i=is-2,ie+2
            ip = i + 1
            b1(i,j,k) =  ( a3(i,jp,k) - a3(i,j,k) - &
                           a2(i,j,kp) + a2(i,j,k) )*real(ijkn)
            b2(i,j,k) =  ( a1(i,j,kp) - a1(i,j,k) - &
                           a3(ip,j,k) + a3(i,j,k) )*real(ijkn)
            b3(i,j,k) =  ( a2(ip,j,k) - a2(i,j,k) - &
                           a1(i,jp,k) + a1(i,j,k) )*real(ijkn)
          enddo
        enddo
      enddo
      deallocate(a1,a2,a3)
      if (myid .eq. 0) then
        write(*,*) 'minimum & maximum velocity:'
        write(*,*) 'min(u_x)=',minval(v1), 'max(u_x)=', maxval(v1)
        write(*,*) 'min(u_y)=',minval(v2), 'max(u_y)=', maxval(v2)
        write(*,*) 'min(u_z)=',minval(v3), 'max(u_z)=', maxval(v3)
        write(*,*) 'minimum & maximum magnetic field:'
        write(*,*) 'min(B_x)=',minval(b1), 'max(B_x)=', maxval(b1)
        write(*,*) 'min(B_y)=',minval(b2), 'max(B_y)=', maxval(b2)
        write(*,*) 'min(B_z)=',minval(b3), 'max(B_z)=', maxval(b3)
        write(*,*) 'minimum & maximum densitiy:'
        write(*,*) 'min(rho)=',minval(d),  'max(rho)=', maxval(d)
      endif
end subroutine benchmark
