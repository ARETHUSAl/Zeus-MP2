!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                   S O D                   \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
!
       subroutine sod
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
!
      implicit NONE
!
      integer  :: i  , j, k, ip1, jp1, &
                  kp1, idirect
!
      real(rl) :: x10, x20, x30, &
                  d0 , p0 , e0 , &
                  d1 , p1 , e1
!
      integer  :: iin (ijkn), iout(ijkn), jin (ijkn), &
                  jout(ijkn), kin (ijkn), kout(ijkn)
!
      real(rl) :: massk(ijkn)
!
      namelist / pgen     / &
                     x10, x20, p0, d0, p1, d1, idirect, x30
!
!-----------------------------------------------------------------------
!
#ifdef IGNORE
       do k = 1, kn
       do j = 1, jn
       do i = 1, in
        v1(i,j,k) = 1.0D10
       enddo
       enddo
       enddo
#endif
       x10  = 0.5D0
       x20  = 0.5D0
       x30  = 0.5D0
       d0   = 1.0D0
       p0   = 1.0D0
       d1   = 0.125D0
       p1   = 0.1D0
!
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = x10 
         buf_in(2) = x20 
         buf_in(3) = d0  
         buf_in(4) = p0  
         buf_in(5) = d1  
         buf_in(6) = p1  
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 6, MPI_FLOAT &
                     , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 1, MPI_INTEGER &
                     , 0, comm3d, ierr )
       if (myid .ne. 0) then
         x10  = buf_in(1)
         x20  = buf_in(2)
         d0   = buf_in(3)
         p0   = buf_in(4)
         d1   = buf_in(5)
         p1   = buf_in(6)
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
!
!      Set up tube.
!
      e0 = p0 / gamm1
      e1 = p1 / gamm1
      do k=1,kn
       do j=1,jn
        do i=1,in
         v1(i,j,k) = 0.0D0
         v2(i,j,k) = 0.0D0
         v3(i,j,k) = 0.0D0
        enddo
       enddo
      enddo
      if(.false.) then
      do k = 1, kn
       do j = 1, jn
        v1oib(j,k,1) = -1.0D0
        v1oib(j,k,2) = -1.0D0
        doib(j,k,1) = d1
        doib(j,k,2) = d1
        eoib(j,k,1) = e1
        eoib(j,k,2) = e1
       enddo
      enddo
      endif
      if(.false.) then
      do k = 1, kn
       do i = 1, in
        v2ojb(i,k,1) = -1.0D0
        v2ojb(i,k,2) = -1.0D0
        dojb(i,k,1) = d1
        dojb(i,k,2) = d1
        eojb(i,k,1) = e1
        eojb(i,k,2) = e1
       enddo
      enddo
      endif
      if(.false.) then
      do j = 1, jn
       do i = 1, in
        v3okb(i,j,1) = -1.0D0
        v3okb(i,j,2) = -1.0D0
        dokb(i,j,1) = d1
        dokb(i,j,2) = d1
        eokb(i,j,1) = e1
        eokb(i,j,2) = e1
       enddo
      enddo
      endif
!
      if(idirect .eq. 1) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x1a(i) .le. x10) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
      if(idirect .eq. 2) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x2a(j) .le. x20) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
      if(idirect .eq. 3) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x3a(k) .le. x30) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
!
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                   S O D                   \\\\\\\\\\
!
!=======================================================================
!
!
      subroutine sodres
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
!
      implicit NONE
!
      integer  :: i, j, k, ip1, jp1, kp1, idirect
      real(rl) :: x10, x20, x30, &
                  d0 , p0 , e0 , &
                  d1 , p1 , e1
!
      integer  :: iin (ijkn), iout(ijkn), jin (ijkn), &
                  jout(ijkn), kin (ijkn), kout(ijkn)
      real(rl) :: massk(ijkn)
!
      namelist / pgen     / &
                     x10, x20, p0, d0, p1, d1, idirect
!
!-----------------------------------------------------------------------
!
       x10  = 0.5D0
       x20  = 0.5D0
       x30  = 0.5D0
       d0   = 1.0D0
       p0   = 1.0D0
       d1   = 0.125D0
       p1   = 0.1D0
!
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = x10 
         buf_in(2) = x20 
         buf_in(3) = d0  
         buf_in(4) = p0  
         buf_in(5) = d1  
         buf_in(6) = p1  
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 6, MPI_FLOAT &
                     , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 1, MPI_INTEGER &
                     , 0, comm3d, ierr )
       if (myid .ne. 0) then
         x10  = buf_in(1)
         x20  = buf_in(2)
         d0   = buf_in(3)
         p0   = buf_in(4)
         d1   = buf_in(5)
         p1   = buf_in(6)
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
!
       return
       end
