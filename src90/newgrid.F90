!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               N E W G R I D               \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
      subroutine newgrid
!
!  PURPOSE: Controls grid motion by computing grid velocities and "new"
!  grid variables (variables at advanced time, to be used in TRANSPRT).
!
!  EXTERNALS: NEWVG
!             NEWX1
!             NEWX2
!             NEWX3
!
!  modified for F90 version of ZEUS-MP by J. Hayes, 5-2003
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use bndry
      use field
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif /* MPI_USED */
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k
!
!=======================================================================
!  return if there is no grid motion in any direction
!
      if ((x1fac .eq. 0.0) .and. (x2fac .eq. 0.0) .and. &
          (x3fac .eq. 0.0)) return
!
!  update "X1" grid
!
      if (x1fac .ne. 0.0) then
       do i = 1, 6
        bvstat(1,i) = 0
        bvstat(3,i) = 0
       enddo
       nreq = 0
       nsub = nsub + 1
       call bvald (3,3,0,0,0,0, d)
       call bvalv1(3,3,0,0,0,0,v1)
!
        call scopy (in,  x1an ,1,  x1a ,1)
        call scopy (in,  x1bn ,1,  x1b ,1)
        call scopy (in, dx1an ,1, dx1a ,1)
        call scopy (in, dx1bn ,1, dx1b ,1)
        call scopy (in, g2an  ,1, g2a  ,1)
        call scopy (in, g2bn  ,1, g2b  ,1)
        call scopy (in, g31an ,1, g31a ,1)
        call scopy (in, g31bn ,1, g31b ,1)
        call scopy (in,dvl1an ,1,dvl1a ,1)
        call scopy (in,dvl1bn ,1,dvl1b ,1)
!
        do 10 i=is-2,ie+2
         x1ai  (i) = 1.0/(x1a  (i)+tiny)
         x1bi  (i) = 1.0/(x1b  (i)+tiny)
         dx1ai (i) = 1.0/(dx1a (i)+tiny)
         dx1bi (i) = 1.0/(dx1b (i)+tiny)
         g2ai  (i) = 1.0/(g2a  (i)+tiny)
         g2bi  (i) = 1.0/(g2b  (i)+tiny)
         g31ai (i) = 1.0/(g31a (i)+tiny)
         g31bi (i) = 1.0/(g31b (i)+tiny)
         dvl1ai(i) = 1.0/(dvl1a(i)+tiny)
         dvl1bi(i) = 1.0/(dvl1b(i)+tiny)
10      continue
!
!  update grid velocities
!
#ifdef MPI_USED
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
#endif /* MPI_USED */
       call newvg
!
       call newx1
!
       do i = is-2, ie+2
        dx1ani(i) = 1.0D0/dx1an(i)
        dx1bni(i) = 1.0D0/dx1bn(i)
       enddo
!
       do k = 1, kn
        do j = 1, jn
         v3oib(j,k,1) = 0.3*x1b(ie+1)*sin(x2b(j))
         v3oib(j,k,2) = 0.3*x1b(ie+2)*sin(x2b(j))
        enddo
       enddo
      endif ! x1fac
!
!  update "X2" grid 
!
      if (x2fac .ne. 0.0) then
        call scopy (jn,  x2an ,1,  x2a ,1)
        call scopy (jn,  x2bn ,1,  x2b ,1)
        call scopy (jn, dx2an ,1, dx2a ,1)
        call scopy (jn, dx2bn ,1, dx2b ,1)
        call scopy (jn, g32an ,1, g32a ,1)
        call scopy (jn, g32bn ,1, g32b ,1)
        call scopy (jn, g4an  ,1, g4a  ,1)
        call scopy (jn, g4bn  ,1, g4b  ,1)
        call scopy (jn,dvl2an ,1,dvl2a ,1)
        call scopy (jn,dvl2bn ,1,dvl2b ,1)
        call scopy (jn,dvl2ani,1,dvl2ai,1)
        call scopy (jn,dvl2bni,1,dvl2bi,1)
        do 20 j=js-2,je+2
         dx2ai(j) = 1.0/dx2a(j)
         dx2bi(j) = 1.0/dx2b(j)
         if(lgeom .eq. 3) then
          dg32ad2(j) = cos(x2a(j))
          dg32bd2(j) = cos(x2b(j))
         endif
20      continue
        call newx2
      endif ! x2fac
!
      if(x3fac .ne. 0)then
        call scopy (kn,  x3an ,1,  x3a ,1)
        call scopy (kn,  x3bn ,1,  x3b ,1)
        call scopy (kn, dx3an ,1, dx3a ,1)
        call scopy (kn, dx3bn ,1, dx3b ,1)
        call scopy (kn,dvl3an ,1,dvl3a ,1)
        call scopy (kn,dvl3bn ,1,dvl3b ,1)
        call scopy (kn,dvl3ani,1,dvl3ai,1)
        call scopy (kn,dvl3bni,1,dvl3bi,1)
        do 30 k=ks-2,ke+2
         dx3ai(k) = 1.0/dx3a(k)
         dx3bi(k) = 1.0/dx3b(k)
30      continue
        call newx3
      endif ! x3fac
!
      return
      end
!
      subroutine  scopy(n,sx,incx,sy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!
      use real_prec
      use param
!
      implicit NONE
!
      real(rl) :: sx(*),sy(*)
      integer  :: i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
