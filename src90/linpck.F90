!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      integer function isamax(n,sx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
      implicit NONE
      real*8 sx(*),smax
      integer i,incx,ix,n
!
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
!
      integer function ismax(n,sx,incx)
!
!     finds the index of element having max. value.
!     jack dongarra, linpack, 3/11/78. (added to this file by j. stone)
!
      implicit NONE
      real*8 sx(*),smax
      integer i,incx,ix,n
!
      ismax = 0
      if( n .lt. 1 ) return
      ismax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = sx(1)
      ix = ix + incx
      do 10 i = 2,n
         if(sx(ix).le.smax) go to 5
         ismax = i
         smax = sx(ix)
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = sx(1)
      do 30 i = 2,n
         if(sx(i).le.smax) go to 30
         ismax = i
         smax = sx(i)
   30 continue
      return
      end
      integer function ismin(n,sx,incx)
!
!     finds the index of element having min. value.
!     jack dongarra, linpack, 3/11/78.  (added to this file by j. stone)
!
      implicit NONE
      real*8 sx(*),smin
      integer i,incx,ix,n
!
      ismin = 0
      if( n .lt. 1 ) return
      ismin = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smin = sx(1)
      ix = ix + incx
      do 10 i = 2,n
         if(sx(ix).ge.smin) go to 5
         ismin = i
         smin = sx(ix)
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smin = sx(1)
      do 30 i = 2,n
         if(sx(i).ge.smin) go to 30
         ismin = i
         smin = sx(i)
   30 continue
      return
      end
      real*8 function sasum(n,sx,incx)
!
!     takes the sum of the absolute values.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      implicit NONE
      real*8 sx(*),stemp
      integer i,incx,m,mp1,n,nincx
!
      sasum = 0.0e0
      stemp = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(sx(i))
   10 continue
      sasum = stemp
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + abs(sx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2)) &
       + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
   60 sasum = stemp
      return
      end
