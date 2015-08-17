!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                  M N M X                  \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine mnmx ( qty, i1, j1, k1, i2, j2, k2 &
                       , qmin, imin, jmin, kmin &
                       , qmax, imax, jmax, kmax )
!
!    dac:zeus3d.mnmx <------------- finds extrema of a regular 3-d array
!    from mln:zeus04.minmax; jms:zeus2d.minmax            february, 1990
!
!    written by: David Clarke
!    modified 1: converted to ZEUS-MP by Mordecai-Mark Mac Low M
!
!  PURPOSE: This subroutine returns the maximum and minimum value of a
!  rectangular 3-D array, along with the coordinates of the extrema.
!
!  INPUT VARIABLES:
!    qty             the 3-D array to be searched for extrema.
!    i1              inner i index.
!    j1              inner j index.
!    k1              inner k index.
!    i2              outer i index.
!    j2              outer j index.
!    k2              outer k index.
!
!  OUTPUT VARIABLES:
!    qmin            minimum value
!    imin,jmin,kmin  coordinates of minimum value
!    qmax            maximum value
!    imax,jmax,kmax  coordinates of maximum value
!
!  LOCAL VARIABLES:
!
!  j-sweep
!    qmaxj (qminj)   vector of maximum (minimum) values of "qty" from
!                    each i-sweep.  This vector is filled during a
!                    j-sweep.
!    imaxj           i-index of each "qmaxj"
!    iminj           i-index of each "qminj"
!
!  k-sweep
!    qmaxk (qmink)   vector of maximum (minimum) values of "qmaxj"
!                    ("qminj") from each j-sweep.
!    imaxk           i-index of each "qmaxk"
!    jmaxk           j-index of each "qmaxk"
!    imink           i-index of each "qmink"
!    jmink           j-index of each "qmink"
!
!  grand maximum
!    qmax (qmin)     maximum (minimum) value of "qmaxk" ("qmink").
!    imax            i-index of "qmax"
!    jmax            j-index of "qmax"
!    kmax            k-index of "qmax"
!    imin            i-index of "qmin"
!    jmin            j-index of "qmin"
!    kmin            k-index of "qmin"
!
!----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use scratch
      use grid
!
      integer  :: j, k, i1, j1, k1, i2, j2, k2, imin, jmin, &
                  kmin, imax, jmax, kmax
      real(rl) :: qmin, qmax
!
      integer  :: imaxj(ijkn), iminj(ijkn), &
                  imaxk(ijkn), imink(ijkn), &
                  jmaxk(ijkn), jmink(ijkn)
!
      real(rl) :: qmaxj(ijkn), qminj(ijkn), &
                  qmaxk(ijkn), qmink(ijkn)
!
      real(rl) :: qty(in,jn,kn)
!
      integer  :: ijkl, ijkx, ijks, i
!
      real(rl) :: q(in*jn*kn)
!
! WARNING: Scratch array wg3d is used by this routine (mnmx)!
!
!
!      External statements
!
       integer       ISMIN   , ISMAX
       external      ISMIN   , ISMAX
!
!-----------------------------------------------------------------------
!
       k = k1
       j = j1
       ijkl = 0
       do 20 k=k1,k2
         do 10 j=j1,j2
           do 5 i=i1,i2
             ijkl = ijkl + 1
             q(ijkl) = qty(i,j,k)
5          continue
10       continue
20     continue
       ijkl = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)
       ijkx = ISMAX(ijkl,q,1)
       qmax = q(ijkx)
       kmax = ijkx/((i2-i1+1)*(j2-j1+1)) + k1
       jmax = (ijkx-(kmax-k1)*(i2-i1+1)*(j2-j1+1))/(i2-i1+1) + j1
       imax = ijkx-(kmax-k1)*(i2-i1+1)*(j2-j1+1) &
            - (jmax-j1)*(i2-i1+1) + i1
       ijks = ISMIN(ijkl,q,1)
       qmin = q(ijks)
       kmin = ijks/((i2-i1+1)*(j2-j1+1)) + k1
       jmin = (ijks-(kmin-k1)*(i2-i1+1)*(j2-j1+1))/(i2-i1+1) + j1
       imin = ijks-(kmin-k1)*(i2-i1+1)*(j2-j1+1) &
            - (jmin-j1)*(i2-i1+1) + i1
!
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                  M N M X                  \\\\\\\\\!
!=======================================================================
