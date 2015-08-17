!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine ggen
!
!  PURPOSE:  Initializes the grid in a new run according to the control
!  parameters in the input deck namelists "ggen1" and "ggen2".  All grid
!  variables are initialized.
!
!  LOCALS:
!
!  Modified 09/01/2006 by John Hayes; implemented Matthias Vigelius's
!  correction to DO loops 36, 136, and 236.  The original versions
!  resulted in incorrect computation of ratioed grids in multi-processor
!  calculations
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use domain
      use bndry
      use grid
      use root
      use scratch
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: nbl,igrid,imin,imax,jmin,jmax,kmin,kmax,iter,i,j,k
      real(rl) :: x1rat,dx1min,dfndx1r,x1r,deltx1r,errx1r, &
                  x2rat,dx2min,dfndx2r,x2r,deltx2r,errx2r, &
                  x3rat,dx3min,dfndx3r,x3r,deltx3r,errx3r, &
                  fn
!
      logical  :: lgrid
!
      integer  :: ibl, mbl, ifirst, ioffst, jfirst, joffst, &
                  kfirst, koffst
      integer  :: nzones(3), nzpt(3)
!
      integer, parameter :: iblmax = 10
!
      real(rl) :: gridd(5,iblmax)
      real(rl) :: remainder
!
      namelist /ggen1/ nbl,x1min,x1max,igrid,x1rat,dx1min,lgrid
      namelist /ggen2/ nbl,x2min,x2max,igrid,x2rat,dx2min,lgrid
      namelist /ggen3/ nbl,x3min,x3max,igrid,x3rat,dx3min,lgrid
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
!=======================================================================
!-----------  GENERATE X1 GRID  ----------------------------------------
!
!  Read in blocks of x1 grid zones.  Note we loop over read statement
!  until all blocks are read (signalled by reading in lgrid = .true.).
!  We can zone within each block completely independently of the others,
!  however we must ensure the starting position of one block (x1min) is
!  the same as the ending position (x1max) of the previous.
!   nbl    is number of active zones in block read in
!   x1min  is x1a(imin) ; bottom position of block
!   x1max  is x1a(imax) ; top    position of block
!   igrid  selects zoning type; we solve the zoning equation:
!            x1max = x1min + SUM OVER N[dx1min*x1rat**n]  ; so we must 
!          input either dx1min or x1rat (the other is calculated)
!      igrid   = 0  => block has already been set (irestart=1)
!              =+1  => (ratioed) use input "x1rat" to compute "dx1min",
!                      where x1a(imin+1) = x1a(imin) + dx1min.
!              =-1  => (ratioed) use input "x1rat" to compute "dx1min",
!                      where x1a(imax) = x1a(imax-1) + dx1min.
!              =+2  => (ratioed) use input "dx1min" to compute "x1rat",
!                      where x1a(imin+1) = x1a(imin) + dx1min.
!              =-2  => (ratioed) use input "dx1min" to compute "x1rat",
!                      where x1a(imax) = x1a(imax-1) + dx1min.
!   lgrid  logical flag for additional blocks ( =.true. reads another)
!   imax,imin  are indices of top and bottom of block
!
      is   = 3
      if (myid .eq. 0) then
        imax = is
        nbl   = 1
        x1min = 0.0
        x1max = 0.0
        igrid = 0
        x1rat = 0.0
        dx1min= 0.0
        lgrid = .false.
        ibl   = 0
        mbl   = 0
!
10      continue
        read (1,ggen1)
        write(2,ggen1)
!
        ibl  = ibl + 1
        if (ibl .gt. iblmax) then
          write(6,"(/1x,'ERROR: number of blocks in 1-direction exceeds' &
          ,' array bounds'/1x,'ibl = ',i4,'  iblmax = ',i4)") ibl,iblmax
          mbl  = 1
          go to 31
        endif
        imin = imax
        imax = imax + nbl
        if (imax .gt. is + (in-4)*ntiles(1)) then
          write(6,"(/1x,'ERROR: number of zones in 1-direction exceeds' &
          ,' array bounds',/1x,'imax = ',i4,'  in = ',i4)") imax,in
          mbl  = 1
          go to 31
        endif
!
!  1)  Compute dx1min from given value of x1rat.
!
        if (abs(igrid) .eq. 1) then
          if (igrid .eq. -1) x1rat = 1.0/x1rat
          if (x1rat .eq. 1.0)  then
            dx1min = (x1max-x1min)/real(nbl)
          else
            dx1min = (x1max-x1min)*(x1rat-1.0)/(x1rat**nbl - 1.0)
          endif
        endif
!
!  2)  Compute x1rat from given value of dx1min.  Newton Raphson 
!      iteration is required to find the root (x1rat) of the function:
!        fn(x1r) = (x1max-x1min) - dx1min*[(x1r)**nbl - 1]/[x1r-1] = 0
!
        if (abs(igrid) .eq. 2) then
          x1r = 1.01
          do 20 iter=1,20
            fn = (x1max - x1min) - dx1min*(x1r**nbl - 1.0)/(x1r - 1.0)
            dfndx1r  =  -nbl*dx1min*x1r**(nbl - 1)/(x1r - 1.0) &
                       + dx1min*(x1r**nbl - 1.0)/(x1r - 1.0)**2
            deltx1r  = -fn/dfndx1r
            errx1r   = abs(deltx1r/x1r)
            x1r = x1r + deltx1r
            if (errx1r .lt. 1.0e-6) goto 30
20        continue
          write(6,"(1x,'GRIDI: Newton-Raphson failed' &
          ,' for x1rat',/1x,'imin = ',i3,' x1r = ',1pe12.5,' deltx1r = ' &
          ,1e12.5,' fn = ',1e12.5)") imin,x1r,deltx1r,fn
          mbl  = 1
          go to 31
!
30        continue
          x1rat = x1r
        endif
        if (igrid .eq. -2) then
          dx1min = dx1min * x1rat**(nbl-1)
          x1rat  = 1.0 / x1rat
        endif
!
! Copy the grid descriptors into a temporary array for broadcast later.
!
        gridd(1,ibl) = float(imin)
        gridd(2,ibl) = float(imax)
        gridd(3,ibl) =  x1min
        gridd(4,ibl) = dx1min
        gridd(5,ibl) =  x1rat
!
!  Go back and read another block of x1 grid zones, if needed.
!
        if (.not. lgrid) go to 10
      endif
31    continue
!
! Broadcast ibl and the grid descriptor array to the others.
!
      if (myid .eq. 0) then
        ibuf_in(1) = mbl 
        ibuf_in(2) = ibl
      endif
      call MPI_BCAST(ibuf_in  ,2,    MPI_INTEGER,0,comm3d,ierr)
      if (myid .ne. 0) then
        mbl  = ibuf_in(1)
        ibl  = ibuf_in(2)
      endif
      if (mbl  .ne. 0) then
        call MPI_FINALIZE(ierr)
        stop
      endif
       call MPI_BCAST(gridd,5*iblmax,MPI_DOUBLE_PRECISION &
                    ,0,comm3d,ierr)
!
! The full grid has gridd(2,ibl) zones.  Find the best number
! of zones per tile.  Assume for now that the number of zones
! is evenly divisible by the number of tiles in that direction.
! Otherwise, load-balancing becomes an issue.
!
      nzones(1) = int(gridd(2,ibl)-gridd(1,1))
! Because of using multigrid Poisson Solver (mgmpi) for gravitation,
! number of zones is not longer restricted to be equal for
! all tiles and some of the following lines are commnted out.
!      nzpt(1)   = nzones(1)/ntiles(1)
!      fn        = real(nzones(1))/real(ntiles(1))
!      if (abs(real(nzpt(1))-fn) .gt. 0.001) then
!        if (myid .eq. 0) write(6,"(1x,'ERROR from GRIDI: nzones(1) = '
!     &  ,i4,' is not evenly divisible by ntiles(1) = ',i3)") nzones(1)
!     &  ,ntiles(1)
!#ifdef 1
!        call MPI_FINALIZE(ierr)
!#endif
!        stop
!      endif
! Relax the assumption of equal number of zones per processor, and
! go through some checks to determine the best zones distribution
! for the zone and tile parameters given.
!
      nzpt(1) = 0
      remainder = mod(nzones(1),ntiles(1))
      if(remainder .eq. 0) then
       nzpt(1) = nzones(1)/ntiles(1)
      else if((ntiles(1)-remainder) .eq. 1) then
       if(coords(1) .ne. ntiles(1)-1) nzpt(1) = nzones(1)/ntiles(1)+1
       if(coords(1) .eq. ntiles(1)-1) nzpt(1) = nzones(1)/ntiles(1)
      else if(((ntiles(1)-remainder) .eq. 3) .and. &
               (ntiles(1) .eq. 4)                 ) then
       if(coords(1) .le. 1) nzpt(1) = nzones(1)/ntiles(1)+1
       if(coords(1) .eq. 2) nzpt(1) = nzones(1)/ntiles(1)
       if(coords(1) .eq. 3) nzpt(1) = nzones(1)/ntiles(1)-1
      else
       if(coords(1) .ne. ntiles(1)-1) nzpt(1) = nzones(1)/ntiles(1)+1
       if(coords(1) .eq. ntiles(1)-1) nzpt(1) = nzones(1) - &
         (ntiles(1)-1)*(nzones(1)/ntiles(1) + 1)
      endif
      if(nzpt(1) .eq. 0) then
       print *,'GGEN1: could not set nzpt(1)'
       call MPI_FINALIZE(IERR)
       stop
      endif
!
! From my coordinates in the 3-D Cartesian virtual topology (coords),
! figure out which zones I own and compute the grid.  This tile's
! first real zone has index ifirst.   Locate the block containing it.
!
!      ioffst = coords(1) * nzpt(1)
      if(remainder .eq. 0) then
       ioffst = coords(1) * nzpt(1)
      else if((ntiles(1)-remainder) .eq. 1) then
       ioffst = coords(1) * (nzones(1)/ntiles(1)+1)
      else if(((ntiles(1)-remainder) .eq. 3) .and. &
               (ntiles(1) .eq. 4)                 ) then
       if(coords(1) .le. 2) &
          ioffst = coords(1)*(nzones(1)/ntiles(1)+1)
       if(coords(1) .eq. 3) &
          ioffst = 2*(nzones(1)/ntiles(1)+1) + nzones(1)/ntiles(1)
      else
       ioffst = coords(1) * (nzones(1)/ntiles(1)+1)
      endif
      ifirst = is + ioffst
      do 32 nbl=1,ibl
        mbl  = nbl
        if (int(gridd(1,mbl)).gt.ifirst) then
          mbl = max(1, mbl-1)
          go to 34
        endif
32    continue
34    continue
        imin = int(gridd(1,mbl))
        imax = int(gridd(2,mbl))
       x1min = gridd(3,mbl)
      dx1min = gridd(4,mbl)
       x1rat = gridd(5,mbl)
!
! The first real zone on this tile is in block mbl.  Compute the
! grid points of zones whose indices are less than ifirst, if any.
!
      if (imin .lt. ifirst) then
        do 36 i=imin+1,ifirst
!JH       dx1min = dx1min * x1rat
           x1min =  x1min + dx1min
          dx1min = dx1min * x1rat
36      continue
      endif
      imin = max(imin,ifirst)
      imax = min(imax,ifirst+nzpt(1))
!
! Now do the rest of the grid points in this block.
! Set up x1a grid lines from i=imin to imax, using known values of
! x1min, x1rat. 
!
       x1a(imin-ioffst) =  x1min
      dx1a(imin-ioffst) = dx1min
      do 40 i=imin+1,imax
        dx1a(i-ioffst) = dx1a(i-ioffst-1) * x1rat
         x1a(i-ioffst) =  x1a(i-ioffst-1) + dx1a(i-ioffst-1)
40    continue
      if (mbl.lt.ibl) then
        mbl = mbl + 1
        go to 34
      endif
!
!  Set up all grid zones and scale factors in x1 direction
!
      ie = imax - ioffst
!
!..................  GRIDI BOUNDARIES  .................................
!
! Get the coords of overlap ghost zones from the neighbors,
! unless I am on a physical boundary.
!
      nreq = 0
      if (coords(1).eq.0 .and. niis(1).ne.4) then
!
! I am on the physical inner boundary.  Compute ghost zones there.
!
! If the BCs are reflecting, the grid should be symmetric.
!
        dx1a(is-1) = dx1a(is  )
        if (niis(1).ne.2 .and. niis(1).ne.3) then
!
! Use this for Reflection-Symmetric grids: 
!
          dx1a(is-2) = dx1a(is+1)
          dx1min     = dx1a(is+2)
        else
!
! The following seems best for Flow In/Out:
!
          dx1a(is-2) = dx1a(is-1)
          dx1min     = dx1a(is-2)
        endif
      else
!
! Receive grid coords of the overlap zones from my "m" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(4),3,MPI_DOUBLE_PRECISION,n1m,n1m &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "m" neighbor.
!
        buf_in(1) = dx1a(is  )
        buf_in(2) = dx1a(is+1)
        buf_in(3) = dx1a(is+2)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(1),3,MPI_DOUBLE_PRECISION,n1m,myid &
                     ,comm3d,req(nreq),ierr)
      endif
!
      if (coords(1).eq.ntiles(1)-1 .and. nois(1).ne.4) then
!
! I am on the physical outer boundary.  Compute ghost zones there.
!
        dx1a(ie  ) = dx1a(ie-1)
        if (nois(1).ne.2 .and. nois(1).ne.3) then
!
! Use this for Reflection-Symmetric grids: 
!
          dx1a(ie+1) = dx1a(ie-2)
          dx1a(ie+2) = dx1a(ie-3)
        else
!
! Use this for Flow In/Out:
!
          dx1a(ie+1) = dx1a(ie  )
          dx1a(ie+2) = dx1a(ie+1)
        endif
      else
!
! Receive grid coords of the overlap zones from my "p" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(1),3,MPI_DOUBLE_PRECISION,n1p,n1p &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "p" neighbor.
!
        buf_in(4) = dx1a(ie-1)
        buf_in(5) = dx1a(ie-2)
        buf_in(6) = dx1a(ie-3)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(4),3,MPI_DOUBLE_PRECISION,n1p,myid &
                     ,comm3d,req(nreq),ierr)
      endif
      if (nreq .gt. 0) call MPI_WAITALL(nreq,req,stat,ierr)
      if (.not.(coords(1).eq.0 .and. niis(1).ne.4)) then
        dx1a(is-1) = buf_out(4)
        dx1a(is-2) = buf_out(5)
        dx1min     = buf_out(6)
      endif
      if (.not.(coords(1).eq.ntiles(1)-1 .and. nois(1).ne.4)) then
        dx1a(ie  ) = buf_out(1)
        dx1a(ie+1) = buf_out(2)
        dx1a(ie+2) = buf_out(3)
      endif
      x1a (is-1) = x1a (is  ) - dx1a(is-1)
      x1a (is-2) = x1a (is-1) - dx1a(is-2)
      x1min      = x1a (is-2) - dx1min
      x1a (ie+1) = x1a (ie  ) + dx1a(ie  )
      x1a (ie+2) = x1a (ie+1) + dx1a(ie+1)
      x1max      = x1a (ie+2) + dx1a(ie+2)
!
!................  GRIDI B-MESH AND GEOM FACTORS  ......................
!
       x1min     = x1min      + 0.5 * dx1min
       x1b(is-2) = x1a (is-2) + 0.5 * dx1a(is-2)
      dx1b(is-2) = x1b (is-2) - x1min
      do 50 i=is-1,ie+2
         x1b(i) = x1a(i) + 0.5*dx1a(i)
        dx1b(i) = x1b(i) - x1b(i-1)
50    continue
      do 60 i=is-2,ie+2
       if(lgeom .eq. 1 .or. lgeom .eq. 2) then
         g2a (i)  = 1.0
         g2b (i)  = 1.0
         g31a (i) = 1.0
         g31b (i) = 1.0
        dg2ad1(i) = 0.0
        dg2bd1(i) = 0.0
        dg31ad1(i) = 0.0
        dg31bd1(i) = 0.0
        vol1a  (i) = x1a(i)
        vol1b  (i) = x1b(i)
       endif ! CARTESIAN OR CYLINDRICAL
       if(lgeom .eq. 3) then
         g2a (i)  = x1a(i)
         g2b (i)  = x1b(i)
         g31a (i) = x1a(i)
         g31b (i) = x1b(i)
        dg2ad1(i) = 1.0
        dg2bd1(i) = 1.0
        dg31ad1(i) = 1.0
        dg31bd1(i) = 1.0
        vol1a  (i) = x1a(i)**3 / 3.0
        vol1b  (i) = x1b(i)**3 / 3.0
       endif ! SPHERICAL
        x1ai  (i) = 1.0 / ( x1a  (i) + tiny )
        x1bi  (i) = 1.0 / ( x1b  (i) + tiny )
        dx1ai (i) = 1.0 / ( dx1a (i) + tiny )
        dx1bi (i) = 1.0 / ( dx1b (i) + tiny )
        g2ai  (i) = 1.0 / ( g2a  (i) + tiny )
        g2bi  (i) = 1.0 / ( g2b  (i) + tiny )
        g31ai (i) = 1.0 / ( g31a (i) + tiny )
        g31bi (i) = 1.0 / ( g31b (i) + tiny )
60    continue
!
      if(lgeom .eq. 1 .or. lgeom .eq. 2) then
       dvl1b (is-2) = vol1b(is-2) - x1min
      endif ! CART OR CYL
!
      if(lgeom .eq. 3) then
       dvl1b (is-2) = vol1b(is-2) - x1min**3 / 3.0
      endif ! SPHERE
!
      dvl1bi(is-2) = 1.0 / ( dvl1b(is-2) + tiny )
      do 70 i=is-2,ie+1
        dvl1a (i  ) = vol1a(i+1) - vol1a(i)
        dvl1b (i+1) = vol1b(i+1) - vol1b(i)
        dvl1ai(i  ) = 1.0 / ( dvl1a(i  ) + tiny )
        dvl1bi(i+1) = 1.0 / ( dvl1b(i+1) + tiny )
70    continue
!
      if(lgeom .eq. 1 .or. lgeom .eq. 2) then
       dvl1a (ie+2) = x1max - vol1a(ie+2)
      endif ! CART OR CYL
!
      if(lgeom .eq. 3) then
       dvl1a (ie+2) = x1max**3 / 3.0 - vol1a(ie+2)
      endif ! SPHERE
      dvl1ai(ie+2) = 1.0 / ( dvl1a(ie+2) + tiny )
      if(xvgrid) then
!
!      Copy grid over to "n+1/2" and "n+1" grids.
!
       do 80 i=1,in
!
         x1ah   (i) = x1a   (i)
         x1bh   (i) = x1b   (i)
         dx1ah  (i) = dx1a  (i)
         dx1bh  (i) = dx1b  (i)
         g2ah   (i) = g2a   (i)
         g2bh   (i) = g2b   (i)
         g31ah  (i) = g31a  (i)
         g31bh  (i) = g31b  (i)
         dvl1ah (i) = dvl1a (i)
         dvl1bh (i) = dvl1b (i)
!
         x1ahi  (i) = x1ai  (i)
         x1bhi  (i) = x1bi  (i)
         dx1ahi (i) = dx1ai (i)
         dx1bhi (i) = dx1bi (i)
         g2ahi  (i) = g2ai  (i)
         g2bhi  (i) = g2bi  (i)
         g31ahi (i) = g31ai (i)
         g31bhi (i) = g31bi (i)
         dvl1ahi(i) = dvl1ai(i)
         dvl1bhi(i) = dvl1bi(i)
!
         x1an   (i) = x1a   (i)
         x1bn   (i) = x1b   (i)
         dx1an  (i) = dx1a  (i)
         dx1bn  (i) = dx1b  (i)
         g2an   (i) = g2a   (i)
         g2bn   (i) = g2b   (i)
         g31an  (i) = g31a  (i)
         g31bn  (i) = g31b  (i)
         dvl1an (i) = dvl1a (i)
         dvl1bn (i) = dvl1b (i)
!
         x1ani  (i) = x1ai  (i)
         x1bni  (i) = x1bi  (i)
         dx1ani (i) = dx1ai (i)
         dx1bni (i) = dx1bi (i)
         g2ani  (i) = g2ai  (i)
         g2bni  (i) = g2bi  (i)
         g31ani (i) = g31ai (i)
         g31bni (i) = g31bi (i)
         dvl1ani(i) = dvl1ai(i)
         dvl1bni(i) = dvl1bi(i)
!
80    continue
      endif ! xvgrid
!
!-----------  X2 GRID GENERATOR  ---------------------------------------
!
!
!  Variable names and values are the same as used in x1 grid generator
!
!
      js   = 3
      if (myid .eq. 0) then
        jmax = js
        nbl   = 1
        x2min = 0.0
        x2max = 0.0
        igrid = 0
        x2rat = 0.0
        dx2min= 0.0
        lgrid = .false.
        ibl   = 0
        mbl   = 0
!
110     continue
        read (1,ggen2)
        write(2,ggen2)
!
        ibl  = ibl + 1
        if (ibl .gt. iblmax) then
          write(6,"(/1x,'ERROR: number of blocks in 2-direction exceeds' &
          ,' array bounds'/1x,'ibl = ',i4,'  iblmax = ',i4)") ibl,iblmax
          mbl  = 1
          go to 131
        endif
        jmin = jmax
        jmax = jmax + nbl
        if (jmax .gt. js + (jn-4)*ntiles(2)) then
          write(6,"(/1x,'ERROR: number of zones in 2-direction exceeds' &
          ,' array bounds',/1x,'jmax = ',i4,'  jn = ',i4)") jmax,jn
          mbl  = 1
          go to 131
        endif
!
!  1)  Compute dx2min from given value of x2rat.
!
        if (abs(igrid) .eq. 1) then
          if (igrid .eq. -1) x2rat = 1.0/x2rat
          if (x2rat .eq. 1.0)  then
            dx2min = (x2max-x2min)/real(nbl)
          else
            dx2min = (x2max-x2min)*(x2rat-1.0)/(x2rat**nbl - 1.0)
          endif
        endif
!
!  2)  Compute x2rat from given value of dx2min.  Newton Raphson 
!      iteration is required to find the root (x2rat) of the function:
!        fn(x2r) = (x2max-x2min) - dx2min*[(x2r)**nbl - 1]/[x2r-1] = 0
!
        if (abs(igrid) .eq. 2) then
          x2r = 1.01
          do 120 iter=1,20
            fn = (x2max - x2min) - dx2min*(x2r**nbl - 1.0)/(x2r - 1.0)
            dfndx2r  =  -nbl*dx2min*x2r**(nbl - 1)/(x2r - 1.0) &
                       + dx2min*(x2r**nbl - 1.0)/(x2r - 1.0)**2
            deltx2r  = -fn/dfndx2r
            errx2r   = abs(deltx2r/x2r)
            x2r = x2r + deltx2r
            if (errx2r .lt. 1.0e-6) goto 130
120       continue
          write(6,"(1x,'GRIDJ: Newton-Raphson failed' &
          ,' for x2rat',/1x,'jmin = ',i3,' x2r = ',1pe12.5,' deltx2r = ' &
          ,1e12.5,' fn = ',1e12.5)") jmin,x2r,deltx2r,fn
          ierr = 1
          go to 131
!
130       continue
          x2rat = x2r
        endif
        if (igrid .eq. -2) then
          dx2min = dx2min * x2rat**(nbl-1)
          x2rat  = 1.0 / x2rat
        endif
!
! Copy the grid descriptors into a temporary array for broadcast later.
!
        gridd(1,ibl) = float(jmin)
        gridd(2,ibl) = float(jmax)
        gridd(3,ibl) =  x2min
        gridd(4,ibl) = dx2min
        gridd(5,ibl) =  x2rat
!
!  Go back and read another block of x2 grid zones, if needed.
!
        if (.not. lgrid) go to 110
      endif
131   continue
!
! Broadcast ibl and the grid descriptor array to the others.
!
      if (myid .eq. 0) then
        ibuf_in(1) = mbl
        ibuf_in(2) = ibl
      endif
      call MPI_BCAST(ibuf_in  ,2,    MPI_INTEGER,0,comm3d,ierr)
      if (myid .ne. 0) then
        mbl  = ibuf_in(1)
        ibl  = ibuf_in(2)
      endif
      if (mbl  .ne. 0) then
       call MPI_FINALIZE(ierr)
       stop
      endif
!
! Is it OK to send a 2-D array like this?
!
      call MPI_BCAST(gridd,5*iblmax,MPI_DOUBLE_PRECISION &
                    ,0,comm3d,ierr)
!
! The full grid has gridd(2,ibl) zones.  Find the best number
! of zones per tile.  Assume for now that the number of zones
! is evenly divisible by the number of tiles in that direction.
! Otherwise, load-balancing becomes an issue.
!
      nzones(2) = int(gridd(2,ibl)-gridd(1,1))
!      nzpt(2)   = nzones(2)/ntiles(2)
!      fn        = real(nzones(2))/real(ntiles(2))
!      if (abs(real(nzpt(2))-fn) .gt. 0.001) then
!        if (myid .eq. 0) write(6,"(1x,'ERROR from GRIDJ: nzones(2) = '
!     &  ,i4,' is not evenly divisible by ntiles(2) = ',i3)") nzones(2)
!     &  ,ntiles(2)
!#ifdef 1
!        call MPI_FINALIZE(ierr)
!#endif
!        stop
!      endif
!
! Relax the assumption of equal number of zones per processor, and
! go through some checks to determine the best zones distribution
! for the zone and tile parameters given.
!
      nzpt(2) = 0
      remainder = mod(nzones(2),ntiles(2))
      if(remainder .eq. 0) then
       nzpt(2) = nzones(2)/ntiles(2)
      else if((ntiles(2)-remainder) .eq. 1) then
       if(coords(2) .ne. ntiles(2)-1) nzpt(2) = nzones(2)/ntiles(2)+1
       if(coords(2) .eq. ntiles(2)-1) nzpt(2) = nzones(2)/ntiles(2)
      else if(((ntiles(2)-remainder) .eq. 3) .and. &
               (ntiles(2) .eq. 4)                 ) then
       if(coords(2) .le. 1) nzpt(2) = nzones(2)/ntiles(2)+1
       if(coords(2) .eq. 2) nzpt(2) = nzones(2)/ntiles(2)
       if(coords(2) .eq. 3) nzpt(2) = nzones(2)/ntiles(2)-1
      else
       if(coords(2) .ne. ntiles(2)-1) nzpt(2) = nzones(2)/ntiles(2)+1
       if(coords(2) .eq. ntiles(2)-1) nzpt(2) = nzones(2) - &
         (ntiles(2)-1)*(nzones(2)/ntiles(2) + 1)
      endif
      if(nzpt(2) .eq. 0) then
       print *,'GGEN2: could not set nzpt(2)'
       call MPI_FINALIZE(IERR)
       stop
      endif
!
! From my coordinates in the 3-D Cartesian virtual topology (coords),
! figure out which zones I own and compute the grid.  This tile's
! first real zone has index ifirst.   Locate the block containing it.
!
!      joffst = coords(2) * nzpt(2)
      if(remainder .eq. 0) then
       joffst = coords(2) * nzpt(2)
      else if((ntiles(2)-remainder) .eq. 1) then
       joffst = coords(2) * (nzones(2)/ntiles(2)+1)
      else if(((ntiles(2)-remainder) .eq. 3) .and. &
               (ntiles(2) .eq. 4)                 ) then
       if(coords(2) .le. 2) &
          joffst = coords(2)*(nzones(2)/ntiles(2)+1)
       if(coords(2) .eq. 3) &
          joffst = 2*(nzones(2)/ntiles(2)+1) + nzones(2)/ntiles(2)
      else
       joffst = coords(2) * (nzones(2)/ntiles(2)+1)
      endif
      jfirst = js + joffst
      do 132 nbl=1,ibl
        mbl  = nbl
        if (int(gridd(1,mbl)).gt.jfirst) then
          mbl = max(1, mbl-1)
          go to 134
        endif
132   continue
134   continue
        jmin = int(gridd(1,mbl))
        jmax = int(gridd(2,mbl))
       x2min = gridd(3,mbl)
      dx2min = gridd(4,mbl)
       x2rat = gridd(5,mbl)
!
! The first real zone on this tile is in block mbl.  Compute the
! grid points of zones whose indeces are less than jfirst, if any.
!
      if (jmin .lt. jfirst) then
        do 136 j=jmin+1,jfirst
!JH       dx2min = dx2min * x2rat
           x2min =  x2min + dx2min
          dx2min = dx2min * x2rat
136     continue
      endif
      jmin = max(jmin,jfirst)
      jmax = min(jmax,jfirst+nzpt(2))
!
! Now do the rest of the grid points in this block.
! Set up x2a grid lines from j=jmin to jmax, using known values of
! x2min, x2rat. 
!
       x2a(jmin-joffst) =  x2min
      dx2a(jmin-joffst) = dx2min
      do 140 j=jmin+1,jmax
        dx2a(j-joffst) = dx2a(j-joffst-1) * x2rat
         x2a(j-joffst) =  x2a(j-joffst-1) + dx2a(j-joffst-1)
140   continue
      if (mbl.lt.ibl) then
        mbl = mbl + 1
        go to 134
      endif
!
!  Set up all grid zones, scale factors in x2 direction
!
      je = jmax - joffst
!
!
!..................  GRIDJ BOUNDARIES  .................................
!
! Get the coords of overlap ghost zones from the neighbors,
! unless I am on a physical boundary.
!
      nreq = 0
      if (coords(2).eq.0 .and. nijs(1).ne.4) then
!
! I am on the physical inner boundary.  Compute ghost zones there.
!
! If the BCs are reflecting, the grid should be symmetric.
!
        dx2a(js-1) = dx2a(js  )
        if (nijs(1).ne.2 .and. nijs(1).ne.3) then
!
! Use this for Reflection-Symmetric grids: 
!
          dx2a(js-2) = dx2a(js+1)
          dx2min     = dx2a(js+1)
        else
!
! The following seems best for Flow In/Out:
!
          dx2a(js-2) = dx2a(js-1)
          dx2min     = dx2a(js-2)
        endif
!
      else
!
! Receive grid coords of the overlap zones from my "m" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(4),3,MPI_DOUBLE_PRECISION,n2m,n2m &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "m" neighbor.
!
        buf_in(1) = dx2a(js  )
        buf_in(2) = dx2a(js+1)
        buf_in(3) = dx2a(js+2)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(1),3,MPI_DOUBLE_PRECISION,n2m,myid &
                     ,comm3d,req(nreq),ierr)
      endif
!
      if (coords(2).eq.ntiles(2)-1 .and. nojs(1).ne.4) then
!
! I am on the physical outer boundary.  Compute ghost zones there.
!
        dx2a(je  ) = dx2a(je-1)
        if (nojs(1).ne.2 .and. nojs(1).ne.3) then
!
! Use this for Reflection-Symmetric and Periodic grids: 
!
          dx2a(je+1) = dx2a(je-2)
          dx2min     = dx2a(je-3)
        else
!
! Use this for Flow In/Out:
!
          dx2a(je+1) = dx2a(je  )
          dx2a(je+2) = dx2a(je+1)
        endif
      else
!
! Receive grid coords of the overlap zones from my "p" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(1),3,MPI_DOUBLE_PRECISION,n2p,n2p &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "p" neighbor.
!
        buf_in(4) = dx2a(je-1)
        buf_in(5) = dx2a(je-2)
        buf_in(6) = dx2a(je-3)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(4),3,MPI_DOUBLE_PRECISION,n2p,myid &
                     ,comm3d,req(nreq),ierr)
      endif
      if (nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
      if (.not.(coords(2).eq.0 .and. nijs(1).ne.4)) then
        dx2a(js-1) = buf_out(4)
        dx2a(js-2) = buf_out(5)
        dx2min     = buf_out(6)
      endif
      if (.not.(coords(2).eq.ntiles(2)-1 .and. nojs(1).ne.4)) then
        dx2a(je  ) = buf_out(1)
        dx2a(je+1) = buf_out(2)
        dx2a(je+2) = buf_out(3)
      endif
      x2a (js-1) = x2a (js  ) - dx2a(js-1)
      x2a (js-2) = x2a (js-1) - dx2a(js-2)
      x2min      = x2a (js-2) - dx2min
      x2a (je+1) = x2a (je  ) + dx2a(je  )
      x2a (je+2) = x2a (je+1) + dx2a(je+1)
      x2max      = x2a (je+2) + dx2a(je+2)
!
!................  GRIDJ B-MESH AND GEOM FACTORS  ......................
!
       x2min     = x2min      + 0.5 * dx2min 
       x2b(js-2) = x2a (js-2) + 0.5 * dx2a(js-2)
      dx2b(js-2) = x2b (js-2) - x2min
      do 150 j=js-1,je+2
         x2b(j) = x2a(j) + 0.5 * dx2a(j)
        dx2b(j) = x2b(j) - x2b(j-1)
150   continue
!
      do 160 j=js-2,je+2
       if(lgeom .eq. 1) then
         g32a (j)  = 1.0
         g32b (j)  = 1.0
        dg32ad2(j) = 0.0
        dg32bd2(j) = 0.0
        vol2a  (j) = x2a(j)
        vol2b  (j) = x2b(j)
       endif ! CART
       if(lgeom .eq. 2) then
         g32a (j)  = x2a(j)
         g32b (j)  = x2b(j)
        dg32ad2(j) = 1.0
        dg32bd2(j) = 1.0
        vol2a  (j) = x2a(j)**2 / 2.0
        vol2b  (j) = x2b(j)**2 / 2.0
       endif ! CYL
       if(lgeom .eq. 3) then
         g32a (j)  = sin( x2a(j) )
         g32b (j)  = sin( x2b(j) )
        dg32ad2(j) = cos( x2a(j) )
        dg32bd2(j) = cos( x2b(j) )
        vol2a (j)  =-cos( x2a(j) )
        vol2b (j)  =-cos( x2b(j) )
       endif ! SPHERE
       x2ai  (j) = 1.0 / ( x2a  (j) + tiny )
       x2bi  (j) = 1.0 / ( x2b  (j) + tiny )
       dx2ai (j) = 1.0 / ( dx2a (j) + tiny )
       dx2bi (j) = 1.0 / ( dx2b (j) + tiny )
!JH -- edit to work around coordinate singularity...
       if(lgeom .eq. 1) then
        g32ai (j) = 1.0 / ( g32a (j) + tiny )
       else
        if(g32a(j) .eq. 0.0) then
         g32ai(j) = 0.0d0
        else
         g32ai (j) = 1.0 / ( g32a (j) + tiny )
        endif
       endif
       g32bi (j) = 1.0 / ( g32b (j) + tiny )
160   continue
!
      if(lgeom .eq. 1 .or. lgeom .eq. 2) then
       dvl2b (js-2) = vol2b(js-2) - x2min
      endif ! CART or CYL
!
      if(lgeom .eq. 3) then
       dvl2b (js-2) = vol2b(js-2) - x2min**2 / 2.0
      endif ! SPHERE
!
      dvl2bi(js-2) = 1.0 / ( dvl2b(js-2) + tiny )
      do 170 j=js-2,je+1
       dvl2a (j  ) = vol2a(j+1) - vol2a(j)
       dvl2b (j+1) = vol2b(j+1) - vol2b(j)
       if(dvl2b(j+1) .eq. 0.0) dvl2b(j+1) = dvl2a(j)
       dvl2ai(j  ) = 1.0 / ( dvl2a(j  ) + tiny )
       dvl2bi(j+1) = 1.0 / ( dvl2b(j+1) + tiny )
170   continue
      if(lgeom .eq. 1 .or. lgeom .eq. 2) then
       dvl2a (je+2) = x2max - vol2a(je+2)
      endif ! CART or CYL
!
      if(lgeom .eq. 3) then
       dvl2a (je+2) = x2max**2 / 2.0 - vol2a(je+2)
      endif ! SPHERE
      dvl2ai(je+2) = 1.0 / ( dvl2a(je+2) + tiny )
!
      if(xvgrid) then
!
!      Copy grid over to "n+1/2" and "n+1" grids.
!
       do 180 j=1,jn
!
         x2ah   (j) = x2a   (j)
         x2bh   (j) = x2b   (j)
         dx2ah  (j) = dx2a  (j)
         dx2bh  (j) = dx2b  (j)
         g32ah  (j) = g32a  (j)
         g32bh  (j) = g32b  (j)
         dvl2ah (j) = dvl2a (j)
         dvl2bh (j) = dvl2b (j)
!
         x2ahi  (j) = x2ai  (j)
         x2bhi  (j) = x2bi  (j)
         dx2ahi (j) = dx2ai (j)
         dx2bhi (j) = dx2bi (j)
         g32ahi (j) = g32ai (j)
         g32bhi (j) = g32bi (j)
         dvl2ahi(j) = dvl2ai(j)
         dvl2bhi(j) = dvl2bi(j)
!
         x2an   (j) = x2a   (j)
         x2bn   (j) = x2b   (j)
         dx2an  (j) = dx2a  (j)
         dx2bn  (j) = dx2b  (j)
         g32an  (j) = g32a  (j)
         g32bn  (j) = g32b  (j)
         dvl2an (j) = dvl2a (j)
         dvl2bn (j) = dvl2b (j)
!
         x2ani  (j) = x2ai  (j)
         x2bni  (j) = x2bi  (j)
         dx2ani (j) = dx2ai (j)
         dx2bni (j) = dx2bi (j)
         g32ani (j) = g32ai (j)
         g32bni (j) = g32bi (j)
         dvl2ani(j) = dvl2ai(j)
         dvl2bni(j) = dvl2bi(j)
!
180    continue
      endif ! xvgrid
!
!-----------  X3 GRID GENERATOR  ---------------------------------------
!
!
!  Variable names and values are the same as used in x1 grid generator
!
      ks   = 3
      if (myid .eq. 0) then
        kmax = ks
        nbl   = 1
        x3min = 0.0
        x3max = 0.0
        igrid = 0
        x3rat = 0.0
        dx3min= 0.0
        lgrid = .false.
        ibl   = 0
        mbl   = 0
!
210     continue
        read (1,ggen3)
        write(2,ggen3)
!
        ibl  = ibl + 1
        if (ibl .gt. iblmax) then
          write(6,"(/1x,'ERROR: number of blocks in 3-direction exceeds' &
          ,' array bounds'/1x,'ibl = ',i4,'  iblmax = ',i4)") ibl,iblmax
          mbl  = 1
          go to 231
        endif
        kmin = kmax
        kmax = kmax + nbl
        if (kmax .gt. ks + (kn-4)*ntiles(3)) then
          write(6,"(/1x,'ERROR: number of zones in 3-direction exceeds' &
          ,' array bounds',/1x,'kmax = ',i4,'  kn = ',i4)") kmax,kn
          mbl  = 1
          go to 231
        endif
!
!  1)  Compute dx3min from given value of x3rat.
!
        if (abs(igrid) .eq. 1) then
          if (igrid .eq. -1) x3rat = 1.0/x3rat
          if (x3rat .eq. 1.0)  then
            dx3min = (x3max-x3min)/real(nbl)
          else
            dx3min = (x3max-x3min)*(x3rat-1.0)/(x3rat**nbl - 1.0)
          endif
        endif
!
!  2)  Compute x3rat from given value of dx3min.  Newton Raphson 
!      iteration is required to find the root (x3rat) of the function:
!        fn(x3r) = (x3max-x3min) - dx3min*[(x3r)**nbl - 1]/[x3r-1] = 0
!
        if (abs(igrid) .eq. 2) then
          x3r = 1.01
          do 220 iter=1,20
            fn = (x3max - x3min) - dx3min*(x3r**nbl - 1.0)/(x3r - 1.0)
            dfndx3r  =  -nbl*dx3min*x3r**(nbl - 1)/(x3r - 1.0) &
                       + dx3min*(x3r**nbl - 1.0)/(x3r - 1.0)**2
            deltx3r  = -fn/dfndx3r
            errx3r   = abs(deltx3r/x3r)
            x3r = x3r + deltx3r
            if (errx3r .lt. 1.0e-6) goto 230
220       continue
          write(6,"(1x,'GRIDK: Newton-Raphson failed' &
          ,' for x3rat',/1x,'kmin = ',i3,' x3r = ',1pe12.5,' deltx3r = ' &
          ,1e12.5,' fn = ',1e12.5)") kmin,x3r,deltx3r,fn
          mbl  = 1
          go to 231
!
230       continue
          x3rat = x3r
        endif
        if (igrid .eq. -2) then
          dx3min = dx3min * x3rat**(nbl-1)
          x3rat  = 1.0 / x3rat
        endif
!
! Copy the grid descriptors into a temporary array for broadcast later.
!
        gridd(1,ibl) = float(kmin)
        gridd(2,ibl) = float(kmax)
        gridd(3,ibl) =  x3min
        gridd(4,ibl) = dx3min
        gridd(5,ibl) =  x3rat
!
!  Go back and read another block of x3 grid zones, if needed.
!
        if (.not. lgrid) go to 210
      endif
231   continue
!
! Broadcast ibl and the grid descriptor array to the others.
!
      if (myid .eq. 0) then
        ibuf_in(1) = mbl 
        ibuf_in(2) = ibl
      endif
      call MPI_BCAST(ibuf_in  ,2,    MPI_INTEGER,0,comm3d,ierr)
      if (myid .ne. 0) then
        mbl  = ibuf_in(1)
        ibl  = ibuf_in(2)
      endif
      if (mbl  .ne. 0) then
        call MPI_FINALIZE(ierr)
        stop
      endif
       call MPI_BCAST(gridd,5*iblmax,MPI_DOUBLE_PRECISION &
                    ,0,comm3d,ierr)
!
! The full grid has gridd(2,ibl) zones.  Find the best number
! of zones per tile.  Assume for now that the number of zones
! is evenly divisible by the number of tiles in that direction.
! Otherwise, load-balancing becomes an issue.
!
      nzones(3) = int(gridd(2,ibl)-gridd(1,1))
!      nzpt(3)   = nzones(3)/ntiles(3)
!      fn        = real(nzones(3))/real(ntiles(3))
!      if (abs(real(nzpt(3))-fn) .gt. 0.001) then
!        if (myid .eq. 0) write(6,"(1x,'ERROR from GRIDK: nzones(3) = '
!     &  ,i4,' is not evenly divisible by ntiles(3) = ',i3)") nzones(3)
!     &  ,ntiles(3)
!#ifdef 1
!        call MPI_FINALIZE(ierr)
!#endif
!        stop
!      endif
!
! Relax the assumption of equal number of zones per processor, and
! go through some checks to determine the best zones distribution
! for the zone and tile parameters given.
!
      nzpt(3) = 0
      remainder = mod(nzones(3),ntiles(3))
      if(remainder .eq. 0) then
       nzpt(3) = nzones(3)/ntiles(3)
      else if((ntiles(3)-remainder) .eq. 1) then
       if(coords(3) .ne. ntiles(3)-1) nzpt(3) = nzones(3)/ntiles(3)+1
       if(coords(3) .eq. ntiles(3)-1) nzpt(3) = nzones(3)/ntiles(3)
      else if(((ntiles(3)-remainder) .eq. 3) .and. &
               (ntiles(3) .eq. 4)                 ) then
       if(coords(3) .le. 1) nzpt(3) = nzones(3)/ntiles(3)+1
       if(coords(3) .eq. 2) nzpt(3) = nzones(3)/ntiles(3)
       if(coords(3) .eq. 3) nzpt(3) = nzones(3)/ntiles(3)-1
      else
       if(coords(3) .ne. ntiles(3)-1) nzpt(3) = nzones(3)/ntiles(3)+1
       if(coords(3) .eq. ntiles(3)-1) nzpt(3) = nzones(3) - &
         (ntiles(3)-1)*(nzones(3)/ntiles(3) + 1)
      endif
      if(nzpt(3) .eq. 0) then
       print *,'GGEN3: could not set nzpt(3)'
       call MPI_FINALIZE(IERR)
       stop
      endif
!
! From my coordinates in the 3-D Cartesian virtual topology (coords),
! figure out which zones I own and compute the grid.  This tile's
! first real zone has index ifirst.   Locate the block containing it.
!
!      koffst = coords(3) * nzpt(3)
      if(remainder .eq. 0) then
       koffst = coords(3) * nzpt(3)
      else if((ntiles(3)-remainder) .eq. 1) then
       koffst = coords(3) * (nzones(3)/ntiles(3)+1)
      else if(((ntiles(3)-remainder) .eq. 3) .and. &
               (ntiles(3) .eq. 4)                 ) then
       if(coords(3) .le. 2) &
          koffst = coords(3)*(nzones(3)/ntiles(3)+1)
       if(coords(3) .eq. 3) &
          koffst = 2*(nzones(3)/ntiles(3)+1) + nzones(3)/ntiles(3)
      else
       koffst = coords(3) * (nzones(3)/ntiles(3)+1)
      endif
      kfirst = ks + koffst
      do 232 nbl=1,ibl
        mbl  = nbl
        if (int(gridd(1,mbl)).gt.kfirst) then
          mbl = max(1, mbl-1)
          go to 234
        endif
232   continue
234   continue
        kmin = int(gridd(1,mbl))
        kmax = int(gridd(2,mbl))
       x3min = gridd(3,mbl)
      dx3min = gridd(4,mbl)
       x3rat = gridd(5,mbl)
!
! The first real zone on this tile is in block mbl.  Compute the
! grid points of zones whose indeces are less than kfirst, if any.
!
      if (kmin .lt. kfirst) then
        do 236 k=kmin+1,kfirst
!JH       dx3min = dx3min * x3rat
           x3min =  x3min + dx3min
          dx3min = dx3min * x3rat
236     continue
      endif
      kmin = max(kmin,kfirst)
      kmax = min(kmax,kfirst+nzpt(3))
!
! Now do the rest of the grid points in this block.
! Set up x3a grid lines from k=kmin to kmax, using known values of
! x3min, x3rat. 
!
       x3a(kmin-koffst) =  x3min
      dx3a(kmin-koffst) = dx3min
      do 240 k=kmin+1,kmax
        dx3a(k-koffst) = dx3a(k-koffst-1) * x3rat
         x3a(k-koffst) =  x3a(k-koffst-1) + dx3a(k-koffst-1)
240   continue
      if (mbl.lt.ibl) then
        mbl = mbl + 1
        go to 234
      endif
!
!  Set up all grid zones, scale factors in x3 direction
!
      ke = kmax - koffst
!
!
!..................  GRIDK BOUNDARIES  .................................
!
! Get the coords of overlap ghost zones from the neighbors,
! unless I am on a physical boundary.
!
      nreq = 0
      if (coords(3).eq.0 .and. niks(1).ne.4) then
!
! I am on the physical inner boundary.  Compute ghost zones there.
!
! If the BCs are reflecting, the grid should be symmetric.
!
        dx3a(ks-1) = dx3a(ks  )
        if (niks(1).ne.2 .and. niks(1).ne.3) then
!
! Use this for Reflection-Symmetric grids: 
!
          dx3a(ks-2) = dx3a(ks+1)
          dx3min     = dx3a(ks+2)
        else
!
! The following seems best for Flow In/Out:
!
          dx3a(ks-2) = dx3a(ks-1)
          dx3min     = dx3a(ks-2)
        endif
      else
!
! Receive grid coords of the overlap zones from my "m" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(4),3,MPI_DOUBLE_PRECISION,n3m,n3m &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "m" neighbor.
!
        buf_in(1) = dx3a(ks  )
        buf_in(2) = dx3a(ks+1)
        buf_in(3) = dx3a(ks+2)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(1),3,MPI_DOUBLE_PRECISION,n3m,myid &
                     ,comm3d,req(nreq),ierr)
      endif
!
      if (coords(3).eq.ntiles(3)-1 .and. noks(1).ne.4) then
!
! I am on the physical outer boundary.  Compute ghost zones there.
!
        dx3a(ke+1) = dx3a(ke  )
        if (noks(1).ne.2 .and. noks(1).ne.3) then
!
! Use this for Reflection-Symmetric and Periodic grids: 
!
          dx3a(ke+1) = dx3a(ke-2)
          dx3a(ke+2) = dx3a(ke-3)
        else
!
! Use this for Flow In/Out:
!
          dx3a(ke+1) = dx3a(ke  )
          dx3a(ke+2) = dx3a(ke+1)
        endif
      else
!
! Receive grid coords of the overlap zones from my "p" neighbor.
!
        nreq = nreq + 1
        call MPI_IRECV(buf_out(1),3,MPI_DOUBLE_PRECISION,n3p,n3p &
                     ,comm3d,req(nreq),ierr)
!
! Send grid coords of the overlap zones to my "p" neighbor.
!
        buf_in(4) = dx3a(ke-1)
        buf_in(5) = dx3a(ke-2)
        buf_in(6) = dx3a(ke-3)
        nreq = nreq + 1
        call MPI_ISEND(buf_in(4),3,MPI_DOUBLE_PRECISION,n3p,myid &
                     ,comm3d,req(nreq),ierr)
      endif
      if (nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
      if (.not.(coords(3).eq.0 .and. niks(1).ne.4)) then
        dx3a(ks-1) = buf_out(4)
        dx3a(ks-2) = buf_out(5)
        dx3min     = buf_out(6)
      endif
      if (.not.(coords(3).eq.ntiles(3)-1 .and. noks(1).ne.4)) then
        dx3a(ke  ) = buf_out(1)
        dx3a(ke+1) = buf_out(2)
        dx3a(ke+2) = buf_out(3)
      endif
      x3a (ks-1) = x3a (ks  ) - dx3a(ks-1)
      x3a (ks-2) = x3a (ks-1) - dx3a(ks-2)
      x3min      = x3a (ks-2) - dx3min
      x3a (ke+1) = x3a (ke  ) + dx3a(ke  )
      x3a (ke+2) = x3a (ke+1) + dx3a(ke+1)
      x3max      = x3a (ke+2) + dx3a(ke+2)
!
!................  GRIDK B-MESH AND GEOM FACTORS  ......................
!
       x3min     = x3min      + 0.5 * dx3min
       x3b(ks-2) = x3a (ks-2) + 0.5 * dx3a(ks-2)
      dx3b(ks-2) = x3b (ks-2) - x3min
      do 250 k=ks-1,ke+2
         x3b(k) = x3a(k) + 0.5 * dx3a(k)
        dx3b(k) = x3b(k) - x3b(k-1)
250   continue
      do 260 k=ks-2,ke+2
        dx3ai(k)  = 1.0 / ( dx3a(k) + tiny )
        dx3bi(k)  = 1.0 / ( dx3b(k) + tiny )
       vol3a  (k) = x3a(k)
       vol3b  (k) = x3b(k)
260   continue
!
      dvl3b (ks-2) = vol3b(ks-2) - x3min
      dvl3bi(ks-2) = 1.0 / ( dvl3b(ks-2) + tiny )
      do 270 k=ks-2,ke+1
        dvl3a (k  ) = vol3a(k+1) - vol3a(k)
        dvl3b (k+1) = vol3b(k+1) - vol3b(k)
        dvl3ai(k  ) = 1.0 / ( dvl3a(k  ) + tiny )
        dvl3bi(k+1) = 1.0 / ( dvl3b(k+1) + tiny )
270   continue
      dvl3a (ke+2) = x3max - vol3a(ke+2)
      dvl3ai(ke+2) = 1.0 / ( dvl3a(ke+2) + tiny )
      if(xvgrid) then
!
!      Copy grid over to "n+1/2" and "n+1" grids.
!
       do 280 k=1,kn
!
         x3ah   (k) = x3a   (k)
         x3bh   (k) = x3b   (k)
         dx3ah  (k) = dx3a  (k)
         dx3bh  (k) = dx3b  (k)
         dvl3ah (k) = dvl3a (k)
         dvl3bh (k) = dvl3b (k)
!
         x3ahi  (k) = x3ai  (k)
         x3bhi  (k) = x3bi  (k)
         dx3ahi (k) = dx3ai (k)
         dx3bhi (k) = dx3bi (k)
         dvl3ahi(k) = dvl3ai(k)
         dvl3bhi(k) = dvl3bi(k)
!
         x3an   (k) = x3a   (k)
         x3bn   (k) = x3b   (k)
         dx3an  (k) = dx3a  (k)
         dx3bn  (k) = dx3b  (k)
         dvl3an (k) = dvl3a (k)
         dvl3bn (k) = dvl3b (k)
!
         x3ani  (k) = x3ai  (k)
         x3bni  (k) = x3bi  (k)
         dx3ani (k) = dx3ai (k)
         dx3bni (k) = dx3bi (k)
         dvl3ani(k) = dvl3ai(k)
         dvl3bni(k) = dvl3bi(k)
!
280    continue
      endif ! xvgrid
!
!-----------------------------------------------------------------------
!
!  is,ie [js,je] are starting and ending indices of ACTIVE i [j] zones,
!  so modify accordingly
!
      ie = ie - 1
      je = je - 1
      ke = ke - 1
      nx1z = ie - is + 1
      nx2z = je - js + 1
      nx3z = ke - ks + 1
!
      return
      end
