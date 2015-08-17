!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 N E W X 1                 \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
      subroutine newx1
!
!  PURPOSE: Computes "new" x1 grid variables (grid variables at advanced
!  timestep) to be used in TRANSPRT.  Grid values are calculated for
!  i=is-2 to ie+2, except for dvl1a (i=is,ie+2) and dvl1b (i=is+1,ie+2).
!
!     implemented in ZEUS-MP by John Hayes
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use grid
      use root
      use scratch
      use bndry
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, ibeg
      real(rl) :: vol1an(in) , vol1bn(in), qa, qb, qc, qd, vfac
      real(rl) :: ibuf(2)
!
!=======================================================================
!
      if(lgeom .eq. 1) vfac = 1.0D0
      if(lgeom .eq. 2) vfac = 1.0D0
      if(lgeom .eq. 3) vfac = 1.0D0/3.0D0
!
!-----------------------------------------------------------------------
!     post receive requests for updated grid values near tile boundaries
!-----------------------------------------------------------------------
!
      nreq = 0
      if(nois(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_irecv(ibuf(1), 1 , MPI_DOUBLE_PRECISION, &
                      n1p, 666, comm3d, req(nreq), ierr)
      endif
!
!-----------------------------------------------------------------------
!     update x1an straddling inner boundary...
!-----------------------------------------------------------------------
!
      x1an(is-2) = x1a(is-2) + vg1(is-2)*dt
      do 10 i=is-1,is+2
         x1an(i  ) = x1a (i) + vg1(i)*dt
        dx1an(i-1) = x1an(i) - x1an(i-1)
10    continue
!
!-----------------------------------------------------------------------
!     ...and if this boundary is a tile boundary, send the updated
!     dx1an to the neighboring processor.
!-----------------------------------------------------------------------
!
      if(niis(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_isend(dx1an(is+1), 1, MPI_DOUBLE_PRECISION, &
                      n1m, 666, comm3d, req(nreq), ierr)
      endif
!
!-----------------------------------------------------------------------
!     finish updating x1an, dx1an away from inner tile boundary...
!-----------------------------------------------------------------------
!
      do 11 i=is+3,ie+2
         x1an(i  ) = x1a (i) + vg1(i)*dt
        dx1an(i-1) = x1an(i) - x1an(i-1)
11    continue
!
!-----------------------------------------------------------------------
!     ...and give all tiles a chance to collect the boundary values
!     and finish up.
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
      if(nois(1) .ne. 0) then
       dx1an(ie+2) = (dx1an(ie+1)/dx1an(ie)) * dx1an(ie+1)
      else
       dx1an(ie+2) = ibuf(1)
      endif
!
!-----------------------------------------------------------------------
!     repeat process for "b" grid
!-----------------------------------------------------------------------
!
      if(niis(1) .eq. 0) then
       ibeg = is-2 
      else
       ibeg = is-1
      endif
      if(niis(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_irecv(ibuf(1), 2, MPI_DOUBLE_PRECISION, &
                      n1m, 667, comm3d, req(nreq), ierr)
      endif
      do 20 i=ie-2, ie+2
         x1bn(i) = x1an(i) + 0.5*dx1an(i)
20    continue
      if(nois(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_isend(x1bn(ie-2), 2, MPI_DOUBLE_PRECISION, &
                      n1p, 667, comm3d, req(nreq), ierr)
      endif
      do 22 i=ibeg, ie-3
         x1bn(i) = x1an(i) + 0.5*dx1an(i)
22    continue
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
!
      if(niis(1) .eq. 0) then
       x1bn(is-2) = ibuf(2)
       do i=is-1, ie+2
         dx1bn(i) = x1bn(i) - x1bn(i-1)
       enddo
       dx1bn(is-2) = x1bn(is-2) - ibuf(1)
      else
       dx1bn(is-2) = dx1an(is-2)
        x1bn(is-2) =  x1an(is-1) - 0.5*dx1an(is-2)
       do i=is-1, ie+2
         dx1bn(i) = x1bn(i) - x1bn(i-1)
       enddo
      endif
!
!-----------------------------------------------------------------------
!     compute inverse coordinate values at advanced time...
!-----------------------------------------------------------------------
!
      do i = is-2,ie+2
       x1ani(i) = 1.0/(x1an(i)+tiny)
       x1bni(i) = 1.0/(x1bn(i)+tiny)
      enddo
!
!-----------------------------------------------------------------------
!     ...and update the time-centered metric arrays AND their inverses
!-----------------------------------------------------------------------
!
      do 30 i=is-2,ie+2
       if(lgeom .eq. 1) then
        g2ah  (i) = 1.0
        g2bh  (i) = 1.0
        g31ah (i) = 1.0
        g31bh (i) = 1.0
        g2an  (i) = 1.0
        g2bn  (i) = 1.0
        g31an (i) = 1.0
        g31bn (i) = 1.0
       endif ! XYZ
       if(lgeom .eq. 2) then
        g2ah  (i) = 1.0
        g2bh  (i) = 1.0
        g31ah (i) = 1.0
        g31bh (i) = 1.0
        g2an  (i) = 1.0
        g2bn  (i) = 1.0
        g31an (i) = 1.0
        g31bn (i) = 1.0
       endif ! ZRP
       if(lgeom .eq. 3) then
        g2ah  (i) = 0.5*(x1a(i) + x1an(i))
        g2bh  (i) = 0.5*(x1b(i) + x1bn(i))
        g31ah (i) = 0.5*(x1a(i) + x1an(i))
        g31bh (i) = 0.5*(x1b(i) + x1bn(i))
        g2an  (i) = x1an (i)
        g2bn  (i) = x1bn (i)
        g31an (i) = x1an (i)
        g31bn (i) = x1bn (i)
       endif ! RTP
30    continue
      do i = is-2, ie+2
       g2ani (i) = 1.0/(g2an (i)+tiny)
       g2bni (i) = 1.0/(g2bn (i)+tiny)
       g31ani(i) = 1.0/(g31an(i)+tiny)
       g31bni(i) = 1.0/(g31bn(i)+tiny)
      enddo
!
!-----------------------------------------------------------------------
!     New volume factors
!-----------------------------------------------------------------------
!
      vol1an(is-2) = g2an(is-2)*g31an(is-2)*x1an(is-2)
      do 40 i=is-2,ie+1
        vol1an(i+1) = g2an(i+1)*g31an(i+1)*x1an(i+1)
        dvl1an(i  ) = vfac*(vol1an(i+1) - vol1an(i))
40    continue
      do i = is-2, ie+1
       dvl1ani(i) = 1.0D0/dvl1an(i)
      enddo
!
      vol1bn(is-2) = g2bn(is-2)*g31bn(is-2)*x1bn(is-2)
      do 50 i=is-2,ie+1
        vol1bn(i+1) = g2bn(i+1)*g31bn(i+1)*x1bn(i+1)
        dvl1bn(i+1) = vfac*(vol1bn(i+1) - vol1bn(i))
50    continue
!
      do i = is-1, ie+2
       dvl1bni(i) = 1.0D0/dvl1bn(i)
      enddo
!
      return
      end
