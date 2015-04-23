!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          L C D M H A L O                  \\\\\\\\\\
!
!=======================================================================
!
      subroutine lcdmhalo
!
!
!  written by: Daniel Whalen 06-07-04
!
!  PURPOSE:  Computes the density field of a primordial LCDM H halo 
!            envelope Brian O'Shea simulated to have surrounded the 
!            first star. This routine also establishes the central 
!            stellar fluxes over the energy ranges relevant to envelope
!            heating and and chemistry and calculates the initial 
!            gravitational potential present in the system.  
!            
!            
! 
!
!  EXTERNALS:  none
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
       integer     :: i, j, k, index       
       real(rl)    :: logr, dntable(127), logrmin, logrmax, logd
       real(rl)    :: r1, r2, dlogr, rdef, rhocrit
       real(rl)    :: rho, rc, rci, q1, xi, n_central, r_core, omega   
       character*3 :: usrtag
       namelist / pgen     / usrtag
       usrtag    = 'usr'
       omega     = -3.0
       n_central = 1.0e06
       r_core    = 2.1e16
       if (myid .eq. 0) then
         open(unit=66, file='lcdmhalotable.dat', status='unknown')
         do i=1,127
           read(66,*) dntable(i)
         enddo
         close(unit=66)
         read  (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
       do i=1,127
         buf_in(i) = dntable(i)
       enddo
#endif /* MPI_USED */
       endif
#ifdef MPI_USED
       call MPI_BCAST( buf_in, 127, MPI_FLOAT &
                     , 0, comm3d, ierr )
       call MPI_BCAST( usrtag, 3, MPI_CHARACTER &
                     , 0, comm3d, ierr )
       if (myid .ne. 0) then
         do i=1,127
           dntable(i) = buf_in(i)
         enddo
       endif
#endif /* MPI_USED */
       logrmin = -5.8893
       logrmax =  3.4278
       dlogr   =  0.076370
       rhocrit =  1.747d-3 * mh 
       do k=ks-3,ke+3
         do j=js,je
         do i=is,ie
           logr = dlog10(x1b(i)/cmpc)
           if (logr .gt. logrmax) logr=logrmax
           if (logr .lt. logrmin) logr=logrmin
           index =min0(127,max0(1,idint((logr-logrmin)/dlogr)+1))
           r1    = logrmin + (index - 1)*dlogr
           r2    = logrmin + (index    )*dlogr
           rdef  = r2 - r1
           logd  = dntable(index)+(logr-r1) &
                 *(dntable(index+1)-dntable(index))/rdef
           d    (i,j,k)= (10.0**logd) * mh
!          d    (i,j,k)= 0.1 * mh
           if (d(i,j,k) .le. rhocrit) d(i,j,k) = rhocrit
           e    (i,j,k)= 1.5 * d(i,j,k) * boltz * 500.0 / mh
         enddo
         enddo
       enddo
!      rho         =  n_central * mh
!      rc          =  r_core
!      rci         =  1.0/r_core
!      do k = ks,ke
!        do j = js,je
!        do i = is,ie
!        xi           = rc - x1a(i)
!        q1           = dsign(0.5D0,xi)
!        d    (i,j,k) = rho*((0.5+q1)+(0.5-q1)*(x1a(i)*rci)**omega)
!        e    (i,j,k) = 1.5 * d(i,j,k) * boltz * 10.0 / mh
!        enddo
!        enddo
!      enddo  
! -- i faces
       nreq = 0
       nsub = nsub+1
       call bvald(3,3,0,0,0,0,d)
       call bvale(3,3,0,0,0,0,e)
       if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
! -- j faces
       nreq = 0
       nsub = nsub+1
       call bvald(0,0,3,3,0,0,d)
       call bvale(0,0,3,3,0,0,e)
       if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
! -- k faces
       nreq = 0
       nsub = nsub+1
       call bvald(0,0,0,0,3,3,d)
       call bvale(0,0,0,0,3,3,e)
       if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
!       call bvald(3,3,3,3,3,3,d)
!       call bvale(3,3,3,3,3,3,e)
      do k=ks-3,ke+3
        do j=js-2,je+2
        do i=is-1,ie-1
         if(d(i,j,k).eq.rhocrit.and.d(i+1,j,k).eq.rhocrit) then
           gp(i+1,j,k) = gp(i,j,k) 
         else
           gp(i+1,j,k) = gp(i,j,k) + gamm1 * 2.0 * &
                         (e(i+1,j,k) - e(i,j,k))/ &
                         (d(i,j,k) + d(i+1,j,k))
         endif 
        enddo
        enddo
      enddo
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////              L C D M H A L O              \\\\\\\\\\
!
!=======================================================================
!
