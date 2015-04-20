#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             L C D M H A L O 2             \\\\\\\\\\
c
c=======================================================================
c
      subroutine lcdmhalo2
c
c
c  written by: Daniel Whalen 06-07-04
c
c  PURPOSE:  Computes the density field of a primordial LCDM H halo 
c            envelope Brian O'Shea simulated to have surrounded the 
c            first star. This routine also establishes the central 
c            stellar fluxes over the energy ranges relevant to envelope
c            heating and and chemistry and calculates the initial 
c            gravitational potential present in the system.  
c            
c            
c  updated:  03.21.05 for BWO hi-res 1st star runs
c  updated:  11.23.06 full 9-species/3D upgrades
c
c  ported to ZEUS-MP 2.1 by DJW 11.29.06
c
c  EXTERNALS:  none
c
c-----------------------------------------------------------------------
c
      use real_prec
      use param
      use cons
      use config
      use root
      use chem
      use field
      use bndry
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE

      integer  :: i, j, k, index       

      real(rl) :: logr, logrmin, logrmax, r1, r2, dlogr,
     .            rdef, rhocrit, ovrdns, logd, omega_b, 
     .            h0,rho_igm,t_halo,r_trans, y

      real(rl), dimension(100)      :: dntable
      real(rl), dimension(in,jn,kn) :: dd

      namelist / pgen     / usrtag, ovrdns, t_halo, r_trans

      usrtag  = 'usr'
      ovrdns  = 1.0
      t_halo  = 500.0
      r_trans = 0.
      h0      = 0.7
      omega_b = 0.04

      if (myid .eq. 0) then
        open(unit=66, file='lcdmhalotable2.dat', status='unknown')
        do i=1,100
          read(66,*) dntable(i)
        enddo
        close(unit=66)
        read  (1, pgen)
        write (2, pgen)
#ifdef MPI_USED
        do i=1,100
          buf_in(i) = dntable(i)
        enddo
        buf_in(101) = ovrdns 
        buf_in(102) = t_halo 
        buf_in(103) = r_trans 
      endif
      call MPI_BCAST( buf_in, 103, MPI_FLOAT
     .              , 0, comm3d, ierr )
      call MPI_BCAST( usrtag, 3, MPI_CHARACTER
     .              , 0, comm3d, ierr )
      if (myid .ne. 0) then
        do i=1,127
          dntable(i) = buf_in(i)
        enddo
        ovrdns    = buf_in(101)
        t_halo    = buf_in(102)
        r_trans   = buf_in(103)
#endif /* MPI_USED */
      endif

      logrmin = -4.1769
      logrmax =  3.3837
      dlogr   =  0.076370
      rhocrit =  1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
      rho_igm =  omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh

      do k=ks,ke
        do j=js,je
        do i=is,ie

          logr = dlog10(x1b(i)/cmpc)

          if (logr .lt. logrmin) logr=logrmin

          index =min0(99,max0(1,idint((logr-logrmin)/dlogr)+1))

          r1    = logrmin + (index - 1)*dlogr
          r2    = logrmin + (index    )*dlogr
          rdef  = r2 - r1
          logd  = dntable(index)+(logr-r1)
     .          *(dntable(index+1)-dntable(index))/rdef


          if (fh .lt. 1.0) then
            d(i,j,k) =  (10.0**logd) * mh * 1.22
          else 
            d(i,j,k) =  (10.0**logd) * mh
          endif

          if (d(i,j,k) .le. rho_igm) d(i,j,k) = rho_igm
          tgas (i,j,k) = t_halo

#ifdef H
        if (fh .gt. 0.8) then
          abun(i,j,k,1) = fh  
          abun(i,j,k,2) = tiny
        else 
          abun(i,j,k,1) = 0.999898 * fh 
          abun(i,j,k,2) = 0.0001   * fh 
        endif
        abun(i,j,k,3) = abun(i,j,k,2)
#endif /* H */
#ifdef He
        abun(i,j,k,4) = (1. - fh + tiny) 
        abun(i,j,k,5) = tiny
        abun(i,j,k,6) = tiny
#endif /* He */
#ifdef H2
        abun(i,j,k,7) = tiny  
        if (fh .gt. 0.8) then
          abun(i,j,k,8) = tiny 
        else 
          abun(i,j,k,8) = 2.0d-06 * fh * 2.0  
        endif
        abun(i,j,k,9) = tiny
#endif /* H2 */
        e(i,j,k) =(abun(i,j,k,1)   + abun(i,j,k,2)    + abun(i,j,k,3)
#ifdef He
     .           + abun(i,j,k,4)/4.+ abun(i,j,k,5)/4. + abun(i,j,k,6)/4.
#endif /* He */
#ifdef H2
     .           + abun(i,j,k,7)   + abun(i,j,k,8)/2. + abun(i,j,k,9)/2.
#endif /* H2 */
     .         ) * boltz * tgas(i,j,k) * d(i,j,k) / (mh * gamm1)
        enddo
        enddo
      enddo
       


      if (ldimen .gt. 1) then
        call random_number(harvest=y)
        call random_number(dd)

        do k = ks-3,ke+3
          do j = js,je
          do i = is-1,ie-1

          if (x1b(i) .ge. r_trans) then
            d(i,j,k) = d(i,j,k) + 0.01 * d(i,j,k) * (dd(i,j,k) - 0.5)
          endif

 20       enddo
          enddo
        enddo 
      endif  ! ldimen = 1
 
c -- i faces
      nreq = 0
      nsub = nsub+1
      call bvald(3,3,0,0,0,0,d)
      call bvale(3,3,0,0,0,0,e)
#ifdef MPI_USED
      if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif /* MPI_USED */

c -- j faces
      nreq = 0
      nsub = nsub+1
      call bvald(0,0,3,3,0,0,d)
      call bvale(0,0,3,3,0,0,e)
#ifdef MPI_USED
      if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif /* MPI_USED */

c -- k faces
      nreq = 0
      nsub = nsub+1
      call bvald(0,0,0,0,3,3,d)
      call bvale(0,0,0,0,3,3,e)
#ifdef MPI_USED
      if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif /* MPI_USED */

      do k=ks-2,ke+2
        do j=js-2,je+2
        do i=is-1,ie-1

         if(d(i,j,k).eq.rho_igm.and.d(i+1,j,k).eq.rho_igm) then

           gp(i+1,j,k) = gp(i,j,k) 

         else

           gp(i+1,j,k) = gp(i,j,k) + gamm1 * 2.0 * 
     .                   (e(i+1,j,k) - e(i,j,k))/
     .                   (d(i,j,k) + d(i+1,j,k))

         endif 

        enddo
        enddo

      enddo

      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////            L C D M H A L O 2              \\\\\\\\\\
c
c=======================================================================
c


