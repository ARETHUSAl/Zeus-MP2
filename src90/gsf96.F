#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 G S F 9 6                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine gsf96
c
c
c  written by: Daniel Whalen 10-31-05
c
c  PURPOSE:    Sets up the the perturbed density field for the I-front
c              enhanced dynamical instability examined by Garcia-Segura
c              and Franco, 1996 (ApJ 469: 171-188).  Initial conditions
c              can be either primordial or galactic with a given 
c              metallicity relative to solar (z_sol, a chemcon pgen
c              variable).  Here, f_H2 = n_H2/n_H.  An initial electron
c              fraction X_e must be assumed.  If we assume the gas has
c              metallicity z_sol then all the C, Si, Fe, and S are 
c              singly-ionized by the background UV.  This determines the
c              initial electron fraction in a galactic environment.  If
c              the environment is instead primordial, we assume X_e = 
c              10^-4, which is expected at redshifts of 20 (e.g. Ricotti
c              et al ApJ 2001).  We can add to this whatever is 
c              contributed by the metals deemed present at that epoch.
c              
c  ported to ZEUS-MP 2.1 by DJW 11.29.06
c
c                        
c  LOCAL VARIABLES:
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

      integer  :: i, j, k 
      real(rl) :: rho, rc, rci, q1, xi, n_central, r_trans,
     .            t_backg, r_core, y, f_H2, C, Si, Fe, S, X_e, anglefac

      real(rl), dimension(in,jn,kn) :: dd 

      namelist / pgen     / usrtag , omega  , n_central, r_core, 
     .                      t_backg, r_trans, f_H2

      usrtag    =  'usr'
      omega     = -2.0
      n_central = 1.0e04
      r_core    = 6.1680e17
      t_backg   = 10.0
      r_trans   = 3.85e17

      C         = 3.75d-4 * z_sol
      Si        = 3.20d-5 * z_sol
      Fe        = 3.20d-5 * z_sol
      S         = 1.40d-5 * z_sol

      X_e       = 1.0e-4 !C + Si + Fe + S + 1.0e-4

      if (myid_w .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
        buf_in(1) = omega 
        buf_in(2) = n_central
        buf_in(3) = r_core
        buf_in(4) = t_backg
        buf_in(5) = r_trans
      endif
      call MPI_BCAST( buf_in, 5, MPI_FLOAT
     .               , 0, comm3d, ierr )
      call MPI_BCAST( usrtag, 3, MPI_CHARACTER
     .               , 0, comm3d, ierr )
      if (myid_w .ne. 0) then
        omega     = buf_in(1)
        n_central = buf_in(2)
        r_core    = buf_in(3)
        t_backg   = buf_in(4)
        r_trans   = buf_in(5)
#endif /* MPI_USED */
      endif

      if (fh .lt. 1.0) then
        rho =  n_central * mh * 1.22
      else 
        rho =  n_central * mh
      endif

      rc  =  r_core
      rci =  1.0/r_core

      do k = ks-2,ke+2
        do j = js,je
        do i = is,ie
 
c        if      (x3a(k) .ge. 0.      .and. x3a(k) .lt.    pi/4) then
c           anglefac = 1.
c        else if (x3a(k) .ge.    pi/4 .and. x3a(k) .lt.    pi/2) then
c           anglefac = 0.1
c        else if (x3a(k) .ge.    pi/2 .and. x3a(k) .lt. 3.*pi/4) then
c           anglefac = 0.01
c        else if (x3a(k) .ge. 3.*pi/4 .and. x3a(k) .lt.    pi  ) then
c           anglefac = 0.001
c        else if (x3a(k) .ge.    pi   .and. x3a(k) .lt. 5.*pi/4) then
c           anglefac = 1.0d-4
c        else if (x3a(k) .ge. 5.*pi/4 .and. x3a(k) .lt. 3.*pi/2) then
c           anglefac = 1.0d-5
c        else if (x3a(k) .ge. 3.*pi/2 .and. x3a(k) .lt. 7.*pi/4) then
c           anglefac = 1.0d-6
c        else if (x3a(k) .ge. 7.*pi/4 .and. x3a(k) .lt. 2.*pi  ) then
c           anglefac = 1.0d-7
c        endif

c        print*,"anglefac is: ",anglefac

        xi           = rc - x1a(i)
        q1           = dsign(0.5D0,xi)

        d    (i,j,k) = rho*((0.5+q1)+(0.5-q1)*(x1a(i)*rci)**omega)
c     .  * anglefac
        tgas (i,j,k) = t_backg
#ifdef H
        if (fh .gt. 0.8) then
          abun(i,j,k,1) = fh  
          abun(i,j,k,2) = tiny
        else 
          abun(i,j,k,1) = fh / (1 + 2.*f_H2 + X_e)
          abun(i,j,k,2) = X_e * abun(i,j,k,1) 
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
          abun(i,j,k,8) = 2.0 * f_H2 * abun(i,j,k,1) 
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


 10     enddo
        enddo
      enddo  

      if (ldimen .gt. 1) then
        call random_number(harvest=y)
        call random_number(dd)

        do k = ks-2,ke+2
          do j = js,je
          do i = is-1,ie-1

          if (x1b(i) .ge. r_trans) then
            d(i,j,k) = d(i,j,k) + 0.01 * d(i,j,k) * (dd(i,j,k) - 0.5)
          endif

 20       enddo
          enddo
        enddo 
      endif   ! ldimen > 1
 
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
        do j=js-3,je+3
        do i=is-1,ie-1
          
         gp(i+1,j,k) = gp(i,j,k) + gamm1 * 2.0 * 
     .                   (e(i+1,j,k) - e(i,j,k))/
     .                   (d(i,j,k) + d(i+1,j,k))
          
        enddo
        enddo
      enddo

      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 G S F 9 6                 \\\\\\\\\\
c
c=======================================================================
c
