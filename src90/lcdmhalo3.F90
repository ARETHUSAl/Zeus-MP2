#include "rtchem.def"
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////             L C D M H A L O 3             \\\\\\\\\\
!
!=======================================================================
!
      subroutine lcdmhalo3
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
!  updated:  03.21.05 for BWO hi-res 1st star runs
!  updated:  11.23.06 full 9-species/3D upgrades
!  updated:  08.27.07 Brian's new halos
!
!  ported to ZEUS-MP 2.1 by DJW 11.29.06
!
!  EXTERNALS:  none
!
!-----------------------------------------------------------------------
!
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
!
      implicit NONE
      integer  :: i, j, k, l, index, nlines(3), nbins, nhalo       
      real(rl) :: logr, logrmin, logrmax, r1, r2, r, a, &
                  rdef, rhocrit, ovrdns, logd, omega_b, &
                  h0,rho_igm,t_halo,r_trans, y, logT, logv
      real(rl), dimension(128) :: dtable,Ttable,vrtable,rtable
      real(rl), dimension(35 ) :: dg,eg,tg,HIg,HIIg,deg,HeIg,HeIIg, &
                                  HeIIIg,HMg,H2Ig,H2IIg,vi
      data nlines(1),nlines(2),nlines(3) / 128,128,127 /
      real(rl), dimension(in,jn,kn) :: dd
      character*11 :: tablefile
      namelist / pgen     / usrtag, ovrdns, r_trans, nhalo
      usrtag  = 'usr'
      ovrdns  = 1.0
      r_trans = 0.
      nhalo   = 1
      h0      = 0.7
      omega_b = 0.04
 78   format(13(1pe11.3e3,1x))
 79   format( 2(1pe11.3e3,1x))
      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)
        nbins = nlines(nhalo)-1
        write(tablefile,"(a5,i2.2,a4)") 'halo0',nhalo,'.dat'
        open(unit=66, file=tablefile, status='unknown')
        read(66,*) logrmin
        do i=1,nbins+1
          read(66,*) rtable(i),dtable(i),Ttable(i),vrtable(i)
        enddo
        read(66,*) logrmax
        close(unit=66)
        if (ldimen .gt. 1) then
        open(unit=77, file='initial_profile.dat',status='unknown')
        do i=1,14
          read(77,78) a,dg(i),eg(i),tg(i),HIg(i),HIIg(i),deg(i), &
                      HeIg(i),HeIIg(i),HeIIIg(i),HMg(i),H2Ig(i), &
                      H2IIg(i)
        enddo
        read(77,* )
        do i=1,14
          read(77,79) a, vi(i)
          vi(i) = vi(i) * cmkm
        enddo
        close(unit=77)
        endif
#ifdef MPI_USED
        buf_in (1) = ovrdns 
        buf_in (2) = r_trans 
        buf_in (3) = logrmin 
        buf_in (4) = logrmax 
        ibuf_in(1) = nbins
        ibuf_in(2) = nhalo
      endif
      call MPI_BCAST( buf_in, 4, MPI_FLOAT &
                    , 0, comm3d, ierr )
      call MPI_BCAST( ibuf_in,2, MPI_INTEGER &
                    , 0, comm3d, ierr )
      call MPI_BCAST( usrtag, 3, MPI_CHARACTER &
                    , 0, comm3d, ierr )
      call MPI_BCAST( dtable , 128, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( Ttable , 128, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( vrtable, 128, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( rtable , 128, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( dg    , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( eg    , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( tg    , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HIg   , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HIIg  , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( deg   , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HeIg  , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HeIIg , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HeIIIg, 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( HMg   , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( H2Ig  , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( H2IIg , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      call MPI_BCAST( vi    , 35, MPI_FLOAT &
                 , 0, comm3d, ierr )
      if (myid .ne. 0) then
        ovrdns    = buf_in (1)
        r_trans   = buf_in (2)
        logrmin   = buf_in (3)
        logrmax   = buf_in (4)
        nbins     = ibuf_in(1)
        nhalo     = ibuf_in(2)
#endif /* MPI_USED */
      endif
      rhocrit =  1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
      rho_igm =  omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh
      do k=ks,ke
        do j=js,je
        do i=is,ie
          do l = 1,nbins
             if (x1b(i)/cmpc .ge. rtable(l  )  .and. &
                 x1b(i)/cmpc .lt. rtable(l+1)) index = l
          enddo
          r = x1b(i)/cmpc
          if (r .lt. rtable(1)            ) index = 1
          if (r .gt. rtable(nlines(nhalo))) index = nlines(nhalo)
!          print*,"index is: ", index
          logr  = dlog10(r)
          r1    = dlog10(rtable(index  ))
          r2    = dlog10(rtable(index+1))
          rdef  = r2 - r1
          logd  = dtable(index  )+(logr-r1) &
                *(dtable(index+1)-dtable(index))/rdef
          logT  = Ttable(index  )+(logr-r1) &
                *(Ttable(index+1)-Ttable(index))/rdef
          logv  = vrtable(index  )+(logr-r1) &
                *(vrtable(index+1)-vrtable(index))/rdef
          d   (i,j,k) = (10.0**logd ) * mh * 1.22
          tgas(i,j,k) =  10.0**logT
          if (d(i,j,k) .le. rho_igm) d(i,j,k) = rho_igm
#ifdef H
          abun(i,j,k,1) = 0.9999   * fh 
          abun(i,j,k,2) = 0.0001   * fh 
          abun(i,j,k,3) = abun(i,j,k,2)
#endif /* H */
#ifdef He
          abun(i,j,k,4) = (1. - fh + tiny) 
          abun(i,j,k,5) = tiny
          abun(i,j,k,6) = tiny
#endif /* He */
#ifdef H2
          abun(i,j,k,7) = tiny  
          abun(i,j,k,8) = tiny !2.0d-06 * fh * 2.0  
          abun(i,j,k,9) = tiny
#endif /* H2 */
        e(i,j,k) =(abun(i,j,k,1)   + abun(i,j,k,2)    + abun(i,j,k,3)
#ifdef He &
                 + abun(i,j,k,4)/4.+ abun(i,j,k,5)/4. + abun(i,j,k,6)/4.
#endif /* He */
#ifdef H2 &
                 + abun(i,j,k,7)   + abun(i,j,k,8)/2. + abun(i,j,k,9)/2.
#endif /* H2 */ &
               ) * boltz * tgas(i,j,k) * d(i,j,k) / (mh * gamm1)
        enddo
        enddo
      enddo
      if (ldimen .gt. 1) then
      do k=ks,ke
      do j=js,je
      do i=is,is+13
         d   (i,j,k  ) = dg    (i-is+1)
         e   (i,j,k  ) = eg    (i-is+1)
         tgas(i,j,k  ) = tg    (i-is+1)
         abun(i,j,k,1) = HIg   (i-is+1)
         abun(i,j,k,2) = HIIg  (i-is+1)
         abun(i,j,k,3) = deg   (i-is+1)
         abun(i,j,k,4) = HeIg  (i-is+1)
         abun(i,j,k,5) = HeIIg (i-is+1)
         abun(i,j,k,6) = HeIIIg(i-is+1)
         abun(i,j,k,7) = HMg   (i-is+1)
         abun(i,j,k,8) = H2Ig  (i-is+1)
         abun(i,j,k,9) = H2IIg (i-is+1)
         v1  (i,j,k  ) = vi    (i-is+1)
      enddo
      enddo
      enddo
      endif
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
! -- i faces
      nreq = 0
      nsub = nsub+1
      call bvald(3,3,0,0,0,0,d)
      call bvale(3,3,0,0,0,0,e)
#ifdef MPI_USED
      if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif /* MPI_USED */
! -- j faces
      nreq = 0
      nsub = nsub+1
      call bvald(0,0,3,3,0,0,d)
      call bvale(0,0,3,3,0,0,e)
#ifdef MPI_USED
      if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif /* MPI_USED */
! -- k faces
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
!    //////////            L C D M H A L O 3              \\\\\\\\\\
!
!=======================================================================
!
