#include "rtchem.def"
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 O I _ 3 D                 \\\\\\\\\\
!
!=======================================================================
!
      subroutine oi_3D
!
!
!  written by: Daniel Whalen 05-25-06
!
!  ported to ZEUS-MP 2.1 by DJW 11.29.06
!
!  PURPOSE:  Computes the density field of a primordial LCDM halo 
!            envelope surrounding the first star, centering it at
!            (xc,yc,zc) in a cartesian box.  A radiation wave enters 
!            the box from the left, sweeps past the 3D halo, and begins 
!            to photoevaporate it. If ihalo=1, the halo is a spherically 
!            averaged radially-symmetric profile computed by Brian 
!            O'Shea with Enzo.  If ihalo=2, a TIS density profile from
!            Shapiro, Iliev, and Raga (MNRAS 1999) is instead used. Note
!            that if we specify an attenuated plane wave (iPWA=1), then
!            r_sep must be greater than or equal to the distance of the
!            halo center from the left edge of the box.
!            
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
      integer      :: i, j, k, index, nbins, ihalo  
      character*20 :: tablefile    
      real(rl)     :: logr, logrmin, logrmax, logd, mu, &
                      r1, r2, dlogr, rdef, rhocrit, r, &
                      ovrdns, omega_b, h0, rho_igm, yc, &
                      t_halo, zc, omega_0, a1, a2, c1, &
                      rho_0, c2, xi, r_trunc, m_halo, t_TIS
       real(rl), dimension(45 ) :: dn1table
       real(rl), dimension(46 ) :: dn2table
       real(rl), dimension(47 ) :: dn3table
       real(rl), dimension(48 ) :: dn4table
       real(rl), dimension(100) :: dntable
       real(rl), dimension(122) :: dn5table
       namelist / pgen     / usrtag, ovrdns, t_halo, xc, yc, zc, nbins, &
                             logrmin, logrmax, mu, ihalo, m_halo
       usrtag    = '002'
       ovrdns    = 1.0
       t_halo    = 400.0
       h0        = 0.7
       omega_b   = 0.043
       omega_0   = 0.27
       xc        = 0.0
       yc        = 0.0
       zc        = 0.0
       nbins     = 99
       mu        = 1.0
       ihalo     = 1
       m_halo    = 4.0d5
       if (myid .eq. 0) then
         read  (1, pgen)
         write (2, pgen)
         write(tablefile,"(a13,a3,a4)") 'lcdmhalotable' &
         ,usrtag,'.dat'
         open(unit=66, file=tablefile, status='unknown')
         if      (nbins .eq. 44) then
         do i=1,nbins+1
           read(66,*) dn1table(i)
         enddo
         else if (nbins .eq. 45) then
         do i=1,nbins+1
           read(66,*) dn2table(i)
         enddo
         else if (nbins .eq. 46) then
         do i=1,nbins+1
           read(66,*) dn3table(i)
         enddo
         else if (nbins .eq. 47) then
         do i=1,nbins+1
           read(66,*) dn4table(i)
         enddo
         else if (nbins .eq. 122) then
         do i=1,nbins+1
           read(66,*) dn5table(i)
         enddo
         else if (nbins .eq. 99) then
         do i=1,nbins+1
           read(66,*) dntable(i)
         enddo
         endif
         close(unit=66)
#ifdef MPI_USED
         buf_in(1) = ovrdns 
         buf_in(2) = t_halo 
         buf_in(3) = xc
         buf_in(4) = yc
         buf_in(5) = zc
         buf_in(6) = logrmin
         buf_in(7) = logrmax
         buf_in(8) = mu
         buf_in(9) = m_halo
         ibuf_in(1) = nbins
         ibuf_in(2) = ihalo
       endif
       call MPI_BCAST( buf_in, 9, MPI_FLOAT &
                     , 0, comm3d, ierr )
       call MPI_BCAST( usrtag, 3, MPI_CHARACTER &
                     , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 2, MPI_INTEGER &
                     , 0, comm3d, ierr )
       if      (nbins .eq. 44) then
        call MPI_BCAST( dn1table, 45, MPI_FLOAT &
                   , 0, comm3d, ierr )
       else if (nbins .eq. 45) then
        call MPI_BCAST( dn2table, 46, MPI_FLOAT &
                   , 0, comm3d, ierr )
       else if (nbins .eq. 46) then
        call MPI_BCAST( dn3table, 47, MPI_FLOAT &
                   , 0, comm3d, ierr )
       else if (nbins .eq. 47) then
        call MPI_BCAST( dn4table, 48, MPI_FLOAT &
                   , 0, comm3d, ierr )
       else if (nbins .eq. 122) then
        call MPI_BCAST( dn5table, 123, MPI_FLOAT &
                   , 0, comm3d, ierr )
       else if (nbins .eq. 99) then
        call MPI_BCAST( dntable, 100, MPI_FLOAT &
                   , 0, comm3d, ierr )
       endif
       if (myid .ne. 0) then
         ovrdns    = buf_in(1)
         t_halo    = buf_in(2)
         xc        = buf_in(3)
         yc        = buf_in(4)
         zc        = buf_in(5)
         logrmin   = buf_in(6)
         logrmax   = buf_in(7)
         mu        = buf_in(8)
         m_halo    = buf_in(9)
         nbins     = ibuf_in(1)
         ihalo     = ibuf_in(2)
#endif /* MPI_USED */
       endif
       rhocrit = 1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
       rho_igm = omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh
       if (ihalo .eq. 2) then
         a1 = 21.38
         a2 =  9.08
         c1 = 19.81
         c2 = 14.62
         r_trunc = 3.1549d20 * (m_halo/2.0d5)**(1./3.)*21./(1+rdshft) &
                 * (0.7/h0)**(2./3.) * (0.27/omega_0)**(1./3.)
         t_TIS   = 593.5 * (mu/1.22) * (h0/0.7)**(2./3.) &
                 * (m_halo/2.0d5)**(2./3.) * (1+rdshft)/21. &
                 * (omega_0/0.27)**(1./3.)
         rho_0   = 4.144d-22 * (omega_b/0.27) * (h0/0.7)**2 &
                 * ((1+rdshft)/21.)**3
!         print*,"r_t, T, rho_0 are: ",r_trunc/cmpc,t_TIS,dlog10(rho_0)
       endif
       do k=ks,ke
         do j=js,je
         do i=is,ie
         if (ihalo .eq. 1) then
           r     = dsqrt((dabs(x1b(i))-xc)**2+(dabs(x2b(j))-yc)**2 &
                 +       (dabs(x3b(k))-zc)**2)/cmpc
           dlogr = (logrmax-logrmin)/nbins
           logr  = dlog10(r)
           if (logr .lt. logrmin) logr=logrmin
           index =min0(nbins,max0(1,idint((logr-logrmin)/dlogr)+1))
           r1    = logrmin + (index - 1)*dlogr
           r2    = logrmin + (index    )*dlogr
           rdef  = r2 - r1
           if      (nbins .eq. 44) then
           logd  = dn1table(index)+(logr-r1) &
                 *(dn1table(index+1)-dn1table(index))/rdef
           else if (nbins .eq. 45) then
           logd  = dn2table(index)+(logr-r1) &
                 *(dn2table(index+1)-dn2table(index))/rdef
           else if (nbins .eq. 46) then
           logd  = dn3table(index)+(logr-r1) &
                 *(dn3table(index+1)-dn3table(index))/rdef
           else if (nbins .eq. 47) then
           logd  = dn4table(index)+(logr-r1) &
                 *(dn4table(index+1)-dn4table(index))/rdef
           else if (nbins .eq. 122) then
           logd  = dn5table(index)+(logr-r1) &
                 *(dn5table(index+1)-dn5table(index))/rdef
           else if (nbins .eq. 99) then
           logd  = dntable(index)+(logr-r1) &
                 *(dntable(index+1)-dntable(index))/rdef
           endif
           d   (i,j,k) = (10.0**logd ) * mh * mu
           tgas(i,j,k) = t_halo
           if (d(i,j,k) .le. rho_igm) d(i,j,k) = rho_igm
         else if (ihalo .eq. 2) then
	 if      (lgeom .eq. 1) then
	   if      (ldimen .eq. 1) then
             r = dsqrt((dabs(x1b(i))-xc)**2)
	   else if (ldimen .eq. 3) then
             r = dsqrt((dabs(x1b(i))-xc)**2+(dabs(x2b(j))-yc)**2 &
               +       (dabs(x3b(k))-zc)**2)
	   endif
	 else if (lgeom .eq. 2) then
	   if      (ldimen .eq. 2) then
             r = dsqrt((dabs(x1b(i)))**2+(dabs(x2b(j)))**2)
	   endif
	 else if (lgeom .eq. 3) then
	     r = x1b(i)
	 endif  
!           if (r .gt. r_trunc) then 
!             d(i,j,k) =  rho_igm
!           else if (r .le. r_trunc) then
             xi       =  r/(r_trunc/29.4)
             d(i,j,k) = (a1/(a2+xi**2) - c1/(c2+xi**2)) * rho_0
!           endif
           tgas(i,j,k) = t_TIS
         endif
         abun(i,j,k,1)= 0.999898 * fh   
         abun(i,j,k,2)= 1.D-4    * fh 
         abun(i,j,k,3)= abun(i,j,k,2)
#ifdef He
         abun(i,j,k,4)= (1.0-fh) + tiny 
         abun(i,j,k,5)= tiny
         abun(i,j,k,6)= tiny
#endif /* He */
#ifdef H2
         abun(i,j,k,7)= tiny  
         abun(i,j,k,8)= 2.0D-6   * fh * 2.0
         abun(i,j,k,9)= tiny  
#endif /* H2 */
         e(i,j,k)=(abun(i,j,k,1)   + abun(i,j,k,2)    + abun(i,j,k,3)
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
!#ifdef MPI_USED
!      call MPI_BARRIER(comm3d, ierr)
!#endif /* MPI_USED */
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 O I _ 3 D                 \\\\\\\\\\
!
!=======================================================================
!
