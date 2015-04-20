#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 O I _ 3 D                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine oi_3D
c
c
c  written by: Daniel Whalen 05-25-06
c
c  ported to ZEUS-MP 2.1 by DJW 11.29.06
c
c  PURPOSE:  Computes the density field of a primordial LCDM halo 
c            satellite to the first star, centering it at
c            (xc,yc,zc) in a cartesian box.  A radiation wave enters 
c            the box from the left, sweeps past the 3D halo, and begins 
c            to photoevaporate it. If ihalo=1, the halo is a spherically 
c            averaged radially-symmetric profile computed by Brian 
c            O'Shea with Enzo.  If ihalo=2, a TIS density profile from
c            Shapiro, Iliev, and Raga (MNRAS 1999) is instead used. Note
c            that if we specify an attenuated plane wave (iPWA=1), then
c            r_sep must be greater than or equal to the distance of the
c            halo center from the left edge of the box.
c            
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

      integer      :: i,j,k,index,nlines(99),nbins,ihalo,nhalo  

      character*11 :: tablefile    
       

      real(rl)     :: logr, logrmin, logrmax, logd, mu,
     .                r1, r2, dlogr, rdef, rhocrit, r, 
     .                ovrdns, omega_b, h0, rho_igm, yc, 
     .                t_halo, zc, omega_0, a1, a2, c1, 
     .                rho_0, c2, xi, r_trunc, m_halo, 
     .                t_TIS, logT, logv, f_H2, X_e

      real(rl), dimension(120) :: dtable,Ttable,vrtable

      data nlines(11),nlines(12),nlines(13),nlines(14),nlines(15),
     .     nlines(16),nlines(17),nlines(18),nlines(19),nlines(20),
     .     nlines(21),nlines(22),nlines(23),nlines(39),nlines(59), 
     .     nlines(73) / 54,62,66,70,73,78,82,84,90,92,98,101,
     .                  50,46,50,54  /
      namelist / pgen     / nhalo, ovrdns, xc, yc, zc,
     .                      mu, ihalo, m_halo, f_H2, X_e

      nhalo   = 10
      ovrdns  = 1.0
      h0      = 0.7
      omega_b = 0.043
      omega_0 = 0.27
      xc      = 0.0
      yc      = 0.0
      zc      = 0.0
      mu      = 1.0
      ihalo   = 1
      m_halo  = 2.0d5
      X_e     = 1.0e-04
      f_H2    = 2.0e-06

 10   format(1p4e11.3)

      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)

	nbins = nlines(nhalo)-1

        write(tablefile,"(a5,i2.2,a4)") 'halo0',nhalo,'.dat'

        open(unit=66, file=tablefile, status='unknown')

	read(66,*) logrmin
        do i=1,nbins+1
          read(66,*) dtable(i),Ttable(i),vrtable(i)
        enddo
	read(66,*) logrmax

        close(unit=66)

#ifdef MPI_USED
        buf_in(1) = ovrdns 
        buf_in(2) = xc
        buf_in(3) = yc
        buf_in(4) = zc
        buf_in(5) = logrmin
        buf_in(6) = logrmax
        buf_in(7) = mu
        buf_in(8) = m_halo
        buf_in(9) = X_e
        buf_in(10)= f_H2

        ibuf_in(1) = nbins
        ibuf_in(2) = ihalo
        ibuf_in(3) = nhalo
      endif
       
      call MPI_BCAST( buf_in,10, MPI_FLOAT
     .              , 0, comm3d, ierr )
      call MPI_BCAST( ibuf_in,3, MPI_INTEGER
     .              , 0, comm3d, ierr )

      if (myid .ne. 0) then
        ovrdns    = buf_in(1)
        xc        = buf_in(2)
        yc        = buf_in(3)
        zc        = buf_in(4)
        logrmin   = buf_in(5)
        logrmax   = buf_in(6)
        mu        = buf_in(7)
        m_halo    = buf_in(8)
        X_e       = buf_in(9)
        f_H2      = buf_in(10)

        nbins     = ibuf_in(1)
        ihalo     = ibuf_in(2)
        nhalo     = ibuf_in(3)

#endif /* MPI_USED */
      endif

#ifdef MPI_USED

      call MPI_BCAST( dtable , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( Ttable , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( vrtable, nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )

#endif /* MPI_USED */

      rhocrit = 1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
      rho_igm = omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh

      if (ihalo .eq. 2) then

        a1 = 21.38
        a2 =  9.08
        c1 = 19.81
        c2 = 14.62

        r_trunc = 3.1549d20 * (m_halo/2.0d5)**(1./3.)*21./(1+rdshft)
     .          * (0.7/h0)**(2./3.) * (0.27/omega_0)**(1./3.)
        t_TIS   = 593.5 * (mu/1.22) * (h0/0.7)**(2./3.)  
     .          * (m_halo/2.0d5)**(2./3.) * (1+rdshft)/21.
     .          * (omega_0/0.27)**(1./3.)
        rho_0   = 4.144d-22 * (omega_b/0.27) * (h0/0.7)**2 
     .          * ((1+rdshft)/21.)**3
c        print*,"r_t, T, rho_0 are: ",r_trunc/cmpc,t_TIS,dlog10(rho_0)

      endif

      do k=ks,ke
        do j=js,je
        do i=is,ie

  	  if      (lgeom .eq. 1) then
	    if      (ldimen .eq. 1) then
              r = dsqrt((dabs(x1b(i))-xc)**2)
	    else if (ldimen .eq. 3) then
              r = dsqrt((dabs(x1b(i))-xc)**2+(dabs(x2b(j))-yc)**2
     .          +       (dabs(x3b(k))-zc)**2)
	    endif
 	  else if (lgeom .eq. 2) then
	    if      (ldimen .eq. 2) then
              r = dsqrt((dabs(x1b(i)))**2+(dabs(x2b(j)))**2)
	    endif
	  else if (lgeom .eq. 3) then
	      r = x1b(i)
	  endif 
 
          if (ihalo .eq. 1) then

            dlogr = (logrmax-logrmin)/nbins
            logr  = dlog10(r/cmpc)

            if (logr .lt. logrmin) logr=logrmin

            index =min0(nbins,max0(1,idint((logr-logrmin)/dlogr)+1))

            r1    = logrmin + (index - 1)*dlogr
            r2    = logrmin + (index    )*dlogr
            rdef  = r2 - r1

            logd  = dtable(index)+(logr-r1)
     .            *(dtable(index+1)-dtable(index))/rdef

            logT  = Ttable(index)+(logr-r1)
     .            *(Ttable(index+1)-Ttable(index))/rdef

            logv  = vrtable(index)+(logr-r1)
     .            *(vrtable(index+1)-vrtable(index))/rdef


            d   (i,j,k) = (10.0**logd ) * mh * mu
            tgas(i,j,k) =  10.0**logT

            if (d(i,j,k) .le. rho_igm) d(i,j,k) = rho_igm

          else if (ihalo .eq. 2) then

c           if (r .gt. r_trunc) then 
c             d(i,j,k) =  rho_igm
c           else if (r .le. r_trunc) then
              xi       =  r/(r_trunc/29.4)
              d(i,j,k) = (a1/(a2+xi**2) - c1/(c2+xi**2)) * rho_0
c            endif
            tgas(i,j,k) = t_TIS

          endif

          abun(i,j,k,1)= fh / (1 + 2.*f_H2 + X_e)   
          abun(i,j,k,2)= X_e * abun(i,j,k,1) 
          abun(i,j,k,3)= abun(i,j,k,2)
#ifdef He
          abun(i,j,k,4)= (1.0-fh) + tiny 
          abun(i,j,k,5)= tiny
          abun(i,j,k,6)= tiny
#endif /* He */
#ifdef H2
          abun(i,j,k,7)= tiny  
          abun(i,j,k,8)= 2.0 * f_H2 * abun(i,j,k,1)
          abun(i,j,k,9)= tiny  
#endif /* H2 */
         e(i,j,k)=(abun(i,j,k,1)   + abun(i,j,k,2)    + abun(i,j,k,3)
#ifdef He
     .           + abun(i,j,k,4)/4.+ abun(i,j,k,5)/4. + abun(i,j,k,6)/4.
#endif /* He */
#ifdef H2
     .           + abun(i,j,k,7)   + abun(i,j,k,8)/2. + abun(i,j,k,9)/2.
#endif /* H2 */
     .           ) * boltz * tgas(i,j,k) * d(i,j,k) / (mh * gamm1)
         enddo
         enddo
       enddo
       
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

c#ifdef MPI_USED
c      call MPI_BARRIER(comm3d, ierr)
c#endif /* MPI_USED */

       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 O I _ 3 D                 \\\\\\\\\\
c
c=======================================================================
c


