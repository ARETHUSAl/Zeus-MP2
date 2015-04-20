#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             B N D R Y F L O W             \\\\\\\\\\
c
c
c           Update time-dependent boundary conditions
c
c
c           Developed by Bob van Veelen
c             Last updated on: 07-11-2007
c             modified for primordial SN by Daniel Whalen on 9-26-07
c=======================================================================
      subroutine bndryflow2
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use chem
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
c
      implicit none
c
      integer  :: i, j, k, z, l, index, nlines(12), nhalo,nbins
      real(rl) :: amtemp, vexpanmax

      real(rl) :: logr, logrmin, logrmax, r1, r2, r, 
     .            rdef, rhocrit, ovrdns, logd, omega_b, 
     .            h0,rho_igm,t_halo,r_trans, y, logT, logv,
     .            loge,logHI,logHII,logel,logHeI,logHeII,logHeIII,
     .            logHM,logH2I,logH2II

      real(rl), dimension(200) :: dtable,Ttable,vrtable,rtable,
     .   etable,HItable,HIItable,eltable,HeItable,HeIItable,HeIIItable,
     .   HMtable,H2Itable,H2IItable,gptable


      character*11 :: tablefile
    
      common /inflow/ vexpanmax,amtemp,dtable,Ttable,vrtable,rtable,
     .                nhalo,nlines,etable,HItable,HIItable,eltable,
     .                HeItable,HeIItable,HeIIItable,HMtable,H2Itable,
     .                H2IItable,gptable
c

      nbins = nlines(nhalo) - 1
      vexpanmax=10.
      if (vexpanmax .ne. 0.) then
c
c      print*,'boundary update'
c
c      call boundflow2D(filename,ni,nk)
c      call boundflow1D(filename,ni)
c
c     You can also specify constant values for the inflow


      h0      = 0.7
      omega_b = 0.04
      rhocrit = 1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
      rho_igm = omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh

      do z=1,2
        do k=ks,ke
          do j=js,je

          do l = 1,nbins
             if (x1b(ie+z)/cmpc .ge. rtable(l  )  .and. 
     .           x1b(ie+z)/cmpc .lt. rtable(l+1)) index = l
          enddo

          r = x1b(ie+z)/cmpc

          if (r .lt. rtable(1)                 ) then
            index = 1
            logd     = dtable    (index)
            loge     = etable    (index)
            logT     = Ttable    (index)
            logHI    = HItable   (index)
            logHII   = HIItable  (index)
            logel    = eltable   (index)
            logHeI   = HeItable  (index)
            logHeII  = HeIItable (index)
            logHeIII = HeIIItable(index)
            logHM    = HMtable   (index)
            logH2I   = H2Itable  (index)
            logH2II  = H2IItable (index)
            logv     = vrtable   (index)
          else if (r .gt. rtable(nlines(nhalo))) then
            index = nlines(nhalo)
            logd     = dtable    (index)
            loge     = etable    (index)
            logT     = Ttable    (index)
            logHI    = HItable   (index)
            logHII   = HIItable  (index)
            logel    = eltable   (index)
            logHeI   = HeItable  (index)
            logHeII  = HeIItable (index)
            logHeIII = HeIIItable(index)
            logHM    = HMtable   (index)
            logH2I   = H2Itable  (index)
            logH2II  = H2IItable (index)
            logv     = vrtable   (index)
          else
            logr  = dlog10(r)
            r1    = dlog10(rtable(index  ))
            r2    = dlog10(rtable(index+1))
            rdef  = r2 - r1
  
            logd     = dtable    (index  ) + (logr-r1)
     .               *(dtable    (index+1) - dtable    (index))/rdef

            loge     = etable    (index  ) + (logr-r1)
     .               *(etable    (index+1) - etable    (index))/rdef

            logT     = Ttable    (index  ) + (logr-r1)
     .               *(Ttable    (index+1) - Ttable    (index))/rdef

            logHI    = HItable   (index  ) + (logr-r1)
     .               *(HItable   (index+1) - HItable   (index))/rdef

            logHII   = HIItable  (index  ) + (logr-r1)
     .               *(HIItable  (index+1) - HIItable  (index))/rdef

            logel    = eltable   (index  ) + (logr-r1)
     .               *(eltable   (index+1) - eltable   (index))/rdef

            logHeI   = HeItable  (index  ) +(logr-r1)
     .               *(HeItable  (index+1) - HeItable  (index))/rdef

            logHeII  = HeIItable (index  ) + (logr-r1)
     .               *(HeIItable (index+1) - HeIItable (index))/rdef

            logHeIII = HeIIItable(index  ) + (logr-r1)
     .               *(HeIIItable(index+1) - HeIIItable(index))/rdef
  
            logHM    = HMtable   (index  ) + (logr-r1)
     .               *(HMtable   (index+1) - HMtable   (index))/rdef

            logH2I   = H2Itable  (index  ) + (logr-r1)
     .               *(H2Itable  (index+1) - H2Itable  (index))/rdef

            logH2II  = H2IItable (index  ) + (logr-r1)
     .               *(H2IItable (index+1) - H2IItable (index))/rdef
  
            logv     = vrtable   (index  ) + (logr-r1)
     .               *(vrtable   (index+1) - vrtable   (index))/rdef
          endif

          if (nhy .eq. 28860) print*,"logv is: ",logv

          doib(j,k,z)    = 10.0**logd 
          eoib(j,k,z)    = 10.0**loge 
          tgas(ie+z,j,k) = 10.0**logT
          v1oib(j,k,z)   = 10.0**logv * cmkm
          v2oib(j,k,z)   = tiny
          v3oib(j,k,z)   = tiny
          aboib(j,k,z,1) = 10.0**logHI
          aboib(j,k,z,2) = 10.0**logHII
          aboib(j,k,z,3) = 10.0**logel
          aboib(j,k,z,4) = 10.0**logHeI
          aboib(j,k,z,5) = 10.0**logHeII
          aboib(j,k,z,6) = 10.0**logHeIII
          aboib(j,k,z,7) = 10.0**logHM
          aboib(j,k,z,8) = 10.0**logH2I
          aboib(j,k,z,9) = 10.0**logH2II

c          if (nhy.eq.9999) then
c             print*,"z,r, and d are: ",z,r,doib(j,k,z),d(ie,j,k)
c          endif
c
c      doib (j,k,z)  = 1.0D-9*rho0
c            doib (j,k,z)  = d(ie,j,k)
c      print*,doib(j,k,z)
c
c

          enddo
        enddo
      enddo

c      print*,"fh is: ",fh
c      print*,"doib is: ",doib(js,ks,1)
c      print*,"doib is: ",doib(js,ks,2)
c      print*,"aboib is: ",aboib(js,ks,1,1),aboib(js,ks,1,2)
c      print*,"aboib is: ",aboib(js,ks,1,4),aboib(js,ks,1,5)

      do k=ks,ke  
        do j=js,je
          do z=1,2
          if(doib(j,k,z).le.0..or.doib(j,k,z).eq.tiny) 
     &         print*,'Error in source terms (doib)',i,j,k,doib(j,k,z)
          if(eoib(j,k,z).le.0..or.eoib(j,k,z).eq.tiny)
     &         print*,'Error in source terms (eoib)',i,j,k,eoib(j,k,z)
          enddo
        enddo
      enddo
c
      endif !vexpanmax
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             B N D R Y F L O W             \\\\\\\\\\
c
c=======================================================================
c
