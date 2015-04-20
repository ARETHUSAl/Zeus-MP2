#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              S U P E R N O V A            \\\\\\\\\\
c
c
c           Routine to set up a supernova explosion problem
c
c
c           Developed by Bob van Veelen
c             Last updated on: 07-11-2007
c             modified for primordial SN by Daniel Whalen on 9-25-07
c=======================================================================
      subroutine supernova2
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use chem
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
      integer  :: i, j, k, ncells, procn, sourceprob
      real(rl) :: Eej,Mej,n,rmax,rmin,vmax
      real(rl) :: rho0, v0,mmwold
      real(rl) :: tset, sntemp, amtemp, vexpanmax

      real(rl)     :: dtemp(in,jn), etemp(in,jn), v1temp(in,jn)
      real(rl)     :: v2temp(in,jn)
      character*99 :: filename
      integer      :: ni, nk

      integer  :: l, index, nlines(12), nbins, nhalo       

      real(rl) :: logr, logrmin, logrmax, r1, r2, r, 
     .            rdef, rhocrit, ovrdns, logd, omega_b, 
     .            h0,rho_igm,t_halo,r_trans, y, logT, logv,
     .            loge,logHI,logHII,logel,logHeI,logHeII,logHeIII,
     .            logHM,logH2I,logH2II,loggp

      real(rl), dimension(200) :: dtable,Ttable,vrtable,rtable,
     .   etable,HItable,HIItable,eltable,HeItable,HeIItable,HeIIItable,
     .   HMtable,H2Itable,H2IItable,gptable

      data nlines(1 ),nlines(2 ),nlines(3 ),nlines(4 ),
     .     nlines(5 ),nlines(6 ),nlines(7 ),nlines(8 ),
     .     nlines(9 ),nlines(10),nlines(11),nlines(12)  
     .     / 200,200,200,200,200,200,200,200,200,200,200,200 /

      character*11 :: tablefile
    
c      real(rl)     :: dtemp(ie), etemp(ie), v1temp(ie)

      namelist / pgen     / Eej, Mej, n, vmax, sntemp, amtemp, nhalo,
     .                      ovrdns

      common /inflow/ vexpanmax,amtemp,dtable,Ttable,vrtable,rtable,
     .                nhalo, nlines,etable,HItable,HIItable,eltable,
     .                HeItable,HeIItable,HeIIItable,HMtable,H2Itable,
     .                H2IItable,gptable
      common /sourceterms/   sourceprob
c
      sourceprob = 1
c


c
c-----Give in supernova parameters and possible input file parameters----
c
      Eej    =  2.0D51
      Mej    =  13.52
      n      =  9                   ! slope of density drop
      nhalo  =  1
      vmax   =  3.0D9               ! max velocity of ejecta
      sntemp =  1.0D6
      amtemp =  1.0D2
c
c     read in data for precomputed circumstellar medium
c
c      filename = '/home/veelen/zeusmp2/inputdata/zu026da'
c      ni       = 1000
c      nk       = 200
c

      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)

	nbins = nlines(nhalo)-1

        write(tablefile,"(a5,i2.2,a4)") 'halo0',nhalo,'.dat'

        open(unit=66, file=tablefile, status='unknown')

        do i=1,nbins+1
          read(66,*) rtable(i),dtable(i),etable(i),Ttable(i),
     .               HItable(i),HIItable(i),eltable(i),HeItable(i),
     .               HeIItable(i),HeIIItable(i),HMtable(i),
     .               H2Itable(i),H2IItable(i),gptable(i),vrtable(i)
        enddo

        close(unit=66)

        logrmin = rtable(1  )
        logrmax = rtable(200)

c        print*,"gptable is: ", (gptable(i), i=1,200)
c        stop

#ifdef MPI_USED
        buf_in(1) = Eej 
        buf_in(2) = Mej
        buf_in(3) = vmax
        buf_in(4) = sntemp
        buf_in(5) = amtemp
        buf_in(6) = ovrdns 
        buf_in(7) = logrmin 
        buf_in(8) = logrmax 

        ibuf_in(1)= n
        ibuf_in(2)= nhalo
        ibuf_in(3)= nbins

      endif
      call MPI_BCAST( buf_in, 8, MPI_FLOAT
     &              , 0, comm3d, ierr )
      call MPI_BCAST( ibuf_in,3, MPI_INTEGER
     .              , 0, comm3d, ierr )
      call MPI_BCAST( dtable    , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( etable    , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( Ttable    , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HItable   , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HIItable  , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( eltable   , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HeItable  , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HeIItable , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HeIIItable, nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( HMtable   , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( H2Itable  , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( H2IItable , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( vrtable   , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )
      call MPI_BCAST( gptable   , nbins+1, MPI_FLOAT
     .           , 0, comm3d, ierr )

      if (myid .ne. 0) then
        Eej     = buf_in (1)
        Mej     = buf_in (2)
	vmax    = buf_in (3)
        sntemp  = buf_in (4)
        amtemp  = buf_in (5)
        ovrdns  = buf_in (6)
        logrmin = buf_in (7)
        logrmax = buf_in (8)

        n       = ibuf_in(1)
        nhalo   = ibuf_in(2)
        nbins   = ibuf_in(3)
#endif /* MPI_USED */
      endif

      Mej = Mej * msol

      ncells =  int(0.8*(ie-is))
      rmin   =  x1b(is)
      rmax   =  x1b(is+ncells-1)

c--------  Call supernova ejecta profile calculator  ------------------
      call ejprof(Eej,Mej,n,rmax,rmin,vmax,ncells,rho0,v0)
c----------------------------------------------------------------------
c
c     assemble circumstellar medium
c
c--------  Call subroutine to calculate ambient medium structure  -----
c      call inter2D_2D(filename,ni,nk,dtemp,etemp,v1temp,v2temp)
c      call inter1D(filename,ni,dtemp,etemp,v1temp)
c----------------------------------------------------------------------
c
c      print*,"rho0 is: ",rho0
      tset         = rmax/vmax
c
c
        do  j=js,je
          do  i=is,ie
          if(dtemp(i,j).eq.0..or.dtemp(i,j).eq.tiny) 
     &           write(20,*)'Error in interpolation-d2',i,j,k
          if(etemp(i,j).eq.0..or.etemp(i,j).eq.tiny)
     &           write(20,*)'Error in interpolation-e2',i,j,k
          enddo
        enddo
c-----------------------------------------------------------------------
c
c
c     Set up atmosphere.
c
c
c      do 30 k=ks,ke
c        do 20 j=js,je
c          do 10 i=is,ie
c            d (i,j,k) = 0.01*rho0
c            v1(i,j,k) = 0.0D0
c            v2(i,j,k) = 0.0D0
c            v3(i,j,k) = 0.0D0
c            e (i,j,k) = (d(i,j,k)*boltz*amtemp) / 
c     &                    (mmw*mh*gamm1)
c            abun(i,j,k,1) = 1.
c            abun(i,j,k,2) = 0.
c10        continue
c20      continue
c30    continue
c
c
c     Set up atmosphere using an input file
c
c
c      
c      do  k=ks,ke
c        do  j=js,je
c          do  i=is,ie
c            d (i,j,k)     = dtemp(i,j)
c            v1(i,j,k)     = v1temp(i,j)
c            v3(i,j,k)     = 0.0D0
c            v2(i,j,k)     = v2temp(i,j)
c            e (i,j,k)     = etemp(i,j)
c
cc            d (i,j,k)     = dtemp(i)
cc            v1(i,j,k)     = v1temp(i)
cc            v2(i,j,k)     = 0.0D0
cc            v3(i,j,k)     = 0.0D0
cc            e (i,j,k)     = etemp(i)
cc
c            abun(i,j,k,1) = 1.
c            abun(i,j,k,2) = 0.
c          enddo
c        enddo
c      enddo

c
c
c     Set up atmosphere of a primordial halo
c
c
c      
      rhocrit  = 1.1314e-05 * h0**2 !(1.8788d-29 [g cm-3] * h^2 / mh)
      rho_igm  = omega_b * rhocrit * (1 + rdshft)**3 * ovrdns * mh

      do k=ks,ke
        do j=js,je
        do i=is-2,ie+2

          do l = 1,nbins
             if (x1b(i)/cmpc .ge. rtable(l  )  .and. 
     .           x1b(i)/cmpc .lt. rtable(l+1)) index = l
          enddo
          r = x1b(i)/cmpc

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
c            loggp    = gptable   (index)
            logr     = dlog10(r)
            r1       = -99.0
            r2       = dlog10(rtable(index))
            rdef     = r2 - r1
            loggp    = (logr-r1)*gptable(index)/rdef
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
            loggp    = gptable   (index)
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

            logHeI   = HeItable  (index  ) + (logr-r1)
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
            loggp    = gptable   (index  ) + (logr-r1)
     .               *(gptable   (index+1) - gptable   (index))/rdef
          endif

          d   (i,j,k  ) = 10.0**logd 
          e   (i,j,k  ) = 10.0**loge 
          v1  (i,j,k  ) = tiny !10.0**logv 
          tgas(i,j,k  ) = 10.0**logT
          abun(i,j,k,1) = 10.0**logHI
          abun(i,j,k,2) = 10.0**logHII
          abun(i,j,k,3) = 10.0**logel
          abun(i,j,k,4) = 10.0**logHeI
          abun(i,j,k,5) = 10.0**logHeII
          abun(i,j,k,6) = 10.0**logHeIII
          abun(i,j,k,7) = 10.0**logHM
          abun(i,j,k,8) = 10.0**logH2I
          abun(i,j,k,9) = 10.0**logH2II
          gp  (i,j,k  ) =-10.0**loggp 

c          if (d(i,j,k) .le. rho_igm) d(i,j,k) = rho_igm

        enddo
        enddo
      enddo
       
c
c
c     Set up supernova ejecta
c
c
      do k=ks,ke
        do j=js,je
          do i=is,is+ncells
            v1(i,j,k) = x1b(i)/tset
            v2(i,j,k) = tiny
            v3(i,j,k) = tiny
c
            if(v1(i,j,k).lt.v0) then
            d (i,j,k) = rho0
            else
            d (i,j,k) = rho0 * (v1(i,j,k)/v0)**(-n)
            endif
c
            abun(i,j,k,1) = 0.9999*fh 
            abun(i,j,k,2) = 0.0001*fh 
            abun(i,j,k,3) = abun(i,j,k,2)
            abun(i,j,k,4) = (1. - fh + tiny) 
            abun(i,j,k,5) = tiny
            abun(i,j,k,6) = tiny
            abun(i,j,k,7) = tiny  
            abun(i,j,k,8) = tiny ! 2.0d-06 * fh * 2.0  
            abun(i,j,k,9) = tiny

            tgas(i,j,k  ) = sntemp

c          mmwold=mmw
c          mmw=2.
c          e (i,j,k) = (d(i,j,k)*boltz*sntemp) / 
c     &                    (mmw*mh*gamm1)
c          mmw=mmwold
        e(i,j,k) =(abun(i,j,k,1)   + abun(i,j,k,2)    + abun(i,j,k,3)
     .           + abun(i,j,k,4)/4.+ abun(i,j,k,5)/4. + abun(i,j,k,6)/4.
     .           + abun(i,j,k,7)   + abun(i,j,k,8)/2. + abun(i,j,k,9)/2.
     .         ) * boltz * sntemp  * d(i,j,k)         / (mh * gamm1)

c            e (i,j,k) = (d(i,j,k)*boltz*sntemp) / 
c     &                    (mmw*mh*gamm1)
c            abun(i,j,k,1) = 0.
c            abun(i,j,k,2) = 1.
          enddo
        enddo
      enddo
       
      do  k=ks,ke
        do  j=js,je
          do  i=is,ie
          if(d(i,j,k).eq.0..or.d(i,j,k).eq.tiny) 
     &           write(20,*)'Error in interpolation-d',i,j,k
          if(e(i,j,k).eq.0..or.e(i,j,k).eq.tiny)
     &           write(20,*)'Error in interpolation-e',i,j,k
          enddo
        enddo
      enddo
c
c
c

      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             S U P E R N O V A             \\\\\\\\\\
c
c=======================================================================
c
