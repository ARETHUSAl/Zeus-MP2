#include "rtchem.def"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c     //////////            S N R E S T A R T 2           \\\\\\\\\\
c
c
c           Routine to restart the supernova explosion problem
c
c
c           Developed by Dan Whalen
c             Last updated on: 10-22-2007
c=======================================================================
      subroutine snrestart2
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
      integer  :: i, j, k, sourceprob
      real(rl) :: Eej,Mej,n,vmax
      real(rl) :: sntemp, amtemp, vexpanmax, logrmin, logrmax, ovrdns

      real(rl)     :: dtemp(in,jn), etemp(in,jn), v1temp(in,jn)
      real(rl)     :: v2temp(in,jn)
      character*99 :: filename
      integer      :: ni, nk

      integer  :: l, index, nlines(12), nbins, nhalo       

      real(rl), dimension(200) :: dtable,Ttable,vrtable,rtable,
     .   etable,HItable,HIItable,eltable,HeItable,HeIItable,HeIIItable,
     .   HMtable,H2Itable,H2IItable,gptable

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
      Mej    =  8.0D33
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


      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////            S N R E S T A R T 2            \\\\\\\\\\
c
c=======================================================================
c
