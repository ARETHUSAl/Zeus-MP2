c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              G E N V E L                  \\\\\\\\\\
c
c=======================================================================
c
       subroutine genvel(nmodes,vrms,idx)
C*********************************************************************
C** DON'T USE INCLUDE STATEMENT WITHIN ZEUS-EDITOR:
C** IT GETS CONFUSED WITH THIS "UNKNOWN" STATEMENT 
C** AND TAKE CARE OF "SYNTAX", I.E. USE SUFFICIENT
C** SPACE BETWEEN DECLARATION AND VARIABLES
C*********************************************************************
c
c
c    written by: Robi Banerjee             09/2001
c    modified 1: ATTENTION                                  !!!!! 
c                genvel must be called before generation    !!!!!
c                of magnetic fields as genvel uses          !!!!!
c                b1,b2 and b3 as help variables             !!!!!
c                Robi Banerjee             01/2002
c
c  PURPOSE: generate stochastic velocity field
c
c           if VELDIVFREE is set the initial velocity field
c           is divergence free
c
c
c  INPUT VARIABLES: nmodes     number of excited modes
c                   vrms  rms velocity
c                   idx   spectral index
c
c  OUTPUT VARIABLES: 
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use vcvars
      implicit NONE   
      real(rl)   ::   ran1
*     external  bperiodic, fouramp, fafotr, normvel, ran1
c     1        , divergence, variance

*
      integer, intent(in)   :: nmodes
      real(rl),intent(in)   :: vrms, idx
c      integer               :: N
      integer, dimension(6) :: st
c
c local variables
c
      real(rl), dimension(3) :: var
      integer :: iseed, i,j,k, ip, jp, kp
c
c allocate local memory (instead of equivalencing)
c
      allocate(vp1r(in, jn, kn), stat=st(1))
      allocate(vp1i(in, jn, kn), stat=st(2))
      allocate(vp2r(in, jn, kn), stat=st(3))
      allocate(vp2i(in, jn, kn), stat=st(4))
      allocate(vp3r(in, jn, kn), stat=st(5))
      allocate(vp3i(in, jn, kn), stat=st(6))
c
c EXTERNALS
c
      iseed = 500

#ifdef MPI_USED
      if (myid_w .eq. 0) then       
        write(*,*) 'GENVEL                 : ------------------ '
        write(*,*) 'GENVEL                 : start generate of  '
        write(*,*) 'GENVEL                 : velocity field     '
        write(*,*) 'GENVEL                 : ------------------ '
        write(*,*) 'GENVEL                 : ', ran1(-iseed)
        write(*,*) 'GENVEL                 : ------------------ '
      endif
#endif /* MPI_USED */

      call fouramp(vp1r, vp1i, nmodes,idx,0) 
      call fouramp(vp2r, vp2i, nmodes,idx,0) 
      call fouramp(vp3r, vp3i, nmodes,idx,0) 

      call fafotr(vp1r, vp1i)
      call fafotr(vp2r, vp2i)
      call fafotr(vp3r, vp3i)
c
c -------------------------------------------------------------------------
c ------- move vector potential to edge centers ---------------------------
c -------------------------------------------------------------------------
      do 330 k=1,kn
        if (k .ne. kn) then
          kp = k + 1 
        else 
          kp = 1
        endif
        do 320 j=1,jn
          if (j .ne. jn) then
            jp = j + 1 
          else 
            jp = 1
          endif
          do 310 i=1,in
            if (i .ne. in) then
              ip = i + 1 
            else 
              ip = 1
            endif

            w3da(i,j,k) = (vp1r(i ,j ,k ) + vp1r(ip,j ,k ) )*0.5d0
            w3db(i,j,k) = (vp2r(i ,j ,k ) + vp2r(i ,jp,k ) )*0.5d0
            w3dc(i,j,k) = (vp3r(i ,j ,k ) + vp3r(i ,j ,kp) )*0.5d0

310       continue
320     continue
330   continue

c   Used the loop from this version of Zeus, mot the one below
      do 430 k=ks-2,ke+2
        do 420 j=js-2,je+2
          do 410 i=is-2,ie+2
            vp1i(i,j,k) = w3da(i,j,k)
            vp2i(i,j,k) = w3db(i,j,k)
            vp3i(i,j,k) = w3dc(i,j,k)
410        continue
420      continue
430    continue
c
c  update all the boundary values 
c
c      nreq=0
c      nsub = nsub + 1
c      call bvalv1(3,3,0,0,0,0,vp1i)
c      call bvalv2(3,3,0,0,0,0,vp2i)
c      call bvalv3(3,3,0,0,0,0,vp3i)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalv1(0,0,3,3,0,0,vp1i)
c      call bvalv2(0,0,3,3,0,0,vp2i)
c      call bvalv3(0,0,3,3,0,0,vp3i)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalv1(0,0,0,0,3,3,vp1i)
c      call bvalv2(0,0,0,0,3,3,vp2i)
c      call bvalv3(0,0,0,0,3,3,vp3i)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

c      call bperiodic( vp1i )
c      call bperiodic( vp2i )
c      call bperiodic( vp3i )

#ifdef VELDIVFREE
c
c compute divergence free velocity field out of vector potential
c

      do 130 k=ks,ke
        kp = k + 1
        do 120 j=js,je
          jp = j + 1
          do 110 i=is,ie
            ip = i + 1

            v1(i,j,k) = 
     &                   ( vp3i(i,jp,k) - vp3i(i,j,k) -
     &                     vp2i(i,j,kp) + vp2i(i,j,k) )*real(N)

            v2(i,j,k) = 
     &                   ( vp1i(i,j,kp) - vp1i(i,j,k) -
     &                     vp3i(ip,j,k) + vp3i(i,j,k) )*real(N)

            v3(i,j,k) = 
     &                   ( vp2i(ip,j,k) - vp2i(i,j,k) -
     &                     vp1i(i,jp,k) + vp1i(i,j,k) )*real(N)

110       continue
120     continue
130   continue
#else /* VELDIVFREE */
      do 130 k=1,kn
        do 120 j=1,jn
          do 110 i=1,in
            v1(i,j,k) = vp1i(i,j,k)
            v2(i,j,k) = vp2i(i,j,k)
            v3(i,j,k) = vp3i(i,j,k)
110       continue
120     continue
130   continue
#endif /* VELDIVFREE */

      deallocate(vp1i)
      deallocate(vp1r)
      deallocate(vp2i)
      deallocate(vp2r)
      deallocate(vp3i)
      deallocate(vp3r)
c
c
c  update all the boundary values of v 
c
c     nreq=0
c     nsub = nsub + 1
c     call bvalv1(3,3,0,0,0,0,v1)
c     call bvalv2(3,3,0,0,0,0,v2)
c     call bvalv3(3,3,0,0,0,0,v3)
c
c     if (nreq .ne. 0) call MPI_WAITALL(nreq, req, stat, ierr)
c     nreq=0
c     nsub = nsub + 1
c     call bvalv1(0,0,3,3,0,0,v1)
c     call bvalv2(0,0,3,3,0,0,v2)
c     call bvalv3(0,0,3,3,0,0,v3)
c
c     if (nreq .ne. 0) call MPI_WAITALL(nreq, req, stat, ierr)
c     nreq=0
c     nsub = nsub + 1
c     call bvalv1(0,0,0,0,3,3,v1)
c     call bvalv2(0,0,0,0,3,3,v2)
c     call bvalv3(0,0,0,0,3,3,v3)
c
c     if (nreq .ne. 0) call MPI_WAITALL(nreq, req, stat, ierr)
c
c   Normalize velocity field to vrms value
c       
      call normvel(vrms)

c   Diagnostic output to screen
c       
      write(*,*) 'GENVEL                 : ------------------ '
      write(*,*) 'GENVEL                 : generation of '
      write(*,*) 'GENVEL                 : velocity field '
      write(*,*) 'GENVEL                 : finished !! '
      write(*,*) 'GENVEL                 : ------------------ '

      end
