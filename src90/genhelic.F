c  generates helical magnetic field
c =======================================================================
c 
c     \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c     //////////             G E N H E L I C               \\\\\\\\\\
c 
c =======================================================================
c 
      subroutine genhelic(nmodes, Brms, idx, helic)
c
      use IFPORT
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
c
      implicit NONE
c
c    rb:zeus3d.helic <----------------------------------------------
c                                                        july, 2001
c
c    written by: Robi Banerjee
c    modified 1: ?
c
c  PURPOSE:          generate vector potential with explicit
c                    helicity
c                    We use the coordinate base {e_+,e_-,k}
c                    to generate the vector potential
c
c  MODIFIED:         reduce of 3D variables        Robi Banerjee,  jan, 2002
c
c
c  INPUT VARIABLES: 
c                   nmodes  :  # of modes to excite
c                   Brms    : rms value of magnetic field
c                   idx     : spectral index
c                   helic   : percentage (0.0 - 1.0) of the maximal helicity
c
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
      real(rl)    ::  normmag
      external ::  bvalb1,bvalb2,bvalb3, fouramp, fafotr, normmag
c 
c input variables
c 
      real(rl), intent(in) ::  Brms, idx, helic
      integer, intent(in)  ::  nmodes 
c  
c EXTERNALS
c 
c 
c local variables
c
      real(rl), dimension(:,:,:), allocatable, save  :: a1, a2, a3
      real(rl), dimension(:,:,:), allocatable, save  :: ah1i, ah2i, ah3i
      real(rl) :: k1, k2, k3 
      real(rl) :: ep1r, ep1i, ep2r, ep2i, ep3r, ep3i
      real(rl) :: em1r, em1i, em2r, em2i, em3r, em3i
      real(rl) :: kmag, en, e11, e12, e13, e21, e22, e23
      real(rl) :: A1r, A1i, A2r, A2i, A3r, A3i
      real(rl) :: var(3), sr2i
      real(rl) :: helicity, norm, frac

      integer ::  iseed
      integer ::  i, j, k
      integer ::  ip, jp, kp
      integer, dimension(6) :: st
 
      allocate(a1(in,jn,kn), stat=st(1))
      allocate(a2(in,jn,kn), stat=st(2))
      allocate(a3(in,jn,kn), stat=st(3))
      allocate(ah1i(in,jn,kn), stat=st(4))
      allocate(ah2i(in,jn,kn), stat=st(5))
      allocate(ah3i(in,jn,kn), stat=st(6))

#ifdef MPI_USED
      if (myid_w .eq. 0) then
        write(6,2030) '----------------------------------------------'
        write(6,2030) 'start generation of magnetic field ... '
        write(6,2030) '----------------------------------------------'
      endif
#endif */ MPI_USED */

      iseed = 200
      sr2i = 1.d0/sqrt(2.d0)
      frac = sqrt((1.d0-helic)/(1.d0+helic))

#ifdef MPI_USED
      if (myid_w .eq. 0) then
        write(6,2040) 'frac = ',frac
        write(6,2040) 'start random number generator'
        call srand(-iseed)
      endif
#endif */ MPI_USED */

c--------------------------------------------------------------------------
c --------------- excite modes --------------------------------------------
c -------------------------------------------------------------------------
      call fouramp(a1,ah1i,nmodes,idx,1) 
      call fouramp(a2,ah2i,nmodes,idx,1) 
      call fouramp(a3,ah3i,nmodes,idx,1) 
c
c generate unit vectors
c origin of coordinate system is at (N/2,N/2,N/2)
c transform amplitudes to new basis
c
      do k=1,kn
        do j=1,jn
          do i=1,in
            kmag = sqrt(real(i-in/2)**2 
     1                + real(j-jn/2)**2 + real(k-kn/2)**2)
            if ( kmag .gt. 0.0 ) then
              k1 = real(i-in/2)/kmag
              k2 = real(j-jn/2)/kmag
              k3 = real(k-kn/2)/kmag
            else
              k1 = 0.0
              k2 = 0.0
              k3 = 0.0
            endif

            e11 = k2-k3
            e12 = k3-k1
            e13 = k1-k2
            en = sqrt(e11**2 + e12**2 + e13**2)
            if ( en > 0.0 ) then 
              e11 = e11/en
              e12 = e12/en
              e13 = e13/en
            endif

            e21 = k2*e13 - k3*e12
            e22 = k3*e11 - k1*e13
            e23 = k1*e12 - k2*e11
            en = sqrt(e21**2 + e22**2 + e23**2)
            if ( en > 0.0) then
              e21 = e21/en
              e22 = e22/en
              e23 = e23/en
            endif

            ep1r = sr2i*e11
            ep2r = sr2i*e12
            ep3r = sr2i*e13

            ep1i = sr2i*e21
            ep2i = sr2i*e22
            ep3i = sr2i*e23

            em1r = ep1r
            em2r = ep2r
            em3r = ep3r

            em1i = -ep1i
            em2i = -ep2i
            em3i = -ep3i

            A1r =   a1(i,j,k)
            A1i = ah1i(i,j,k)
      
            A2r = frac * A1r
            A2i = frac * A1i

            A3r = 0.0
            A3i = 0.0
c
c   1-vector potential
c
            a1(i,j,k)    =   A1r*ep1r  
     &                     - A1i*ep1i  
     &                     + A2r*em1r  
     &                     - A2i*em1i  
     &                     + A3r*k1
            ah1i(i,j,k)  =   A1i*ep1r  
     &                     + A1r*ep1i  
     &                     + A2i*em1r  
     &                     + A2r*em1i  
     &                     + A3i*k1
c 
c   2-vector potential
c
            a2(i,j,k)    =   A1r*ep2r 
     &                     - A1i*ep2i 
     &                     + A2r*em2r 
     &                     - A2i*em2i 
     &                     + A3r*k2 
            ah2i(i,j,k)  =   A1i*ep2r 
     &                     + A1r*ep2i 
     &                     + A2i*em2r  
     &                     + A2r*em2i 
     &                     + A3i*k2
c 
c   3-vector potential
c
            a3(i,j,k)    =   A1r*ep3r 
     &                     - A1i*ep3i 
     &                     + A2r*em3r  
     &                     - A2i*em3i 
     &                     + A3r*k3
            ah3i(i,j,k)  =   A1i*ep3r 
     &                     + A1r*ep3i 
     &                     + A2i*em3r 
     &                     + A2r*em3i 
     &                     + A3i*k3

           enddo
         enddo
       enddo

c -------------------------------------------------------------------------
c ------- find real space variables by fft --------------------------------
c -------------------------------------------------------------------------
      call fafotr(a1, ah1i)
      call fafotr(a2, ah2i)
      call fafotr(a3, ah3i)
c -------------------------------------------------------------------------
c ------- move A to edge centers ------------------------------------------
c -------------------------------------------------------------------------
      do k=1,kn
        if (k .NE. kn) then
          kp = k + 1 
        else 
          kp = 1
        endif
        do j=1,jn
          if (j .NE. jn) then
            jp = j + 1 
          else 
            jp = 1
          endif
          do i=1,in
            if (i .NE. in) then
              ip = i + 1 
             else 
              ip = 1
             endif

            w3dd(i,j,k) = ( a1(i ,j ,k ) + a1(ip,j ,k ) )*0.5d0
            w3de(i,j,k) = ( a2(i ,j ,k ) + a2(i ,jp,k ) )*0.5d0
            w3df(i,j,k) = ( a3(i ,j ,k ) + a3(i ,j ,kp) )*0.5d0

          enddo   
        enddo 
      enddo   

      do 430 k=ks-2,ke+2
        do 420 j=js-2,je+2
          do 410 i=is-2,ie+2
            a1(i,j,k) = w3dd(i,j,k)
            a2(i,j,k) = w3de(i,j,k)
            a3(i,j,k) = w3df(i,j,k)
410       enddo
420     enddo
430   enddo

c
c  update all the boundary values of a
c
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(3,3,0,0,0,0,a1)
c      call bvalb2(3,3,0,0,0,0,a2)
c      call bvalb3(3,3,0,0,0,0,a3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(0,0,3,3,0,0,a1)
c      call bvalb2(0,0,3,3,0,0,a2)
c      call bvalb3(0,0,3,3,0,0,a3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(0,0,0,0,3,3,a1)
c      call bvalb2(0,0,0,0,3,3,a2)
c      call bvalb3(0,0,0,0,3,3,a3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

c ------------------------------------------------------------------------
c ---- compute b out of curl A -------------------------------------------
c ------------------------------------------------------------------------
      do k=ks,ke
        kp = k + 1
        do j=js,je
          jp = j + 1
          do i=is,ie
            ip = i + 1

            b1(i,j,k) = 
     &                   ( a3(i,jp,k) - a3(i,j,k) -
     &                     a2(i,j,kp) + a2(i,j,k) )*real(ijkn)

            b2(i,j,k) = 
     &                   ( a1(i,j,kp) - a1(i,j,k) -
     &                     a3(ip,j,k) + a3(i,j,k) )*real(ijkn)

            b3(i,j,k) = 
     &                   ( a2(ip,j,k) - a2(i,j,k) -
     &                     a1(i,jp,k) + a1(i,j,k) )*real(ijkn)

          enddo
        enddo
      enddo

c
c  update all the boundary values of b 
c
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(3,3,0,0,0,0,b1)
c      call bvalb2(3,3,0,0,0,0,b2)
c      call bvalb3(3,3,0,0,0,0,b3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(0,0,3,3,0,0,b1)
c      call bvalb2(0,0,3,3,0,0,b2)
c      call bvalb3(0,0,3,3,0,0,b3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
c      nreq=0
c      nsub = nsub + 1
c      call bvalb1(0,0,0,0,3,3,b1)
c      call bvalb2(0,0,0,0,3,3,b2)
c      call bvalb3(0,0,0,0,3,3,b3)
c
c      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

      norm = normmag(Brms)

      a1 = a1/norm
      a2 = a2/norm
      a3 = a3/norm
 
      helicity = 0.0

      do k=ks,ke
        kp = k+1
        do j=js,je
          jp = j+1
          do i=is,ie
            ip = i+1
            helicity = helicity + 
     &                ( ( a1(i ,j ,k ) + a1(i ,j ,kp) +
     &                    a1(i ,jp,k ) + a1(i ,jp,kp) ) * 
     &                  ( b1(i ,j ,k ) + b1(ip,j ,k ) ) +
     &                  ( a2(i ,j ,k ) + a2(i ,j ,kp) +
     &                    a2(ip,j ,k ) + a2(ip,j ,kp) ) *
     &                  ( b2(i ,j ,k ) + b2(i ,jp,k ) ) +
     &                  ( a3(i ,j ,k ) + a3(i ,jp,k ) +
     &                    a3(ip,j ,k ) + a3(ip,jp,k ) ) * 
     &                  ( b3(i ,j ,k ) + b3(i ,j ,kp) ) ) * 
     &                dvl1a(i) * dvl2a(j) * dvl3a(k) * 0.125
          enddo
        enddo
      enddo

      deallocate(a1, stat=st(1))
      deallocate(a2, stat=st(2))
      deallocate(a3, stat=st(3))
      deallocate(ah1i, stat=st(4))
      deallocate(ah2i, stat=st(5))
      deallocate(ah3i, stat=st(6))

#ifdef MPI_USED
      if (myid_w .eq. 0) then      
        write(6,2030) '----------------------------------------------'
        write(6,2040) 'helicity    :', helicity
        write(6,2030) 'generation of magnetic field finished  '
        write(6,2030) '----------------------------------------------'
      endif
#endif */ MPI_USED */

2030  format('GENHELIC : ',a50)
2040  format('GENHELIC : ',a25,1pe13.6)
2050  format('GENHELIC : ',a25,i7)

      end subroutine genhelic
