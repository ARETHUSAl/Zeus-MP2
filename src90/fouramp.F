c =======================================================================
c 
c     \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c     //////////              F O U R A M P                \\\\\\\\\\
c 
c =======================================================================
c
c    generates Fourier amplitudes for different k-modes;
c    have to only generate modes for half k-space;
c    choose no modes in minus z direction, if in x-y-plane no modes in
c    minus y direction, no mode for k=0 
C
C    modified 1: full preparation for fft now in A field
C                03/01 Robi Banerjee
C
C    modified 2: power spectrum with spectral index idx
C                09/01 Robi Banerjee
C
C    modified 3: reduce of 3D variables
C
C    INPUT VARIABLES:
C                
C                nmodes  number of excited modes
C                idx     spectral index
C                cc = 0: A*(k) = +A(-k)     
C                   = 1: A*(k) = -A(-k) (for a helical field)
C
      subroutine fouramp(ar,ai,nmodes,idx,cc)
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
      implicit NONE
      real(rl), intent(in)  :: idx
      integer,intent(in)    :: nmodes, cc
      real(rl), dimension(in, jn, kn), intent(inout) :: ar, ai
c
c LOCAL VARIABLES
c
      real(rl) :: kmag, kmin
      real(rl) :: help1 , help2
      integer  :: iseed = 1
      integer  :: i,j,k
c
c EXTERNALS
c
      real(rl)  ::    sigma
      external  sigma
c      real(rl)  ::   norm
c      external norm
c      A is vector potential 
c      i = 0: k=-N(pi/L), i=N: k=N(pi/L) Nyquist frequencies
c      exclude Nyquist frequencies

      if (nmodes.gt.ijkn/2) then
        write(*,*) 'FOURAMP: Two many modes !'
        stop
      endif

      kmin = 0.0d0

c initialize vector field
#ifdef MPI_USED
      if(myid_w.eq.0)
      write(*,fmt=123) 'FOURAMP: Excited k-mode range: kmin:',kmin 
     & ,'   kmax:',two*pi*sqrt(3.0)*real(nmodes)
#endif /* MPI_USED */
123   format(A,F9.2,A,F9.2)
      ar = 0.0
      ai = 0.0
c
c*nopar
      do 30 k=kn/2, kn/2+nmodes-1
        do 20 j=jn/2-nmodes+1, jn/2+nmodes-1         
          do 10 i=in/2-nmodes+1, in/2+nmodes-1 
            if (k.eq.kn/2) then    ! lies in xy-plane, make sure only half-space
              if (j.gt.jn/2) go to 1
              if (i.eq.in/2.and.i.gt.in/2) go to 1                 
              go to 4
            endif

1           continue


            kmag = two*pi*sqrt(float((i-in/2)**2)+          ! L - box size =1
     &               float((j-jn/2)**2)+float((i-in/2)**2))

            if ( kmag.ge.kmin ) then
c              help1 = sigma(kmag,idx)*norm(1)
c              help2 = two*pi*rand()
              call gasdev(iseed, help1)
              call ran1(iseed, help2)
              help1 = help1*sigma(kmag,idx)
              help2 = help2*two*pi
c             write(*,*) kmag,help1, help2
              ar(i,j,k) = help1*cos(help2)
              ai(i,j,k) = help1*sin(help2)

              if ( cc.eq.1 ) then
c
c the condition for A(x) to be real in this case
c is (A(k))* = -A(-k)
c
                ar(in-i,jn-j,kn-k) = - ar(i,j,k)
                ai(in-i,jn-j,kn-k) =   ai(i,j,k)

              else
                ar(in-i,jn-j,kn-k) =   ar(i,j,k)
                ai(in-i,jn-j,kn-k) = - ai(i,j,k)

              endif
3             continue          
 
            endif

4           continue

10        continue ! ix
20      continue ! iy
30    continue ! iz
c
c   zero mode
c
      ar(in/2,jn/2,kn/2) = 0.0
      ai(in/2,jn/2,kn/2) = 0.0
c      write(*,fmt=234) maxval(ar), minval(ar), maxval(ai), minval(ar)
c234   format(4(E10.3))
      return
      end
