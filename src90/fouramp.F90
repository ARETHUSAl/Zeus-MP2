! =======================================================================
! 
!     \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!     //////////              F O U R A M P                \\\\\\\\\\
! 
! =======================================================================
!
!    generates Fourier amplitudes for different k-modes;
!    have to only generate modes for half k-space;
!    choose no modes in minus z direction, if in x-y-plane no modes in
!    minus y direction, no mode for k=0 
!
!    modified 1: full preparation for fft now in A field
!                03/01 Robi Banerjee
!
!    modified 2: power spectrum with spectral index idx
!                09/01 Robi Banerjee
!
!    modified 3: reduce of 3D variables
!
!    INPUT VARIABLES:
!                
!                nmodes  number of excited modes
!                idx     spectral index
!                cc = 0: A*(k) = +A(-k)     
!                   = 1: A*(k) = -A(-k) (for a helical field)
!
subroutine fouramp(ar,ai,nmodes,idx,cc)
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
!
! LOCAL VARIABLES
!
      real(rl) :: kmag, kmin
      real(rl) :: help1 , help2
      integer,save  :: iseed = 192
      integer  :: i,j,k
!
! EXTERNALS
!
      real(rl)  ::    sigma
      external  sigma
!      real(rl)  ::   norm
!      external norm
!      A is vector potential 
!      i = 0: k=-N(pi/L), i=N: k=N(pi/L) Nyquist frequencies
!      exclude Nyquist frequencies
      if (nmodes.gt.ijkn/2) then
        write(*,*) 'FOURAMP: Two many modes !'
        stop
      endif
      kmin = 0.0d0
! initialize vector field
#ifdef MPI_USED
      if(myid_w.eq.0) then
        write(*,fmt=123) 'FOURAMP: Excited k-mode range: kmin:',kmin &
         ,'   kmax:',two*pi*sqrt(3.0)*real(nmodes)
      endif
#endif /* MPI_USED */
123   format(A,F9.2,A,F9.2)
      ar = 0.0
      ai = 0.0
!
!*nopar
      do k=kn/2, kn/2+nmodes-1
        do j=jn/2-nmodes+1, jn/2+nmodes-1         
          do i=in/2-nmodes+1, in/2+nmodes-1 
            if (k.eq.kn/2) then    ! lies in xy-plane, make sure only half-space
              if (j.gt.jn/2) go to 1
              if (i.eq.in/2.and.i.gt.in/2) go to 1                 
              go to 4
            endif
1           continue
            kmag = two*pi*sqrt(float((i-in/2)**2)+  &      ! L - box size =1
                     float((j-jn/2)**2)+float((i-in/2)**2))
            if ( kmag.ge.kmin ) then
!              help1 = sigma(kmag,idx)*norm(1)
!              help2 = two*pi*rand()
              call gasdev(iseed, help1)
              call ran1(iseed, help2)
              help1 = help1*sigma(kmag,idx)
              help2 = help2*two*pi
!             write(*,*) kmag,help1, help2
              ar(i,j,k) = help1*cos(help2)
              ai(i,j,k) = help1*sin(help2)
              if ( cc.eq.1 ) then
!
! the condition for A(x) to be real in this case
! is (A(k))* = -A(-k)
!
                ar(in-i,jn-j,kn-k) = - ar(i,j,k)
                ai(in-i,jn-j,kn-k) =   ai(i,j,k)
              else
                ar(in-i,jn-j,kn-k) =   ar(i,j,k)
                ai(in-i,jn-j,kn-k) = - ai(i,j,k)
              endif
3             continue          
            endif
4           continue
          enddo ! i
        enddo ! j
      enddo ! k
!
!   zero mode
!
      ar(in/2,jn/2,kn/2) = 0.0
      ai(in/2,jn/2,kn/2) = 0.0
!      write(*,fmt=234) maxval(ar), minval(ar), maxval(ai), minval(ar)
!234   format(4(E10.3))
      return
end subroutine fouramp
