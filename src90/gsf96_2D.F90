!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                 G S F 9 6                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
       subroutine gsf96
!
!
!  PURPOSE:  Sets up the the 2D perturbed density field for the I-front
!            enhanced dynamical instability examined by Garcia-Segura
!            and Franco, 1996 (ApJ 469: 171-188)
!
!  LOCAL VARIABLES:
!
!  EXTERNALS:  none
!
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use field
      use bndry
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
!
      implicit none
!
      integer     :: i, j
      real(rl)    :: rho, rc, rci, xi, n_central, r_core &
                     t_backg, r_trans, dd(in,jn,kn), y, q1 &
                     omega
!
      namelist / pgen     / omega, n_central, r_core, t_backg, r_trans
!
!-----------------------------------------------------------------------
!
      usrtag    =  'usr'
      omega     = -2.0
      n_central = 1.0e05
      r_core    = 2.1e16
      t_backg   = 100.0
      r_trans   = 3.084e17
      if (myid .eq. 0) then
        read  (1, pgen)
        write (2, pgen)
#ifdef MPI_USED
        buf_in(1) = omega 
        buf_in(2) = n_central
        buf_in(3) = r_core
        buf_in(4) = t_backg
        buf_in(5) = r_trans
      endif
      call MPI_BCAST( buf_in, 5, MPI_FLOAT &
                     , 0, comm3d, ierr )
      if (myid .ne. 0) then
        omega     = buf_in(1)
        n_central = buf_in(2)
        r_core    = buf_in(3)
        t_backg   = buf_in(4)
        r_trans   = buf_in(5)
#endif /* MPI_USED */
      endif
      rho = n_central * mh
      rc  = r_core
      rci = 1.0/r_core
      do j = js,je
      do i = is,ie
        xi = rc - x1a(i)
        q1 = dsign(0.5D0,xi)
        d(i,j,k) = rho*((0.5+q1)+(0.5-q1)*(x1a(i)*rci)**omega)
        e(i,j,k) = d(i,j,k) * boltz * t_backg / (mh * gamm1)
      enddo
      enddo
      call random_number(harvest=y)
      call random_number(dd)
      do j = js,je
      do i = is-1,ie-1
        if (x1b(i) .ge. r_trans) then
          d(i,j,k) = d(i,j,k) + 0.01 * d(i,j,k) * (dd(i,j,k) - 0.5)
        endif
      enddo
      enddo
! -- i faces
       nreq = 0
       nsub = nsub+1
       call bvald(3,3,0,0,0,0,d)
       call bvale(3,3,0,0,0,0,e)
       if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
! -- j faces
       nreq = 0
       nsub = nsub+1
       call bvald(0,0,3,3,0,0,d)
       call bvale(0,0,3,3,0,0,e)
       if(nreq .eq. 0) call mpi_waitall(nreq, req, stat, ierr)
       do j=js-3,je+3
       do i=is-1,ie-1
         gp(i+1,j,k) = gp(i,j,k) + gamm1 * 2.0 * &
                         (e(i+1,j,k) - e(i,j,k))/ &
                         (d(i,j,k) + d(i+1,j,k))
       enddo
       enddo
      return
      end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////                 G S F 9 6                 \\\\\\\\\\
!
!=======================================================================
!
