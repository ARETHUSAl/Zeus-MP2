!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               C G S O L V E               \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine cgsolve(dd, ddp1, &
                             ddp2, &
                             ddp3, &
                         x,rhs,toler,error)
!
      use real_prec
      use config
      use param
      use root
      use grid
      use impsoln
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez, &
                  ierror
!
      real(rl) ::   dd(neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
!
      real(rl) :: x   (neqm,in,jn,kn), &
                  rhs (neqm,in,jn,kn)
      real(rl) :: r   (neqm,in,jn,kn)
      real(rl) :: p   (neqm,in,jn,kn)
      real(rl) :: z   (neqm,in,jn,kn)
!
      real(rl) :: toler, error, glberror
!
      integer  :: i, jx, jy, jz
!
      integer  :: isxm1, isxp1, iexm1, iexp1, isym1, isyp1, ieym1, &
                  ieyp1, iszm1, iszp1, iezm1, iezp1
!
      real(rl) :: rnorm, beta_k, beta_ko, b_k
      real(rl) :: denom, alpha_k
!
      real(rl) :: ls_dprd
      real(rl) :: ls_dprd2
!
      error = 0.D0
!
!     set loop bounds
!
      isx = is
      iex = ie
      isy = js
      iey = je
      isz = ks
      iez = ke
      if(ldimen .eq. 1) then
       isym1 = isy
       isyp1 = isy
       ieym1 = isy
       ieyp1 = isy
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif 
      if(ldimen .eq. 2) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif
      if(ldimen .eq. 3) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz-1
       iszp1 = isz+1
       iezm1 = iez-1
       iezp1 = iez+1
      endif
!
!
!                          zero out the working arrays
!                          rbar, pbar, and zbar
!
      do i = 1, neqm
       do jz = isz-1, iez+1
        do jy = isy-1, iey+1
         do jx = isx-1, iex+1
          r   (i,jx,jy,jz) = 0.D0
          p   (i,jx,jy,jz) = 0.D0
          z   (i,jx,jy,jz) = 0.D0
         enddo
        enddo
       enddo
      enddo
!
!
!                                  set the iteration counter to zero
      nits=0
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!                                  multiply jacobian*x & store in r
      call sym_mul_bnd(isx,iex,isy,iey,isz,iez, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                      x,r)
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    r)
#endif
      call sym_mul_int(isx,iex,isy,iey,isz,iez, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                      x,r)
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!
!                                  now calculate residiuals
!
!     "x-faces":
!
      do jz = isz, iez, 1
       do jy = isy, iey, 1
        do i = 1,neqm,1
         r   (i,isx,jy,jz) = rhs(i,isx,jy,jz)-r(i,isx,jy,jz)
         r   (i,iex,jy,jz) = rhs(i,iex,jy,jz)-r(i,iex,jy,jz)
        enddo
       enddo
      enddo
      if(ldimen .gt. 1) then
!
!     "y-faces"
!
       do jz = isz, iez, 1
        do jx = isx+1, iex-1, 1
         do i = 1,neqm,1
          r   (i,jx,isy,jz) = rhs(i,jx,isy,jz)-r(i,jx,isy,jz)
          r   (i,jx,iey,jz) = rhs(i,jx,iey,jz)-r(i,jx,iey,jz)
         enddo
        enddo
       enddo
       if(ldimen .eq. 3) then
!
!     "z-faces":
!
        do jy = isyp1, ieym1, 1
         do jx = isx+1, iex-1, 1
          do i = 1,neqm,1
           r   (i,jx,jy,isz) = rhs(i,jx,jy,isz)-r(i,jx,jy,isz)
           r   (i,jx,jy,iez) = rhs(i,jx,jy,iez)-r(i,jx,jy,iez)
          enddo
         enddo
        enddo
       endif ! ldimen = 3
      endif ! ldimen > 1
!
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    r   )
#endif
!
      do jz = iszp1,iezm1,1
       do jy = isyp1,ieym1,1
        do jx = isx+1,iex-1,1
         do i = 1,neqm,1
            r   (i,jx,jy,jz) = rhs(i,jx,jy,jz)-r(i,jx,jy,jz)
         enddo
        enddo
       enddo
      enddo
!
!
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!
!                                 calculate norm of rhs
      rnorm = ls_dprd(isx,iex,isy,iey,isz,iez, &
                      rhs,rhs)
      rnorm = sqrt(rnorm)
!
!                                 precondition r & store in z
      call sym_prcn_bnd(isx,iex,isy,iey,isz,iez, &
                       ipcflag,0, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                      z,r)
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    z)
#endif
      call sym_prcn_int(isx,iex,isy,iey,isz,iez, &
                       ipcflag,0, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                      z,r)
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!     >------------------------------------------------------------
!                                 main loop
10    if(nits.le.maxitr) then
!
!                                 increment interation number
        nits=nits+1
!
!                                 calculate beta_k
      beta_k = ls_dprd(isx,iex,isy,iey,isz,iez, &
                       z,r)
!
!                                 if this is the first iteration
!                                 then set p and pbar to z and zbar
!                                 otherwise, calculate p & pbar by
!                                 the recursion relationship
!
!     "x-faces"
!
      if(nits.eq.1) then
       do jz = isz, iez, 1
        do jy = isy, iey, 1
         do i = 1,neqm,1
          p   (i,isx,jy,jz) = z   (i,isx,jy,jz)
          p   (i,iex,jy,jz) = z   (i,iex,jy,jz)
         enddo
        enddo
       enddo
      else
       b_k = beta_k/(beta_ko+tiny)
       do jz = isz, iez, 1
        do jy = isy, iey, 1
         do i = 1,neqm,1
          p   (i,isx,jy,jz) = b_k*p   (i,isx,jy,jz) + z   (i,isx,jy,jz)
          p   (i,iex,jy,jz) = b_k*p   (i,iex,jy,jz) + z   (i,iex,jy,jz)
         enddo
        enddo
       enddo
      endif
      if(ldimen .gt. 1)then
!
!     "y-faces"
!
       if(nits.eq.1) then
        do jz = isz, iez, 1
         do jx = isx+1, iex-1, 1
          do i = 1,neqm,1
           p   (i,jx,isy,jz) = z   (i,jx,isy,jz)
           p   (i,jx,iey,jz) = z   (i,jx,iey,jz)
          enddo
         enddo
        enddo
       else
        b_k = beta_k/(beta_ko+tiny)
        do jz = isz, iez, 1
         do jx = isx+1, iex-1, 1
          do i = 1,neqm,1
           p(i,jx,isy,jz) = b_k*p(i,jx,isy,jz) + z(i,jx,isy,jz)
           p(i,jx,iey,jz) = b_k*p(i,jx,iey,jz) + z(i,jx,iey,jz)
          enddo
         enddo
        enddo
       endif
       if(ldimen .gt. 2) then
!
!     "z-faces"
!
        if(nits.eq.1) then
         do jy = isyp1, ieym1, 1
          do jx = isx+1, iex-1, 1
           do i = 1,neqm,1
            p   (i,jx,jy,isz) = z   (i,jx,jy,isz)
            p   (i,jx,jy,iez) = z   (i,jx,jy,iez)
           enddo
          enddo
         enddo
        else
         b_k = beta_k/(beta_ko+tiny)
         do jy = isyp1, ieym1, 1
          do jx = isx+1, iex-1, 1
           do i = 1,neqm,1
            p(i,jx,jy,isz) = b_k*p(i,jx,jy,isz) + z(i,jx,jy,isz)
            p(i,jx,jy,iez) = b_k*p(i,jx,jy,iez) + z(i,jx,jy,iez)
           enddo
          enddo
         enddo
        endif
       endif ! ldimen = 3
      endif ! ldimen > 1
!
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    p   )
#endif
!
!     "interior"
!
      if(nits .eq. 1) then
       do jz = iszp1,iezm1,1
        do jy = isyp1,ieym1,1
         do jx = isx+1,iex-1,1
          do i=1,neqm,1
             p   (i,jx,jy,jz) = z   (i,jx,jy,jz)
          enddo
         enddo
        enddo
       enddo
      else
       b_k = beta_k/(beta_ko+tiny)
       do jz = iszp1,iezm1,1
        do jy = isyp1,ieym1,1
         do jx = isx+1,iex-1,1
          do i = 1,neqm,1
             p   (i,jx,jy,jz) = b_k*p   (i,jx,jy,jz)+z   (i,jx,jy,jz)
          enddo
         enddo
        enddo
       enddo
      endif
!
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!
!                                 save beta_k for use on next pass
      beta_ko = beta_k
!
!                                 multiply jacobian*p & store in z
      call sym_mul_bnd(isx,iex,isy,iey,isz,iez, &
                      dd, ddp1, &
                          ddp2, &
                          ddp3, &
                      p,z)
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    z)
#endif
      call sym_mul_int(isx,iex,isy,iey,isz,iez, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                      p,z)
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!                                 now calculate denominator of a_k
      denom = ls_dprd(isx,iex,isy,iey,isz,iez, &
                      z,p)
!
!                                 calculate alpha_k
      alpha_k = beta_k/(denom+tiny)
!
!
!                                 calculate solution vector &
!                                 update residuals
!
!     "x-faces"
!
      do jz = isz,iez,1
       do jy = isy,iey,1
        do i = 1,neqm,1
         x   (i,isx,jy,jz) = x   (i,isx,jy,jz) + alpha_k* &
                             p   (i,isx,jy,jz)
         r   (i,isx,jy,jz) = r   (i,isx,jy,jz) - alpha_k* &
                             z   (i,isx,jy,jz)
         x   (i,iex,jy,jz) = x   (i,iex,jy,jz) + alpha_k* &
                             p   (i,iex,jy,jz)
         r   (i,iex,jy,jz) = r   (i,iex,jy,jz) - alpha_k* &
                             z   (i,iex,jy,jz)
        enddo
       enddo
      enddo
      if(ldimen .gt. 1) then
!
!     "y-faces"
!
       do jz = isz,iez,1
        do jx = isx+1,iex-1,1
         do i=1,neqm,1
          x   (i,jx,isy,jz) = x   (i,jx,isy,jz) + alpha_k* &
                              p   (i,jx,isy,jz)
          r   (i,jx,isy,jz) = r   (i,jx,isy,jz) - alpha_k* &
                              z   (i,jx,isy,jz)
          x   (i,jx,iey,jz) = x   (i,jx,iey,jz) + alpha_k* &
                              p   (i,jx,iey,jz)
          r   (i,jx,iey,jz) = r   (i,jx,iey,jz) - alpha_k* &
                              z   (i,jx,iey,jz)
         enddo
        enddo
       enddo
       if(ldimen .gt. 2) then
!
!     "z-faces"
!
        do jy = isyp1,ieym1,1
         do jx = isx+1,iex-1,1
          do i = 1,neqm,1
           x   (i,jx,jy,isz) = x   (i,jx,jy,isz) + alpha_k* &
                               p   (i,jx,jy,isz)
           r   (i,jx,jy,isz) = r   (i,jx,jy,isz) - alpha_k* &
                               z   (i,jx,jy,isz)
           x   (i,jx,jy,iez) = x   (i,jx,jy,iez) + alpha_k* &
                               p   (i,jx,jy,iez)
           r   (i,jx,jy,iez) = r   (i,jx,jy,iez) - alpha_k* &
                               z   (i,jx,jy,iez)
          enddo
         enddo
        enddo
       endif ! ldimen > 2
      endif ! ldimen > 1
!
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    x   )
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    r   )
#endif
!
      do jz = iszp1,iezm1,1
       do jy = isyp1,ieym1,1
        do jx = isx+1,iex-1,1
         do i=1,neqm,1
            x   (i,jx,jy,jz) = x   (i,jx,jy,jz)+alpha_k* &
                               p   (i,jx,jy,jz)
            r   (i,jx,jy,jz) = r   (i,jx,jy,jz)-alpha_k* &
                               z   (i,jx,jy,jz)
         enddo
        enddo
       enddo
      enddo
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!
!                                 solve for new z based on
!                                 updated residual
      call sym_prcn_bnd(isx,iex,isy,iey,isz,iez, &
                       ipcflag,0, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                       z,r)
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call sendrec_bnd( &
                    isx,iex,isy,iey,isz,iez, &
                    z)
#endif
      call sym_prcn_int(isx,iex,isy,iey,isz,iez, &
                       ipcflag,0, &
                       dd, ddp1, &
                           ddp2, &
                           ddp3, &
                       z,r)
!
#ifdef MPI_USED
      if(nreq .ne. 0) call mpi_waitall(nreq, req, stat, ierr)
#endif
!
!                                 estimate the error
!
      if(mod(nits,2) .ne. 0) go to 10
!
      if(cgerrcrit .eq. 1) then
       error = ls_dprd(isx,iex,isy,iey,isz,iez,r,r)
       error = sqrt(error)/(rnorm+tiny)
      else if(cgerrcrit .eq. 2) then
       error = 0.D0
       do i = 1, neqm
        do jz = isz, iez
         do jy = isy, iey
          do jx = isx, iex
           error = max(error,abs(r(i,jx,jy,jz)/(rhs(1,jx,jy,jz)+tiny)))
          enddo
         enddo
        enddo
       enddo
      else
       print *,'CGSOLVE: cgerrcrit has invalid value: ',cgerrcrit
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif /* MPI_USED */
       stop
      endif
#ifdef MPI_USED
      call MPI_ALLREDUCE(error,glberror,1,MPI_FLOAT,mpi_max, &
                         comm3d,ierr)
#else
      glberror = error
#endif /* MPI_USED */
!
!                                 if error  is still too large, then
!                                 iterate again
      if(glberror.gt.toler) goto 10
!
      endif ! 10
!
!
!            ok, we're done so return to calling routine
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               C G S O L V E               \\\\\\\\\\
!
!=======================================================================
