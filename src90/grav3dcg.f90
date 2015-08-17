!=======================================================================
!
!                            Developed by
!                Laboratory for Computational Astrophysics
!                  University of California at San Diego
!
      subroutine grav3D_CG
!
!     3-D gravitational potential solver using John Hayes's
!     CGSOLVE package for the covariant form of Poisson's equation
!
      use real_prec
      use config
      use param
      use cons
      use domain
      use root
      use field
      use grid
      use bndry
      use impsoln
      use gravmod
      use mpiyes
      use mpipar
!
      implicit NONE
!
      real(rl), dimension(:,:,:,:,:), allocatable, save :: &
                ad, a0, a1, a2
      real(rl), dimension(:,:,:,:), allocatable, save :: &
                x, b
!
      real(rl), dimension(:,:), allocatable :: &
                gpiibc, gpoibc, gpijbc, gpojbc, gpikbc, gpokbc
!
      real(rl) :: toler, error, glbtoler, eps, glbeps, epstol, epsmax, &
                  epsmin, glbemin, phip(in,jn,kn)
!
      integer :: i, j, k, prev_tot, imax, jmax, kmax
!
!-----------------------------------------------------------------------
!     First timestep only: initialize everything required for CGSOLVE
!-----------------------------------------------------------------------
!
      if(igcall .eq. 0) then
!
!-----------------------------------------------------------------------
!      allocate matrix, solution, and RHS array storage
!-----------------------------------------------------------------------
!
!
       if (.not. allocated(ad)) allocate(ad(neqm,neqm,IN,JN,KN))
       if (.not. allocated(a0)) allocate(a0(neqm,neqm,IN,JN,KN))
       if (.not. allocated(a1)) allocate(a1(neqm,neqm,IN,JN,KN))
       if (.not. allocated(a2)) allocate(a2(neqm,neqm,IN,JN,KN))
       if (.not. allocated(x)) allocate(x (     neqm,IN,JN,KN))
       if (.not. allocated(b)) allocate(b (     neqm,IN,JN,KN))
      endif ! igcall = 0
!
!-----------------------------------------------------------------------
!     compute BC arrays
!
!     NOTE: see comments in grav2dcg.F regarding use of subroutine
!           gpbv vs BC's specified through f**s(13) arrays.
!-----------------------------------------------------------------------
!
      call gpbv
!
!-----------------------------------------------------------------------
!     Scale Gravitational Constant out of BC arrays
!
!     NOTE: use only if gpbv is NOT called.  fiis(13), fois(13), etc.
!           must be in CGS units in zmp_inp.
!-----------------------------------------------------------------------
!
!     do k = ks, ke
!      do j = js, je
!       gpiib(j,k,1) = fiis(13) / guniv
!       gpiib(j,k,2) = gpiib(j,k,1)
!       gpoib(j,k,1) = fois(13) / guniv
!       gpoib(j,k,2) = gpoib(j,k,1)
!      enddo
!     enddo
!     do k = ks, ke
!      do i = is, ie
!       gpijb(i,k,1) = fijs(13) / guniv
!       gpijb(i,k,2) = gpijb(i,k,1)
!       gpojb(i,k,1) = fojs(13) / guniv
!       gpojb(i,k,2) = gpojb(i,k,1)
!      enddo
!     enddo
!     do j = js, je
!      do i = is, ie
!       gpikb(i,j,1) = fiks(13) / guniv
!       gpikb(i,j,2) = gpikb(i,j,1)
!       gpokb(i,j,1) = foks(13) / guniv
!       gpokb(i,j,2) = gpokb(i,j,1)
!      enddo
!     enddo
!
!-----------------------------------------------------------------------
!     Load arrays for matrix
!-----------------------------------------------------------------------
!
      do k = 1, kn
      do j = 1, jn
      do i = 1, in
       ad(1,1,i,j,k) = 0.0D0
       a0(1,1,i,j,k) = 0.0D0
       a1(1,1,i,j,k) = 0.0D0
       a2(1,1,i,j,k) = 0.0D0
       b (  1,i,j,k) = 0.0D0
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!     Compute diagonals along 1-faces
!-----------------------------------------------------------------------
!
      call compute_3d_diagonals(is, is, js, je, ks, ke, &
                                ad, a0, a1, a2)
      call compute_3d_diagonals(ie, ie, js, je, ks, ke, &
                                ad, a0, a1, a2)
!
!-----------------------------------------------------------------------
!     initiate MPI send/recv's for 1-face data
!-----------------------------------------------------------------------
!
      nreq = 0
      nsub = 0
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, ad)
!
!-----------------------------------------------------------------------
!     in the meantime, compute 2-face data
!-----------------------------------------------------------------------
!
      call compute_3d_diagonals(is+1, ie-1, js, js, ks, ke, &
                                ad, a0, a1, a2)
      call compute_3d_diagonals(is+1, ie-1, je, je, ks, ke, &
                                ad, a0, a1, a2)
!
!-----------------------------------------------------------------------
!     verify that 1-face data send/recv's have finished...
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
!
!-----------------------------------------------------------------------
!     and fire off the 2-face data
!-----------------------------------------------------------------------
!
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, ad)
!
!-----------------------------------------------------------------------
!     compute 3-face data
!-----------------------------------------------------------------------
!
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ks, ks, &
                                ad, a0, a1, a2)
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ke, ke, &
                                ad, a0, a1, a2)
!
!-----------------------------------------------------------------------
!     verify that 2-face data has been exchanged...
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
!
!-----------------------------------------------------------------------
!     ...and fire off the 3-face data
!-----------------------------------------------------------------------
!
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, ad)
!
!-----------------------------------------------------------------------
!     compute the matrix elements in all zones interior to the domain
!     faces
!-----------------------------------------------------------------------
!
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ks+1, ke-1, &
                                ad, a0, a1, a2)
!
!-----------------------------------------------------------------------
!     verify that 3-face data exchanges are finished
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
!
!-----------------------------------------------------------------------
!     compute RHS
!-----------------------------------------------------------------------
!
      DO k=ks,ke
         DO j=js,je
            DO i=is,ie
               b(1,i,j,k) = 4.0D0*d(i,j,k)* &
                             PI*guniv*dvl1a(i)*dvl2a(j)* &
                                      dvl3a(k)
            END DO
         END DO
      END DO
!
      epstol = 1.0D-12 ! * sqrt(float((ie-is+1)*(je-js+1)*(ke-ks+1)*
!     .                             ntiles(1)*ntiles(2)*ntiles(3) ) )
!
!-----------------------------------------------------------------------
!     modify matrix elements for boundary conditions as needed
!-----------------------------------------------------------------------
!
      if(coords(1) .eq. 0) then
       do k = ks, ke
        do j = js, je
         if((niis(3) .eq. 1) .or. (niis(3) .eq. 2)) then
           ad(1,1,is,j,k) = ad(1,1,is,j,k) + &
                            dvl2a(j)*dvl3a(k)*g2a(is  )* &
                            g31a(is  )*dx1bi(is  )
         endif
         if(niis(3) .eq. 3) then 
           b(1,is,j,k) = b(1,is,j,k) + guniv* &
           dvl2a(j)*dvl3a(k)*g2a(is  )*g31a(is)*dx1bi(is)*gpiib(j,k,1)
         endif
        enddo ! j
       enddo ! k
      endif ! coords(1)
      if(coords(1) .eq. ntiles(1)-1) then
       do k = ks, ke
        do j = js, je
         a0(1,1,ie,j,k) = 0.0D0
         if((nois(3) .eq. 1) .or. (nois(3) .eq. 2)) then
           ad(1,1,ie,j,k) = ad(1,1,ie,j,k) + &
                            dvl2a(j)*dvl3a(k)*g2a(ie+1)* &
                            g31a(ie+1)*dx1bi(ie+1)
         endif
         if(nois(3) .eq. 3) then
           b(1,ie,j,k) = b(1,ie,j,k) + guniv* &
           dvl2a(j)*dvl3a(k)*g2a(ie+1)*g31a(ie+1)*dx1bi(ie+1)* &
           gpoib(j,k,1)
         endif
        enddo ! j
       enddo ! k
      endif ! coords(1)
      if(coords(2) .eq. 0) then
       do k = ks, ke
        do i = is, ie
         if((nijs(3) .eq. 1) .or. (nijs(3) .eq. 2)) then
           ad(1,1,i,js,k) = ad(1,1,i,js,k) + &
                            dvl1a(i)*dvl3a(k)* &
                            g2bi(i)**2*g32a(js  )*dx2bi(js  )
         endif ! nijs
         if(nijs(3) .eq. 3) then
           b(1,i,js,k) = b(1,i,js,k) + guniv* &
           dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(js  )*dx2bi(js  )* &
           gpijb(i,k,1)
         endif ! nijs
        enddo ! i
       enddo ! k
      endif ! coords(2)
      if(coords(2) .eq. ntiles(2)-1) then
       do k = ks, ke
        do i = is, ie
         if(nojs(3) .ne. 4) a1(1,1,i,je,k) = 0.0D0
         if((nojs(3) .eq. 1) .or. (nojs(3) .eq. 2)) then
           ad(1,1,i,je,k) = ad(1,1,i,je,k) + &
                            dvl1a(i)*dvl3a(k)* &
                            g2bi(i)**2*g32a(je+1)*dx2bi(je+1)
         endif ! nojs
         if(nojs(3) .eq. 3) then
           b(1,i,je,k) = b(1,i,je,k) + guniv* &
           dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(je+1)*dx2bi(je+1)* &
           gpojb(i,k,1)
         endif ! nojs
        enddo ! i
       enddo ! k
      endif ! coords(2)
      if(coords(3) .eq. 0) then
       do j = js, je
        do i = is, ie
        if((niks(3) .eq. 1) .or. (niks(3) .eq. 2)) then
          ad(1,1,i,j,ks) = ad(1,1,i,j,ks) + &
                           dvl1a(i)*dvl2a(j)* &
                           g31bi(i)**2*g32bi(j)**2*dx3bi(ks  )
         endif
         if(niks(3) .eq. 3) b(1,i,j,ks) = b(1,i,j,ks) + guniv* &
          dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2*dx3bi(ks  )* &
          gpikb(i,j,1)
        enddo ! i
       enddo ! j
      endif ! coords(3)
      if(coords(3) .eq. ntiles(3)-1) then
       do j = js, je
        do i = is, ie
         if(noks(3) .ne. 4) a2(1,1,i,j,ke) = 0.0D0
         if((noks(3) .eq. 1) .or. (noks(3) .eq. 2)) then
          ad(1,1,i,j,ke) = ad(1,1,i,j,ke) + &
                           dvl1a(i)*dvl2a(j)* &
                           g31bi(i)**2*g32bi(j)**2*dx3bi(ke+1)
         endif
         if(noks(3) .eq. 3) b(1,i,j,ke) = b(1,i,j,ke) + guniv* &
          dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2*dx3bi(ke+1)* &
          gpokb(i,j,1)
        enddo ! i
       enddo ! k
      endif ! coords(3)
!
      call phicheck(ad, a0, a1, a2, b, gp, phip, epsmax, epsmin, &
                    imax, jmax, kmax)
!
      call mpi_allreduce(epsmax, glbeps, 1, MPI_DOUBLE_PRECISION, &
                         mpi_max, comm3d, ierr)
      call mpi_allreduce(epsmin, glbemin, 1, MPI_DOUBLE_PRECISION, &
                         mpi_min, comm3d, ierr)
!
      epsmax = glbeps
      epsmin = glbemin
!
!      if(myid .eq. 0) then
!       write(*,"('PHICHECK: max errphi = ',1pd12.4)")epsmax
!       write(*,"('imax, jmax, kmax = ',3i4)")imax,jmax,kmax
!       write(*,"('phip, gp = ',1p2d16.8)")phip(imax,jmax,kmax),
!     .                                    gp  (imax,jmax,kmax)
!      endif
!
      if(epsmax .lt. epstol) then
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       return
      else
       ncgcall = ncgcall + 1
      endif
!
!-----------------------------------------------------------------------
!     initialize potential
!-----------------------------------------------------------------------
!
      if(igcall .eq. 0) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
           x(1,i,j,k) = 0.0D0
         enddo
        enddo
       enddo
      else
       do k = ks, ke
        do j = js, je
         x(1,is,j,k) = -gp(is,j,k)
         x(1,ie,j,k) = -gp(ie,j,k)
        enddo
       enddo
       nreq = 0
       nsub = 0
       call updt_vec_bnd_1(is, ie, js, je, ks, ke, x)
       do k = ks, ke
        do i = is+1, ie-1
         x(1,i,js,k) = -gp(i,js,k)
         x(1,i,je,k) = -gp(i,je,k)
        enddo
       enddo
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
       call updt_vec_bnd_2(is, ie, js, je, ks, ke, x)
       do j = js+1, je-1
        do i = is+1, ie-1
         x(1,i,j,ks) = -gp(i,j,ks)
         x(1,i,j,ke) = -gp(i,j,ke)
        enddo
       enddo
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
       call updt_vec_bnd_3(is, ie, js, je, ks, ke, x)
       do k = ks+1, ke-1
        do j = js+1, je-1
         do i = is+1, ie-1
          x(1,i,j,k) = -gp(i,j,k)
         enddo
        enddo
       enddo
      endif ! igcall
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
!
!     -------------------------------------------------------
!     solve the system
!     -------------------------------------------------------
!
      igcall   = 1
!
      prev_tot = totcgit
      maxitr   = (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  * &
                 ntiles(1)  *  ntiles(2)  *  ntiles(3)
      toler    = 1.0D-6
!JH      toler    = 1.0D-8 * sqrt(float(
!JH     .           (ie-is+1)*(je-js+1)*(ke-ks+1)*
!JH     .           ntiles(1)*ntiles(2)*ntiles(3) ) )
!#ifdef 1
!      call MPI_ALLREDUCE(toler,glbtoler,1,MPI_DOUBLE_PRECISION,
!     .                   mpi_min,comm3d,ierr)
!#else
      glbtoler = toler
!#endif 
!
      call cgsolve(ad, a0, a1, a2, x, b, &
                   glbtoler, error)
      if(nits .ge. maxitr .and. error .gt. glbtoler &
             .and. myid_w .eq. 0) then
       write(*,"('**********  CGSOLVE did not converge with dt=', &
        1pe12.5,'  **********',/1x,'(nhy,error,toler,time,nits) =', &
        i7,3e14.5,i5)")dt,nhy,error,glbtoler,time,nits
       call mpi_finalize(ierr)
       stop
      else
!       if(myid .eq. 0) write(*,"('GRAVITY: cgsolve converged',
!     .    ' in ',i5,' iterations')")nits
       totcgit = totcgit + nits
      endif
!
!     -------------------------------------------------------
!     Finish up
!     -------------------------------------------------------
!
!     positive potential, consistent with ZEUS convention
!
      do k = ks, ke
       do j = js, je
        do i = is, ie
         gp(i,j,k) = -x(1,i,j,k)
        enddo
       enddo
      enddo
!
      nreq = 0
      nsub = nsub + 1
      call updt_gp_bnd_1(is,ie,js,je,ks,ke,gp)
!
!     scale BC arrays to CGS
!
      do k = ks, ke
       do j = js, je
        gpiib(j,k,1) = guniv*gpiib(j,k,1)
        gpiib(j,k,2) =       gpiib(j,k,1)
        gpoib(j,k,1) = guniv*gpoib(j,k,1)
        gpoib(j,k,2) =       gpoib(j,k,1)
       enddo
      enddo
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr )
       nreq = 0
      endif
      nsub = nsub + 1
      call updt_gp_bnd_2(is,ie,js,je,ks,ke,gp)
!
!     scale BC arrays to CGS
!
      do k = ks, ke
       do i = is, ie
        gpijb(i,k,1) = guniv*gpijb(i,k,1)
        gpijb(i,k,2) =       gpijb(i,k,1)
        gpojb(i,k,1) = guniv*gpojb(i,k,1)
        gpojb(i,k,2) =       gpojb(i,k,1)
       enddo
      enddo
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr )
       nreq = 0
      endif
      nsub = nsub + 1
      call updt_gp_bnd_3(is,ie,js,je,ks,ke,gp)
!
!     scale BC arrays to CGS
!
      do j = js, je
       do i = is, ie
        gpikb(i,j,1) = guniv*gpikb(i,j,1)
        gpikb(i,j,2) =       gpikb(i,j,1)
        gpokb(i,j,1) = guniv*gpokb(i,j,1)
        gpokb(i,j,2) =       gpokb(i,j,1)
       enddo
      enddo
!
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr )
       nreq = 0
      endif
      if(coords(1) .eq. 0) then
       do k = ks, ke
        do j = js, je
         if(niis(3) .eq. 1) gp(is-1,j,k) = gp(is,j,k)
         if(niis(3) .eq. 2) gp(is-1,j,k) = gp(is,j,k)
         if(niis(3) .eq. 3) gp(is-1,j,k) = gpiib(j,k,1)
        enddo
       enddo
      endif ! coords(1)
      if(coords(1) .eq. ntiles(1)-1) then
       do k = ks, ke
        do j = js, je
         if(nois(3) .eq. 1) gp(ie+1,j,k) = gp(ie,j,k)
         if(nois(3) .eq. 2) gp(ie+1,j,k) = gp(ie,j,k)
         if(nois(3) .eq. 3) gp(ie+1,j,k) = gpoib(j,k,1)
        enddo
       enddo
      endif ! coords(1)
      if(coords(2) .eq. 0) then
       do k = ks, ke
        do i = is, ie
         if(nijs(3) .eq. 1) gp(i,js-1,k) = gp(i,js,k)
         if(nijs(3) .eq. 2) gp(i,js-1,k) = gp(i,js,k)
         if(nijs(3) .eq. 3) gp(i,js-1,k) = gpijb(i,k,1)
        enddo
       enddo
      endif
      if(coords(2) .eq. ntiles(2)-1) then
       do k = ks, ke
        do i = is, ie
         if(nojs(3) .eq. 1) gp(i,je+1,k) = gp(i,je,k)
         if(nojs(3) .eq. 2) gp(i,je+1,k) = gp(i,je,k)
         if(nojs(3) .eq. 3) gp(i,je+1,k) = gpojb(i,k,1)
        enddo
       enddo
      endif
      if(coords(3) .eq. 0) then
       do j = js, je
        do i = is, ie
         if(niks(3) .eq. 1) gp(i,j,ks-1) = gp(i,j,ks)
         if(niks(3) .eq. 2) gp(i,j,ks-1) = gp(i,j,ks)
         if(niks(3) .eq. 3) gp(i,j,ks-1) = gpikb(i,j,1)
        enddo
       enddo
      endif
      if(coords(3) .eq. ntiles(3)-1) then
       do j = js, je
        do i = is, ie
         if(noks(3) .eq. 1) gp(i,j,ke+1) = gp(i,j,ke)
         if(noks(3) .eq. 2) gp(i,j,ke+1) = gp(i,j,ke)
         if(noks(3) .eq. 3) gp(i,j,ke+1) = gpokb(i,j,1)
        enddo
       enddo
      endif
!
!      if(myid .eq. 0)
!     .write(*,"('Phi(i): ',i3,1pd16.8)")(i,gp(i,je,ks),i=is,ie)
!
      return
      end
!
      subroutine compute_3d_diagonals(ibeg,iend,jbeg,jend,kbeg,kend, &
                                      ad, a0, a1, a2)
!
      use real_prec
      use param
      use grid
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
!
      real(rl) :: ad(1,1,in,jn,kn), a0(1,1,in,jn,kn), &
                  a1(1,1,in,jn,kn), a2(1,1,in,jn,kn)
!
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
!JH         a0(1,1,i,j,k) = dvl2a(j)*dvl3a(k)*g2a(i+1)*g31a(i+1)
!JH     .                 * dx1bi(i+1)
!JH         a1(1,1,i,j,k) = dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(j+1)
!JH     .                 * dx2bi(j+1)
         a0(1,1,i,j,k) = dvl2a(j)*dvl3a(k)*(g2a(i+1)*g31a(i+1))**2 &
                       * dvl1bi(i+1)
         a1(1,1,i,j,k) = dvl1a(i)*dvl3a(k)*(g2bi(i)*g32a(j+1))**2 &
                       * dvl2bi(j+1)
         a2(1,1,i,j,k) = dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2 &
                       * dx3bi(k+1)
         ad(1,1,i,j,k) = &
!     .           -dvl2a(j)*dvl3a(k)*( g2a(i+1)*g31a(i+1)*
!     .                                 dx1bi(i+1)+
!     .                                g2a(i  )*g31a(i  )*
!     .                                 dx1bi(i  ) )
!     .           -dvl1a(i)*dvl3a(k)*g2bi(i)**2*
!     .              (g32a(j)*dx2bi(j) + g32a(j+1)*dx2bi(j+1)) &
                 -dvl2a(j)*dvl3a(k)*( (g2a(i+1)*g31a(i+1))**2* &
                                       dvl1bi(i+1)+ &
                                      (g2a(i  )*g31a(i  ))**2* &
                                       dvl1bi(i  ) ) &
                 -dvl1a(i)*dvl3a(k)*g2bi(i)**2* &
                  (g32a(j)**2*dvl2bi(j) + g32a(j+1)**2*dvl2bi(j+1)) &
                 -dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2* &
                                    (dx3bi(k) + dx3bi(k+1))
        enddo
       enddo
      enddo
!
      return
      end
!
      subroutine phicheck(ad, a0, a1, a2, b, gp, phip, epsmax, &
                          epsmin, imax, jmax, kmax)
!
      use real_prec
      use param
      use grid
!
      implicit NONE
!
      integer  :: i, j, k, imax, jmax, kmax
!
      real(rl) :: ad(1,1,in,jn,kn), a0(1,1,in,jn,kn), &
                  a1(1,1,in,jn,kn), a2(1,1,in,jn,kn), &
                  b (1,  in,jn,kn)
!
      real(rl) :: gp(in,jn,kn), phip(in,jn,kn)
!
      real(rl) :: eps, epsmin, epsmax, glbeps, glbemin
!
      do k = ks, ke
       do j = js, je
        do i = is, ie
         phip(i,j,k) = -(1.0/ad(1,1,i,j,k))* &
                       ( b (1,  i  ,j  ,k  )                 + &
                         a2(1,1,i  ,j  ,k-1)*gp(i  ,j  ,k-1) + &
                         a2(1,1,i  ,j  ,k  )*gp(i  ,j  ,k+1) + &
                         a1(1,1,i  ,j-1,k  )*gp(i  ,j-1,k  ) + &
                         a1(1,1,i  ,j  ,k  )*gp(i  ,j+1,k  ) + &
                         a0(1,1,i-1,j  ,k  )*gp(i-1,j  ,k  ) + &
                         a0(1,1,i  ,j  ,k  )*gp(i+1,j  ,k  ) &
                        )
        enddo
       enddo
      enddo
!
      epsmax = -1.0d30
      epsmin =  1.0d30
      do k = ks, ke
       do j = js, je
        do i = is, ie
         eps = abs(phip(i,j,k)-gp(i,j,k))/(gp(i,j,k)+tiny)
         if(eps .gt. epsmax) then
          epsmax = eps
          imax   = i
          jmax   = j
          kmax   = k
         endif
         if(eps .lt. epsmin) then
          epsmin = eps
         endif
        enddo
       enddo
      enddo
!
      return
      end
