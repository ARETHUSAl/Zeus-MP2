!=======================================================================
!
!                            Developed by
!                Laboratory for Computational Astrophysics
!                  University of California at San Diego
!
      subroutine grav2D_CG
!
!     2-D gravitational potential solver using John Hayes's
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
      real(rl) :: toler, error, glbtoler
!
      integer :: i, j, k, prev_tot
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
       allocate(ad(neqm,neqm,IN,JN,KN))
       allocate(a0(neqm,neqm,IN,JN,KN))
       allocate(a1(neqm,neqm,IN,JN,KN))
       allocate(a2(neqm,neqm,IN,JN,KN))
       allocate(x (     neqm,IN,JN,KN))
       allocate(b (     neqm,IN,JN,KN))
!
      endif ! igcall
!
!-----------------------------------------------------------------------
!     compute BC arrays
!
!     Uncomment the "call gpbv" line if external BC's are to be computed
!     from a multipole approximation.  NOTE: gpbv uses a low-order
!     approximation which is not accurate if the external boundary is
!     not far from the mass distribution.  If the external boundary
!     potential can be computed analytically (as in the case of spherical
!     collapse), then those values should be adopted as fois(13) and
!     input in zmp_inp in the BC namelists as appropriate.
!-----------------------------------------------------------------------
!
!      call gpbv
      do k = 1, kn
       do j = 1, jn
        gpoib(j,k,1) = fois(13)
       enddo
      enddo
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
!     compute diagonals along i faces
!-----------------------------------------------------------------------
!
      call compute_diagonals(ad,a0,a1,is,is,js,je,ks,ke)
      call compute_diagonals(ad,a0,a1,ie,ie,js,je,ks,ke)
!
!-----------------------------------------------------------------------
!     exchange i-face values with neighboring processors
!-----------------------------------------------------------------------
!
      nreq = 0
      nsub = 0
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, ad)
!
!-----------------------------------------------------------------------
!     compute diagonals at interior zones
!-----------------------------------------------------------------------
!
      call compute_diagonals(ad,a0,a1,is+1,ie-1,js+1,je-1,ks,ke)
!
!-----------------------------------------------------------------------
!     wait for i-face MPI exchanges to finish
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
!
!-----------------------------------------------------------------------
!     compute diagonals along j faces
!-----------------------------------------------------------------------
!
      call compute_diagonals(ad,a0,a1,is+1,ie-1,js,js,ks,ke)
      call compute_diagonals(ad,a0,a1,is+1,ie-1,je,je,ks,ke)
!
!-----------------------------------------------------------------------
!     send j-face values to j-neighbor processes
!-----------------------------------------------------------------------
!
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, ad)
!
!-----------------------------------------------------------------------
!     compute RHS
!-----------------------------------------------------------------------
!
      DO k=ks,ke
         DO j=js,je
            DO i=is,ie
               b(1,i,j,k) = 4.0D0*d(i,j,k)* &
                             PI*guniv*dvl1a(i)*dvl2a(j)
            END DO
         END DO
      END DO
!
!-----------------------------------------------------------------------
!     wait for j-face exchanges to finish
!-----------------------------------------------------------------------
!
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
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
                            dvl2a(j)*g2a(is  )* &
                            g31a(is  )*dx1bi(is  )
         endif
         if(niis(3) .eq. 3) then 
           b(1,is,j,k) = b(1,is,j,k) + guniv* &
           dvl2a(j)*g2a(is  )*g31a(is)*dx1bi(is)*gpiib(j,k,1)
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
                            dvl2a(j)*g2a(ie+1)* &
                            g31a(ie+1)*dx1bi(ie+1)
         endif
         if(nois(3) .eq. 3) then
           b(1,ie,j,k) = b(1,ie,j,k) + guniv* &
           dvl2a(j)*g2a(ie+1)*g31a(ie+1)*dx1bi(ie+1)* &
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
                            dvl1a(i)* &
                            g2bi(i)**2*g32a(js  )*dx2bi(js  )
         endif ! nijs
         if(nijs(3) .eq. 3) then
           b(1,i,js,k) = b(1,i,js,k) + guniv* &
           dvl1a(i)*g2bi(i)**2*g32a(js  )*dx2bi(js  )* &
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
                            dvl1a(i)* &
                            g2bi(i)**2*g32a(je+1)*dx2bi(je+1)
         endif ! nojs
         if(nojs(3) .eq. 3) then
           b(1,i,je,k) = b(1,i,je,k) + guniv* &
           dvl1a(i)*g2bi(i)**2*g32a(je+1)*dx2bi(je+1)* &
           gpojb(i,k,1)
         endif ! nojs
        enddo ! i
       enddo ! k
      endif ! coords(2)
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
      prev_tot = totcgit
      maxitr   = (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  * &
                 ntiles(1)  *  ntiles(2)  *  ntiles(3)
      toler    = 1.0D-8 ! *
!     .           (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
!     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)  
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
       totcgit = totcgit + nits
       ncgcall = ncgcall + 1
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
!
      igcall = 1
!
      return
      end
!
      subroutine compute_diagonals(ad,a0,a1,is,ie,js,je,ks,ke)
!
      use real_prec
      use param, ONLY: in,jn,kn
      use grid,  ONLY: g2a, g31a, g2bi, g32a, &
                       dx1bi, dx2bi, dvl1a, dvl2a
!
      implicit NONE
!
      integer  :: i, j, k, is, ie, js, je, ks, ke
!
      real(rl) :: ad(1,1,in,jn,kn), &
                  a0(1,1,in,jn,kn), &
                  a1(1,1,in,jn,kn)
!
      do k = ks, ke
       do j = js, je
        do i = is, ie
         a0(1,1,i,j,k) = dvl2a(j)*g2a(i+1)*g31a(i+1) &
                       * dx1bi(i+1)
         a1(1,1,i,j,k) = dvl1a(i)*g2bi(i)**2*g32a(j+1) &
                       * dx2bi(j+1)
         ad(1,1,i,j,k) = &
                 -dvl2a(j)*(g2a(i+1)*g31a(i+1)*dx1bi(i+1)+ &
                                     g2a(i  )*g31a(i  )*dx1bi(i  )) &
                 -dvl1a(i)*g2bi(i)**2* &
                            (g32a(j)*dx2bi(j)+g32a(j+1)*dx2bi(j+1))
        enddo
       enddo
      enddo
!
      return
      end
