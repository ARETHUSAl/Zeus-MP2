subroutine power_randomphase_hel(ampl,initpower,initpower2, &
      cutoff,ncutoff,kpeak,relhel,kgaussian, &
      lskip_projection,lno_second_ampl,lscale_tobox)
!
!  Produces helical k^initpower*exp(-k**2/cutoff**2) spectrum.
!  The relative helicity is relhel.
!
!  initpower=0 -> k^2  (if used for vector potential)
!           =2 -> k^4
!          =-3 -> k^{-1}
!          =-3.67 k^{-5/3}
!          =-5 -> k^{-3}
!
!  08-sep-14/axel: adapted from power_randomphase
!
      use param
      use root
      use real_prec
      use field
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
!
      logical, intent(in), optional :: lscale_tobox
      logical :: lscale_tobox1, lskip_projection, lno_second_ampl, lroot
      integer :: i,i1,i2,ikx,iky,ikz,stat
      real(rl), dimension (:,:,:,:), allocatable :: u_re, u_im, v_re, v_im
      real(rl), dimension (:,:,:), allocatable :: k2, r
      real(rl), dimension (:), allocatable :: kx, ky, kz
!
      if (present(lscale_tobox)) then
        lscale_tobox1 = lscale_tobox
      else
        lscale_tobox1 = .false.
      endif

      lroot = .false.
#ifdef MPI_USED
      if (myid_w .eq. 0) lroot = .true.
#endif
!
!  Allocate memory for arrays.
!
      allocate(k2(in,jn,kn),stat=stat)
      allocate(r(in,jn,kn),stat=stat)
!
      allocate(u_re(in,jn,kn,3),stat=stat)
      allocate(u_im(in,jn,kn,3),stat=stat)
!
      allocate(v_re(in,jn,kn,3),stat=stat)
      allocate(v_im(in,jn,kn,3),stat=stat)
!
      allocate(kx(in),stat=stat)
      allocate(ky(jn),stat=stat)
      allocate(kz(kn),stat=stat)
!
      if (ampl==0) then
        f(:,:,:,:) = 0.0
        if (lroot) print*,'power_randomphase: set variable to zero;'
      else
!
!  calculate k^2
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lx
        kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Ly
        ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lz
        kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*scale_factor
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!
        if (lroot ) &
             print*,'power_randomphase:fft done; now integrate over shells...'
        do iky=1,kn
          do ikx=1,jn
            do ikz=1,in
              k2(ikz,ikx,iky)=kx(ikx+ipy*ny)**2+ky(iky+ipz*nz)**2+kz(ikz+ipx*nx)**2
            enddo
          enddo
        enddo
        if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
!
!  generate flat spectrum with random phase (between -pi and pi)
!
        do i=1,3
          call random_number_wrapper(r)
          u_re(:,:,:,i)=ampl*k2**mhalf*cos(pi*(2*r-1))
          u_im(:,:,:,i)=ampl*k2**mhalf*sin(pi*(2*r-1))
        enddo !i
!
!  To get the shell integrated power spectrum E ~ k^n, we need u ~ k^m
!  and since E(k) ~ u^2 k^2 we have n=2m+2, so m=n/2-1.
!  Further, since we operate on k^2, we need m/2 (called mhalf below)
!
        mhalf=.5*(.5*initpower-1)
!
!  generate all 3 velocity components separately
!  generate k^n spectrum with random phase (between -pi and pi)
!
        nexp1=.25*nfact*(initpower-initpower2)
        nexp2=1./nfact
        kpeak21=1./kpeak**2
!
!  By mistake, a second ampl factor got introduced, so
!  the resulting amplitudes depended quadratically on ampl.
!  To keep continuity, we must keep this, but we can avoid
!  this by setting lno_second_ampl=T in the call.
!
        if (lno_second_ampl) then
          r=((k2*kpeak21)**mhalf)/(1.+(k2*kpeak21)**nexp1)**nexp2
        else
          r=ampl*((k2*kpeak21)**mhalf)/(1.+(k2*kpeak21)**nexp1)**nexp2
        endif
!
!  cutoff (changed to hyperviscous cutoff filter)
!
        if (cutoff /= 0.) r=r*exp(-(k2/cutoff**2.)**ncutoff)
!
!  apply Gaussian on top of everything
!
        if (kgaussian /= 0.) r=r*exp(-.5*(k2/kgaussian**2.-1.))
!
!  scale with r
!
        do i=1,3
          u_re(:,:,:,i)=r*u_re(:,:,:,i)
          u_im(:,:,:,i)=r*u_im(:,:,:,i)
        enddo !i
!
!  Apply projection operator
!  Use r=1/k^2 for normalization in khat_i * khat_j = ki*kj/k2.
!  Remember that for the return transform, data have to be
!  arranged in the order (kz,kx,ky).
!
        if (lskip_projection) then
          v_re=u_re
          v_im=u_im
        else
        do iky=1,kn
          do ikx=1,jn
            do ikz=1,in
!
!  Real part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
              r(ikz,ikx,iky)=(kx(ikx+ipy*ny)*u_re(ikz,ikx,iky,1) &
                             +ky(iky+ipz*nz)*u_re(ikz,ikx,iky,2) &
                             +kz(ikz+ipx*nx)*u_re(ikz,ikx,iky,3))/k2(ikz,ikx,iky)
              v_re(ikz,ikx,iky,1)=u_re(ikz,ikx,iky,1)-kx(ikx+ipy*ny)*r(ikz,ikx,iky)
              v_re(ikz,ikx,iky,2)=u_re(ikz,ikx,iky,2)-ky(iky+ipz*nz)*r(ikz,ikx,iky)
              v_re(ikz,ikx,iky,3)=u_re(ikz,ikx,iky,3)-kz(ikz+ipx*nx)*r(ikz,ikx,iky)
!
!  Imaginary part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
              r(ikz,ikx,iky)=(kx(ikx+ipy*ny)*u_im(ikz,ikx,iky,1) &
                             +ky(iky+ipz*nz)*u_im(ikz,ikx,iky,2) &
                             +kz(ikz+ipx*nx)*u_im(ikz,ikx,iky,3))/k2(ikz,ikx,iky)
              v_im(ikz,ikx,iky,1)=u_im(ikz,ikx,iky,1)-kx(ikx+ipy*ny)*r(ikz,ikx,iky)
              v_im(ikz,ikx,iky,2)=u_im(ikz,ikx,iky,2)-ky(iky+ipz*nz)*r(ikz,ikx,iky)
              v_im(ikz,ikx,iky,3)=u_im(ikz,ikx,iky,3)-kz(ikz+ipx*nx)*r(ikz,ikx,iky)
!
            enddo
          enddo
        enddo
        endif
!
!  Make it helical, i.e., multiply by delta_ij + epsilon_ijk ik_k*sigma.
!  Use r=sigma/k for normalization of sigma*khat_i = sigma*ki/sqrt(k2).
!
        r=relhel/sqrt(k2)
        do iky=1,kn
          do ikx=1,jn
            do ikz=1,in
!
!  (vx, vy, vz) -> ux
!
              u_re(ikz,ikx,iky,1)=v_re(ikz,ikx,iky,1) &
                  +kz(ikz+ipx*nx)*v_im(ikz,ikx,iky,2)*r(ikz,ikx,iky) &
                  -ky(iky+ipz*nz)*v_im(ikz,ikx,iky,3)*r(ikz,ikx,iky)
              u_im(ikz,ikx,iky,1)=v_im(ikz,ikx,iky,1) &
                  -kz(ikz+ipx*nx)*v_re(ikz,ikx,iky,2)*r(ikz,ikx,iky) &
                  +ky(iky+ipz*nz)*v_re(ikz,ikx,iky,3)*r(ikz,ikx,iky)
!
!  (vx, vy, vz) -> uy
!
              u_re(ikz,ikx,iky,2)=v_re(ikz,ikx,iky,2) &
                  +kx(ikx+ipy*ny)*v_im(ikz,ikx,iky,3)*r(ikz,ikx,iky) &
                  -kz(ikz+ipx*nx)*v_im(ikz,ikx,iky,1)*r(ikz,ikx,iky)
              u_im(ikz,ikx,iky,2)=v_im(ikz,ikx,iky,2) &
                  -kx(ikx+ipy*ny)*v_re(ikz,ikx,iky,3)*r(ikz,ikx,iky) &
                  +kz(ikz+ipx*nx)*v_re(ikz,ikx,iky,1)*r(ikz,ikx,iky)
!
!  (vx, vy, vz) -> uz
!
              u_re(ikz,ikx,iky,3)=v_re(ikz,ikx,iky,3) &
                  +ky(iky+ipz*nz)*v_im(ikz,ikx,iky,1)*r(ikz,ikx,iky) &
                  -kx(ikx+ipy*ny)*v_im(ikz,ikx,iky,2)*r(ikz,ikx,iky)
              u_im(ikz,ikx,iky,3)=v_im(ikz,ikx,iky,3) &
                  -ky(iky+ipz*nz)*v_re(ikz,ikx,iky,1)*r(ikz,ikx,iky) &
                  +kx(ikx+ipy*ny)*v_re(ikz,ikx,iky,2)*r(ikz,ikx,iky)
!
            enddo
          enddo
        enddo
!
!  back to real space
!
        do i=1,3
          call fourier_transform(u_re(:,:,:,i),u_im(:,:,:,i),linv=.true.)
        enddo !i
        f(l1:l2,m1:m2,n1:n2,i1:i2)=u_re
!
!  notification
!
        if (lroot) then
          if (cutoff==0) then
            print*,'powern: k^',initpower,' spectrum : var  i=',i
          else
            print*,'powern: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                initpower,', k0 =',cutoff,' : var  i=',i
          endif
        endif
!
      endif !(ampl==0)
!
!  Deallocate arrays.
!
      if (allocated(k2))   deallocate(k2)
      if (allocated(r))  deallocate(r)
      if (allocated(u_re)) deallocate(u_re)
      if (allocated(u_im)) deallocate(u_im)
      if (allocated(v_re)) deallocate(v_re)
      if (allocated(v_im)) deallocate(v_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
endsubroutine power_randomphase_hel
!***********************************************************************
