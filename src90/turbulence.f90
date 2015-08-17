!
! JReppin 04/08/15
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////            T U R B U L E N C E            \\\\\\\\\!=======================================================================
!
subroutine turbulence
!
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use mpiyes
      use mpipar
      implicit none
!
      integer :: i,j,k
!      real(rl), dimension(in,jn,kn) :: a1 , a2, a3
      real(rl)  :: t_unit , box_L , l_unit, rho_unit , v_unit , b_unit
      real(rl)  ::  rho0 , b0, v0, e0, p0
      integer   :: nmodes_B, nmodes_V
      real(rl)  ::  idx_A, idx_V, helic
      !real(rl)  ::  diff_eta, diff_nu 
      namelist / pgen  / &
            l_unit, v_unit, rho_unit, &
!            diff_eta, diff_nu, &
            rho0, b0, v0, e0, p0, &
            nmodes_B, idx_A, nmodes_V, idx_V, helic
!
!      default values
!      fix all units according to the values of rho, c_s at z = 0 
!
      nmodes_B     = 2
      nmodes_V     = 2
      idx_A        = 0.0
      idx_V        = 0.0
      helic        = 0.0
      rho0         = 1.0
      v0           = 0.1
      e0           = 1.0
      p0           = 1.0
      b0           = 0.1
      l_unit       = 3.086e20
      v_unit       = 19875.8   ! = sqrt(1./mu_mh)*15005.9         
      rho_unit     = 2.3485e-31
!      diff_nu      = 1.e-3     ! kinetic diffusion constant
!      diff_eta     = 1.e-3     ! magnetic diffusion constant
!
!     read in namelist and write it to log file
!
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
        ibuf_in( 1) = nmodes_B 
        ibuf_in( 2) = nmodes_V
        buf_in( 1)  = idx_A
        buf_in( 2)  = idx_V
        buf_in( 3)  = helic  
        buf_in( 4)  = rho0  
        buf_in( 5)  = v0   
        buf_in( 6)  = e0   
        buf_in( 7)  = p0       
        buf_in( 8)  = b0        
        buf_in( 9)  = l_unit    
        buf_in(10)  = v_unit   
        buf_in(11)  = rho_unit
!        buf_in(12)  = diff_nu  
!        buf_in(13)  = diff_eta
      endif
      call MPI_BCAST( ibuf_in, 2, MPI_INTEGER &
                     , 0, comm3d, ierr )
      call MPI_BCAST(ibuf_in, 13, MPI_DOUBLE_PRECISION &
                     , 0, comm3d, ierr )
      if (myid .ne. 0) then
        nmodes_B = ibuf_in( 1)
        nmodes_v = ibuf_in( 2)
        idx_A    = buf_in( 1)
        idx_V    = buf_in( 2)
        helic    = buf_in( 3)
        rho0     = buf_in( 4)
        v0       = buf_in( 5)
        e0       = buf_in( 6)
        p0       = buf_in( 7)
        b0       = buf_in( 8)
        l_unit   = buf_in( 9)
        v_unit   = buf_in(10)
        rho_unit = buf_in(11)
!        diff_nu  = buf_in(12)
!        diff_eta = buf_in(13)
      endif
!
!   Set other physical constants
!
      t_unit    =  l_unit/v_unit
      b_unit    =  sqrt(4*pi*rho_unit*v_unit**2)       ! in Gauss !! 
      if(myid.eq.0)then 
        write(6,2030) 'setup for MHD turbulence'
        write(6,2030) ' '
        write(6,2040) 'b0         : ',b0
        write(6,2040) 'rho0       : ',rho0
        write(6,2040) 'v0         : ',v0
        write(6,2040) 'Starttime              : ',time
        write(6,2040) 'Endtime                : ',tlim
        write(6,2030) '------------------------------------------'
        write(6,2030) 'Units:'
        write(6,2040) 'rho_u                  : ',rho_unit
        write(6,2040) 'L_u                    : ',l_unit
        write(6,2040) 'V_u                    : ',v_unit
        write(6,2040) 't_u                    : ',t_unit
        write(6,2040) 'b_u                    : ',b_unit
        write(6,2030) '-------------------------------------------'
        write(6,2050) 'excited modes B-field  : ', nmodes_B
        write(6,2050) 'excited modes V-field  : ', nmodes_V
        write(6,2040) 'Spectral idx   A-field : ', idx_A
        write(6,2040) 'Spectral idx vel-field : ', idx_V
        write(6,2040) 'helic                  : ', helic
!        write(6,2040) 'kinetic diff.-const.   : ',diff_nu
!        write(6,2040) 'magnetic diff.-const   : ',diff_eta
        write(6,2030) '  '
      endif
!
! set up density and other fields
!
      b0 = b0*  1.0D0 / sqrt(4.0D0*pi)
      do k=1,kn
        do j=1,jn
          do i=1,in
            d(i,j,k) = rho0
            e(i,j,k) = e0
            !v1(i,j,k) = 0
            !v2(i,j,k) = 0
            !v3(i,j,k) = 0
          enddo 
        enddo
      enddo
!
! generate stochastic velocity field
!
      call genvel(nmodes_V, v0, idx_V)
!
! generate stochastic B-field
!        b0 = b0/sqrt(mu_0c2)
      call genhelic(nmodes_B, b0, idx_A, helic)
   
2010  format('TURB   : Initialization complete.')
2020  format('TURB   : csiso(1) = ',1pe13.6)
2030  format('TURB   : ',a50)
2040  format('TURB   : ',a25,1pe13.6)
2050  format('TURB   : ',a25,i7)
      write(*,*) maxval(v1), maxval(b1)
      write(*,*) maxval(v2), maxval(b2)
      write(*,*) maxval(v3), maxval(b3)
      return
end subroutine turbulence
