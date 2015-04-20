c
c JReppin 04/08/15
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            T U R B U L E N C E            \\\\\\\\\\
c=======================================================================
c
      subroutine turbulence
c
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
      implicit none
c
      integer :: i,j,k
c      real(rl), dimension(in,jn,kn) :: a1 , a2, a3
      real(rl)  :: t_unit , box_L , l_unit, rho_unit , v_unit , b_unit
      real(rl)  ::  rho0 , b0, v0, e0, p0
      integer   :: nmodes_B, nmodes_V
      real(rl)  ::  idx_A, idx_V, helic
      !real(rl)  ::  diff_eta, diff_nu 

      namelist / pgen  /
     &      l_unit, v_unit, rho_unit, 
c     &      diff_eta, diff_nu, 
     &      rho0, b0, v0, e0, p0,
     &      nmodes_B, idx_A, nmodes_V, idx_V, helic
c
c      default values
c      fix all units according to the values of rho, c_s at z = 0 
c
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
c      diff_nu      = 1.e-3     ! kinetic diffusion constant
c      diff_eta     = 1.e-3     ! magnetic diffusion constant
c
c     read in namelist and write it to log file
c
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
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
c        buf_in(12)  = diff_nu  
c        buf_in(13)  = diff_eta
      endif
      call MPI_BCAST( ibuf_in, 2, MPI_INTEGER
     &               , 0, comm3d, ierr )
      call MPI_BCAST(ibuf_in, 13, MPI_FLOAT
     &               , 0, comm3d, ierr )
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
c        diff_nu  = buf_in(12)
c        diff_eta = buf_in(13)
#endif /* MPI_USED */
      endif
c
c   Set other physical constants
c
      t_unit    =  l_unit/v_unit
      b_unit    =  sqrt(4*pi*rho_unit*v_unit**2)       ! in Gauss !! 

#ifdef MPI_USED
      if(myid.eq.0)then 
        write(6,2030) 'setup for MHD turbulence'
        write(6,2030) ' '
        write(6,2040) 'b0         : ',b0
        write(6,2040) 'rho0       : ',rho0
        write(6,2040) 'v0         : ',v0
#ifdef ISO
        write(6,2040) 'e0         : ',e0
        write(6,2040) 'p0         : ',p0
#endif /* ISO */
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
c        write(6,2040) 'kinetic diff.-const.   : ',diff_nu
c        write(6,2040) 'magnetic diff.-const   : ',diff_eta
        write(6,2030) '  '
      endif
#endif /* MPI_USED */
c
c set up density and other fields
c
      d(:,:,:) = rho0
      e(:,:,:) = e0
      v0 = v0 + tiny
      b0 = b0 + tiny
c
c generate stochastic velocity field
c
      call genvel(nmodes_V, v0, idx_V)
c
c generate stochastic B-field
c        b0 = b0/sqrt(mu_0c2)
      call genhelic(nmodes_B, b0, idx_A, helic)
   

2010  format('TURB   : Initialization complete.')
2020  format('TURB   : csiso(1) = ',1pe13.6)
2030  format('TURB   : ',a50)
2040  format('TURB   : ',a25,1pe13.6)
2050  format('TURB   : ',a25,i7)


c       return
      end subroutine
