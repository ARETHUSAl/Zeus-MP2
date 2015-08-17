!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               O P T I O N S               \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine options
!
      use param
      use config
      use mpiyes
      use mpipar
!
      implicit NONE
!
      character*39 :: GEOMETRY_OPTION
      character*39 :: GRID_OPTION      
      character*39 :: PT_MASS_OPTION   
      character*39 :: GRAVITY_OPTION   
      character*39 :: PS_OPTION      
      character*39 :: MHD_OPTION      
      character*39 :: RAD_OPTION   
      character*39 :: CHEM_OPTION
      character*39 :: AV_OPTION   
      character*39 :: MPI_OPTION      
      character*39 :: RES_DUMPS_OPTION 
      character*39 :: HDF_DUMP_OPTION 
      character*39 :: HST_DUMP_OPTION
      character*39 :: TSL_DUMPS_OPTION
      character*39 :: TXT_DUMPS_OPTION
!
      if(lgeom .eq. 1) then
      GEOMETRY_OPTION  = "              * geometry:        XYZ   "
      endif
      if(lgeom .eq. 2) then
      GEOMETRY_OPTION  = "              * geometry:        ZRP   "
      endif
      if(lgeom .eq. 3) then
      GEOMETRY_OPTION  = "              * geometry:        RTP   "
      endif
!
      if(xvgrid) then
      GRID_OPTION      = "              * moving grid      ON    "
      else
      GRID_OPTION      = "              * moving grid      OFF   "
      endif
!
      if(xptmass) then
      PT_MASS_OPTION   = "              * point masses     ON    "
      else
      PT_MASS_OPTION   = "              * point masses     OFF   "
      endif
!
      if(xgrav .or. xgrvfft) then
      GRAVITY_OPTION   = "              * self-gravity     ON    "
      if(xgrvfft) then
       PS_OPTION       = "              * Poisson solver: FFTW   "
      else
       PS_OPTION       = "              * Poisson solver: MGMPI  "
      endif
      else
      GRAVITY_OPTION   = "              * self-gravity     OFF   "
      endif
!
      if(xmhd) then
      MHD_OPTION       = "              * magnetic fields  ON    "
      else
      MHD_OPTION       = "              * magnetic fields  OFF   "
      endif
!
      if(lrad .eq. 0) then
      RAD_OPTION       = "              * rad transport    OFF   "
      endif
      if(lrad .eq. 1) then
      RAD_OPTION       = "              * rad transport = greyFLD"
      endif
      if(lrad .eq. 2) then
      RAD_OPTION       = "              * rad transport = MGFLD  "
      endif
      if(lrad .eq. 3) then
      RAD_OPTION       = "              * rad transport = MGVTEF "
      endif
      if(xchem) then
      CHEM_OPTION      = "              * chemistry        ON    "
      else
      CHEM_OPTION      = "              * chemistry        OFF   "
      endif
!
      if(xsubav) then
       AV_OPTION       = "              * A.V. sub-cycling ON    "
      else
       AV_OPTION       = "              * A.V. sub-cycling OFF   "
      endif
!
      MPI_OPTION       = "              * message passing  ON    "
!
      if(xrestart) then
      RES_DUMPS_OPTION = "              * restart dumps    ON    "
      else
      RES_DUMPS_OPTION = "              * restart dumps    OFF   "
      endif
!
      if(xhdf) then
      HDF_DUMP_OPTION =  "              * HDF5 VIZ dumps   ON    "
      else
      HDF_DUMP_OPTION =  "              * HDF VIZ dumps    OFF   "
      endif
!
      if(xhst) then
      HST_DUMP_OPTION = "              * history dump     ON    "
      else
      HST_DUMP_OPTION = "              * history dump     OFF   "
      endif
!
      if(xtsl) then
      TSL_DUMPS_OPTION = "              * TSL dumps        ON    "
      else
      TSL_DUMPS_OPTION = "              * TSL dumps        OFF   "
      endif
!
      if(xascii) then
      TXT_DUMPS_OPTION = "              * text dumps       ON    "
      else
      TXT_DUMPS_OPTION = "              * text dumps       OFF   "
      endif
!
!  write greeting.
!
      write(6,"(///10x,'ZZZZZ EEEEE U   U  SSSS     M   M PPPP ')")
      write(6,   "(10x,'   Z  E     U   U S         MM MM P   P')")
      write(6,   "(10x,'  Z   EEEE  U   U  SSS  === M M M PPPP ')")
      write(6,   "(10x,' Z    E     U   U     S     M   M P    ')")
      write(6,   "(10x,'ZZZZZ EEEEE  UUU  SSSS      M   M P    ')")
      write(6,"()")
      write(6,"()")
      write(6,   "(10x,'    DEVELOPED BY ROBERT A. FIEDLER')")
      write(6,   "(10x,'      ZEUS-MP V2.1.2 - 1/25/07')")
      write(6,"()")
      write(6,   "(10x,'    RUNNING AS ',  i4,' PROCESS(ES)')")nprocs_w
      write(6,   "(10x,'   WITH THE FOLLOWING CONFIGURATION:')")
      write(6,"()")
      write(6,1) GEOMETRY_OPTION
      write(6,1) GRID_OPTION  
      write(6,1) PT_MASS_OPTION   
      write(6,1) GRAVITY_OPTION   
      if(xgrvfft) then
       write(6,1) PS_OPTION   
      endif
      write(6,1) MHD_OPTION       
      write(6,1) RAD_OPTION 
      write(6,1) CHEM_OPTION
      write(6,1) MPI_OPTION 
      write(6,1) RES_DUMPS_OPTION 
      write(6,1) HDF_DUMP_OPTION 
      write(6,1) HST_DUMP_OPTION 
      write(6,1) TSL_DUMPS_OPTION 
      write(6,1) TXT_DUMPS_OPTION 
      write(6,"()")
!
1     format(a39)
!
      return
      end
