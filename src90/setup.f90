!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine setup
!
!  PURPOSE: Sets up execution of a new run by initializing all variables
!  according to the flags and values in the input deck.  Calls benchmark
!  a CPP macro) to initialize field variables for the particular
!  problem to be studied, otherwise field variables are set to "tiny"
!  (parameter defined to be smallest number possible on machine).  Note
!  user must define benchmark to be the appropriate subroutine name in the
!  file zeusmp.def or in the Make_zeusmp command line.
!
!  Order in which namelists are read has been changed from ZEUS-2D so 
!  that the boundary conditions are known before the MPI virtual 
!  topology is defined and before the grid is computed.  MPI needs to 
!  know if the grid is periodic, while the grid should be symmetric 
!  across reflecting boundaries but constant across flow in/out and 
!  periodic for periodic ones.
!
!  EXTERNALS:
!     ggen    -- initializes grid according to input deck
!     bval*   -- boundary value routines
!     nudt    -- computes initial timestep
!     newgrid -- computes new grid position for moving grid
!     scopy
!
!  LOCALS:
!
!  LAST MODIFIED: JCH 3/13/97.
!
!               : DJW 11-28-06 for 9-species primordial chemistry/RT
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use bndry
      use radiation
      use opac
      use opac_law
      use scratch
      use gravmod
      use mpiyes
      use mpipar
      use cons
      use impsoln
      use chem
!
      implicit NONE
!
      integer  :: i, j, k, n, jone, kone
      integer  :: iord,istp
      real(rl) :: dtrat
!
      integer imax,ISMAX, iost
      real(rl) :: dmax
!
      namelist /pcon/ nlim,tlim,cpulim,tsave,mbatch
      namelist /hycon/ &
       qcon,qlin,courno,dtrat,iord,istp,iCMA, &
       iordd,iorde,iords1,iords2,iords3,iordb1,iordb2,iordb3,iorder, &
       istpd,istpe,istps1,istps2,istps3,istpb1,istpb2,istpb3,istper, &
       dfloor,efloor,v1floor,v2floor,v3floor,b1floor,b2floor,b3floor, &
       emf1floor,emf2floor,emf3floor,erfloor,gpfloor
!
      namelist /rtcon/ ifreq, ispct, iOTS, iPWA, iLW, IHM, ibkgnd, &
                       iextinct, nphdot, nsrc, ephot, t_on, t_off, &
                       r_sep, t_star, hnumax, alpha, nLWdot, L_star, &
                       J_21
!
      namelist /chemcon/ iceco, icico, ireco, ibrco, ipiht, ih2co, &
                         icmpt, iDMco, iHDco, icycle, t_cutoff, rdshft, &
                         fh   , idust, z_sol
!
      namelist /iib/     niis, fiis
      namelist /oib/     nois, fois
      namelist /ijb/     nijs, fijs
      namelist /ojb/     nojs, fojs
      namelist /ikb/     niks, fiks
      namelist /okb/     noks, foks
!
      namelist /mpitop/ ntiles, periodic
      namelist /grvcon/ guniv,tgrav,ptmass,x1ptm,x2ptm,x3ptm,xwedge, &
                        gsup
      namelist /radcon/ ifld,epsme,demax,dermax,nmeiter,radth,epsrad, &
           cnvcrit,ernom,ennom,epsmaxd,cgerrcrit, &
           ipcflag,rmfp0,xnu,powr,rho0,t_0,kpfrac
      namelist /eqos/ gamma, ciso, mmw
      namelist /gcon/ x1fac,x2fac,x3fac,iga,jga,kga,igcon
      external ctable
!=======================================================================
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
!
!------------------------  benchmark CONTROL  ---------------------------
!
!   nlim   = cycles to run                               (default    1M)
!   tlim   = physical (problem) time to stop calculation (default  0.0 )
! cpulim   = CPU time in seconds to stop the calculation (default  3.6M)
!  tsave   = CPU time to reserve for terminating the run (default 30.0s)
! mbatch   = 0 interactive mode                          (default     0)
!          = 1 batch mode (does not scan for keyboard input)
!
      nlim   = 1000000
      tlim   = 0.0
      cpulim = 3600000.0
      tsave  = 30.0
      mbatch = 0
!
      nred   = 0
      if (myid_w .eq. 0) then
        open (1,file='zmp_inp', status='old', iostat=iost)
        read (1,pcon)
        write(2,pcon)
      endif
!
!------------------------  HYDRO CONTROL  ------------------------------
!
!  qcon   = quadratic artificial viscosity (q) constant (default 2.0)
!  qlin   = linear    artificial viscosity (q) constant (default 0.0)
!  courno = courant number                              (default 0.5)
!  dtrat  = ratio of initial dt to dtmin (used to compute dtmin below)
!  iord   = default order of advection scheme for all variables
!  iCMA   = use consistent advection scheme of Plewa & Mueller, A&A, 342
!           179 to correct abundance conservation problems not addressed
!           by the MCONSRV module in some 9-species I-front problems. 
!           iCMA=0 is no consistent advection; iCMA=1 means H and He are
!           separately enforced to be fh and 1-fh; iCMA=2 means the sum
!           of the H, He, and metal abundances are required to be 1.
!  istp   = default steepening flag for all variables
!  iord** = order of advection scheme to be used for variable **
!  iostp**= steepening flag for 3rd order advection.  When istp**=1,
!           use the discontinuity detection to steepen shocks during
!           interpolation for variable ** in X1INT,X1INTFC,X2INT,X2INTFC
!  **floor = smallest value desired for variable ** on grid
!            Can also be used to set a default value for initialization.
!            Note that no attempt is made to ensure that actual values 
!            stay above the floor values.  
!
      if (myid_w .eq. 0) then
        qcon   = 2.0
        qlin   = 0.0
        courno = 0.5
        dtrat  = 1.0e-3
        iord   = 2
        iCMA   = 0
        iordd  = 0
        iorde  = 0
        iords1 = 0
        iords2 = 0
        iords3 = 0
        iordb1 = 0
        iordb2 = 0
        iordb3 = 0
        iorder = 0
        istp   = 0
        istpd  = 2
        istpe  = 2
        istps1 = 2
        istps2 = 2
        istps3 = 2
        istpb1 = 2
        istpb2 = 2
        istpb3 = 2
        istper = 2
        dfloor = tiny
        efloor = tiny
        v1floor = 0.0
        v2floor = 0.0
        v3floor = 0.0
        b1floor = 0.0
        b2floor = 0.0
        b3floor = 0.0
        emf1floor= 0.0
        emf2floor= 0.0
        emf3floor= 0.0
        erfloor = tiny
        gpfloor = 0.0
!
        read (1,hycon)
        write(2,hycon)
!
! Set flags to default values unless they were set in the input deck.
!
        if(iordd  .eq. 0) iordd  = iord
        if(iorde  .eq. 0) iorde  = iord
        if(iords1 .eq. 0) iords1 = iord
        if(iords2 .eq. 0) iords2 = iord
        if(iords3 .eq. 0) iords3 = iord
        if(iordb1 .eq. 0) iordb1 = iord
        if(iordb2 .eq. 0) iordb2 = iord
        if(iordb3 .eq. 0) iordb3 = iord
        if(iorder .eq. 0) iorder = iord
!
        if(istpd  .eq. 2) istpd  = istp
        if(istpe  .eq. 2) istpe  = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istper .eq. 2) istper = istp
!
! copy input flags to a buffer for later use and broadcasting.
!
         buf_in( 1) = qcon
         buf_in( 2) = qlin
         buf_in( 3) = courno
         buf_in( 4) = dtrat
         buf_in( 5) = dfloor
         buf_in( 6) = efloor
         buf_in( 7) = v1floor
         buf_in( 8) = v2floor
         buf_in( 9) = v3floor
         buf_in(10) = b1floor
         buf_in(11) = b2floor
         buf_in(12) = b3floor
         buf_in(13) = emf1floor
         buf_in(14) = emf2floor
         buf_in(15) = emf3floor
         buf_in(16) = erfloor
         buf_in(17) = gpfloor
         buf_in(18) = tlim
         buf_in(19) = cpulim
         buf_in(20) = tsave
!
        ibuf_in( 1) = iordd 
        ibuf_in( 2) = iorde 
        ibuf_in( 3) = iords1
        ibuf_in( 4) = iords2
        ibuf_in( 5) = iords3
        ibuf_in( 6) = iordb1
        ibuf_in( 7) = iordb2
        ibuf_in( 8) = iordb3
        ibuf_in( 9) = iorder
        ibuf_in(10) = istpd 
        ibuf_in(11) = istpe 
        ibuf_in(12) = istps1
        ibuf_in(13) = istps2
        ibuf_in(14) = istps3
        ibuf_in(15) = istpb1
        ibuf_in(16) = istpb2
        ibuf_in(17) = istpb3
        ibuf_in(18) = istper
        ibuf_in(19) = nlim
        ibuf_in(20) = mbatch
        ibuf_in(21) = iCMA
      endif
!
! Broadcast pcon and hycon to the others (use arrays).
!
      call MPI_BCAST( buf_in,20,MPI_DOUBLE_PRECISION  ,0 &
                                             ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibuf_in,21,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (myid_w .ne. 0) then
        qcon    =  buf_in( 1)
        qlin    =  buf_in( 2)
        courno  =  buf_in( 3)
        dtrat   =  buf_in( 4)
        dfloor  =  buf_in( 5)
        efloor  =  buf_in( 6)
        v1floor =  buf_in( 7)
        v2floor =  buf_in( 8)
        v3floor =  buf_in( 9)
        b1floor =  buf_in(10)
        b2floor =  buf_in(11)
        b3floor =  buf_in(12)
        emf1floor=  buf_in(13)
        emf2floor=  buf_in(14)
        emf3floor=  buf_in(15)
        erfloor =  buf_in(16)
        gpfloor =  buf_in(17)
        tlim    =  buf_in(18)
        cpulim  =  buf_in(19)
        tsave   =  buf_in(20)
!
        iordd   = ibuf_in( 1)
        iorde   = ibuf_in( 2)
        iords1  = ibuf_in( 3)
        iords2  = ibuf_in( 4)
        iords3  = ibuf_in( 5)
        iordb1  = ibuf_in( 6)
        iordb2  = ibuf_in( 7)
        iordb3  = ibuf_in( 8)
        iorder  = ibuf_in( 9)
        istpd   = ibuf_in(10)
        istpe   = ibuf_in(11)
        istps1  = ibuf_in(12)
        istps2  = ibuf_in(13)
        istps3  = ibuf_in(14)
        istpb1  = ibuf_in(15)
        istpb2  = ibuf_in(16)
        istpb3  = ibuf_in(17)
        istper  = ibuf_in(18)
        nlim    = ibuf_in(19)
        mbatch  = ibuf_in(20)
        iCMA    = ibuf_in(21)
      endif
!------------------------------ RT CONTROL -----------------------------
!
!  ichem     = activate 9-species chemistry solver (off=0, on=1)--read
!              in by configure.F in the rchmconf namelist, not by setup 
!              in rtcon
!  iRT       = specify which RT scheme to apply (none=0, raytracing=1,
!              photon packet=2)--this is read in by configure.F in the 
!              rchmconf namelist, not by setup in rtcon
!  nnu1,nnu2 = number of photon bins below and above 13.6 eV, respectively.
!              These are also read into the rchmconf namelist rather than
!              in rtcon because they are used to allocate arrays at runtime,
!              as is iRT
!  ifreq     = specify monochromatic (=1), multifrequency (=2), or multigroup
!              (=3) transport 
!  ispct     = type of source spectrum (monochromatic=1, black body=2, QSO=3
!              time-dependent=4 --see spectrum.F). ispct=4 should only be 
!              used if ifreq=2 and for nspec > 6 (full H_2 chemistry); this
!              option requires that photon rates for the freqency bins be
!              read in from external tables (can be spectra for individual 
!              stars or spectral energy distributions (SED's) for stellar
!              populations or galaxies (but the data must be photons/sec for
!              each energy bin)
!  iOTS      = use the on-the-spot approximation (in other words, do not 
!              explicitly transport diffuse (recombination) radiation) (no=0,
!              yes=1)
!  iLW       = include Lyman-Werner radiation (no=0, single band=1, multiline
!              =2)
!  ibkgnd    = include uniform (metagalactic) Lyman-Werner background (no=0, 
!              yes=1).  Only a frequency-averaged flux is currently imple-
!              mented.
!  iHM       = include H- photodetachment radiation (no=0, yes=1)
!  iPWA      = in XYZ or ZRP geometries, attenuate the intensity of the plane 
!              wave according to its distance from its source by 1/r**2 (no=0,
!              yes=1)
!  iextinct  = include dust extinction in photon-conservative transfer (need
!	       user-supplied tables of tau vs hnu)(no=0, yes=1)
!  nsrc      = the number of stars comprising the UV source
!  nphdot    = emission rate of ionizing source photons (s^-1)
!  ephot     = energy (in eV) of source photons, if monochromatic
!  t_on      = time in the problem at which the radiation source should be
!              turned on (in sec)
!  t_off     = time in the problem at which the radiation source should be
!              turned off (in sec)
!  r_sep     = if we are assuming a plane wave UV source, specify r_sep (in
!              pc) to compute photon number flux at far left yz-face (by 
!              nphdot / (4 pi r_sep^2) if XYZ or ZRP are defined together
!              with RT in zeusmp.def)
!  t_star    = effective temperature of black body source (in K)
!  hnumax    = cutoff energy (in eV) for the normalization integral we 
!              compute to obtain an n(nu)dnu that is consistent with nphdot 
!  alpha     = specify QSO power-law spectrum (according to F(nu) being 
!              proportional to nu^(-alpha)
!  nLWdot    = emission rate of photons in the Lyman-Werner band
!  L_star    = log total luminosity (erg/s) of the star (needed to normalize  
!              the H- photodetachment spectrum, obtainable from Schaerer, 2002)
!  J_21      = frequency-averaged flux assumed for the uniform LW background
!              (in units of 10^-21 erg cm^-2 Hz^-1 str^-1, eg Haiman et al 2001)
!              Recall that 4*pi*J_nu = F_nu.
!
        ifreq    = 1
        ispct    = 1
        iOTS     = 1
        iPWA     = 0
        iLW      = 0
        ibkgnd   = 0
        iHM      = 0
        iextinct = 0
        nsrc     = 1
        nphdot   = 1.391d50
        ephot    = 15.0
        t_on     = 0.
        t_off    = 7.95e13
        r_sep    = 272.4
        t_star   = 9.572e04
        hnumax   = 8616.0
        alpha    = -1.8
        nLWdot   = 0.
        L_star   = 6.243
        J_21     = 0
      if (myid_w .eq. 0) then
        read (1,rtcon)
        write(2,rtcon)
        r_sep   = r_sep * cmpc
        ibuf_in( 1) = ifreq
        ibuf_in( 2) = ispct
        ibuf_in( 3) = iOTS
        ibuf_in( 4) = iPWA
        ibuf_in( 5) = iLW
        ibuf_in( 6) = iHM
        ibuf_in( 7) = nsrc
        ibuf_in( 8) = ibkgnd
        ibuf_in( 9) = iextinct
        buf_in ( 1) = nphdot
        buf_in ( 2) = t_off
        buf_in ( 3) = r_sep
        buf_in ( 4) = t_star
        buf_in ( 5) = hnumax
        buf_in ( 6) = alpha
        buf_in ( 7) = nLWdot
        buf_in ( 8) = L_star
        buf_in ( 9) = ephot
        buf_in (10) = t_on
        buf_in (11) = J_21
      endif
!
! Broadcast rtcon to the others.
!
      call MPI_BCAST(ibuf_in, 9,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( buf_in,11,MPI_DOUBLE_PRECISION  , &
                      0,MPI_COMM_WORLD,ierr)
      if (myid_w .ne. 0) then
        ifreq    = ibuf_in( 1)
        ispct    = ibuf_in( 2)
        iOTS     = ibuf_in( 3)
        iPWA     = ibuf_in( 4)
        iLW      = ibuf_in( 5)
        iHM      = ibuf_in( 6)
        nsrc     = ibuf_in( 7)
        ibkgnd   = ibuf_in( 8)
        iextinct = ibuf_in( 9)
        nphdot   = buf_in ( 1)
        t_off    = buf_in ( 2)
        r_sep    = buf_in ( 3)
        t_star   = buf_in ( 4)
        hnumax   = buf_in ( 5)
        alpha    = buf_in ( 6)
        nLWdot   = buf_in ( 7)
        L_star   = buf_in ( 8)
        ephot    = buf_in ( 9)
        t_on     = buf_in (10)
        J_21     = buf_in (11)
      endif
!
!------------------------ CHEMISTRY CONTROL  ---------------------------
!
!  nratec = number of bins in which temperatures are interpolated in
!           the cooling and rate tables computed by CTABLE.  Read in
!           by configure in the rchmconf namelist rather than here in
!           setup because it is used to allocate array sizes at run
!           time
!
!  iceco  = include collisional excitation cooling (yes=1, no=0)
!  icico  = include collisional ionization cooling (yes=1, no=0)
!  ireco  = include recombination cooling          
!  ibrco  = include bremmstrahlung cooling
!  ipiht  = include photoionization heating
!  ih2co  = include ro-vibrational cooling by H2 (0=none, 1=Lepp & Shull, 
!           2=Galli & Palla)
!  icmpt  = include Compton cooling          
!  iDMco  = include Dalgarno McCray (1972) cooling 
!  iHDco  = include HD cooling
!  icycle = choice of timestep hierarchy (1 for original, which is slightly
!           more accurate but much slower at transporting R-type I-fronts or
!           2 for optimized/compressed subcycling, with speedups of 2 - 15
!           for R-type fronts (and is a little faster with D-type fronts too))
!  idust  = include H_2 production by dust grains in reaction network (0=no,
!	    1=yes) (Jura, M. ApJ 1975, 197, 575.)  Gas cooling due to dust
!           is not included in this option.
!  rdshft = redshift (needed to specify Compton cooling)          
!  t_cutoff = temperature below which cooling isn't included in edot
!  fh     = mass fraction of hydrogen in the gas (fh = 1.0 for pure hydrogen
!           and 0.76 if primordial abundances for H and He are assumed)
!  z_sol  = fraction of solar metallicity assumed for gas (set to tiny for
!           zero-metallicity gas to avoid floating point exceptions)
!
!  Note:  electron collisional excitation of H and He is contained in the
!         Dalgarno & McCray cooling curves, so iceco should be set to 0 if
!         iDMco = 1.  Apparently, electron ionizational cooling of H and He
!         are not contained in DM so icico should remain on.
      if (myid_w .eq. 0) then
        iceco  = 0
        icico  = 0
        ireco  = 0
        ibrco  = 0
        ipiht  = 0
        ih2co  = 0
        icmpt  = 0
        iDMco  = 0
        iHDco  = 0
        idust  = 0
        icycle = 1
        rdshft = 0.
        t_cutoff = 0.
        fh     = 1.0
        z_sol  = tiny
        read (1,chemcon)
        write(2,chemcon)
        ibuf_in( 1) = iceco
        ibuf_in( 2) = icico
        ibuf_in( 3) = ireco
        ibuf_in( 4) = ibrco
        ibuf_in( 5) = ipiht
        ibuf_in( 6) = ih2co
        ibuf_in( 7) = icmpt
        ibuf_in( 8) = iDMco
        ibuf_in( 9) = iHDco
        ibuf_in(10) = icycle
        ibuf_in(11) = idust
        buf_in ( 1) = rdshft
        buf_in ( 2) = t_cutoff
        buf_in ( 3) = fh
        buf_in ( 4) = z_sol
      endif
!
! Broadcast chemcon to the others.
!
!      call MPI_BCAST(ibuf_in, 9,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibuf_in,11,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( buf_in,4,MPI_DOUBLE_PRECISION &
      ,0,MPI_COMM_WORLD,ierr)
      if (myid_w .ne. 0) then
!
        iceco  = ibuf_in( 1)
        icico  = ibuf_in( 2)
        ireco  = ibuf_in( 3)
        ibrco  = ibuf_in( 4)
        ipiht  = ibuf_in( 5)
        ih2co  = ibuf_in( 6)
        icmpt  = ibuf_in( 7)
        iDMco  = ibuf_in( 8)
        iHDco  = ibuf_in( 9)
        icycle = ibuf_in(10)
        idust  = ibuf_in(11)
        rdshft   = buf_in ( 1)
        t_cutoff = buf_in ( 2)
        fh       = buf_in ( 3)
        z_sol    = buf_in ( 4)
      endif
!
!  Call ctable to set up chemistry coefficient tables
! 
      if ( xchem )  then
        write(*,*) 'Doing some chemistry'
        call ctable
      endif
!
!------------------------  BOUNDARY CONTROL ----------------------------
!
!  The following points describecoords(1), how boundaries are handled:
!
!  1)  Any of 6 mhd boundary conditions may be specified independently
!  at every zone on the physical problem boundary.  The boundary type is
!  specified by nflo, where
!
!      nflo = 0  =>  interior boundary (get data from neighboring tile)
!           = 1  =>  reflecting (v(normal) = b(normal) = 0)
!           =-1  =>  reflecting (XYZ: same as 1; ZRP: same as 1 with
!                    inversion of 3-components at ijb; RTP: same as 1
!                    with inversion of 2- and 3-components at iib and
!                    inversion of 3-components at ijb and ojb.)
!           = 2  =>  flow out
!           = 3  =>  flow in
!           = 4  =>  periodic
!           = 5  =>  reflecting (v(normal) = 0, b(tangential) = 0)
!
!  Note that in ZRP and RTP, some boundary conditions are implied by
!  the choice of limits.  e.g., if 0 .le. x3a .le. 2*pi in either ZRP
!  or RTP, then periodic boundary conditions should be imposed.
!  Set "niib" to -1 (reflecting with inversion of 2- and 3-components) 
!  if the inner i boundary is at the   origin (RTP).
!  Set "nijb" to -1 (reflecting with inversion of 3-components) if
!     the inner j boundary is on the "Z" axis (ZRP or RTP).
!  Set "nojb" to -1 (reflecting with inversion of 3-components) if
!     the outer j boundary is on the "Z" axis (RTP).
!
!  Since the grid is staggered, the boundary conditions apply over
!  slightly different regions, depending on the centering of the 
!  variable c  and the type of boundary condition used.  Thus, "niib" 
!  is applied to zone or 1-face centered quantities (d, v1, e, gp, b1),
!  "niib2" is applied to 2-face centred quantities (v2, b2), 
!  "niib3" is applied to 3-face centred quantities (v3, b3), and 
!  "niib23" is appled to corner centred quantities (emf1, emf2, emf3).
!  Note that the secondary boundary integer flags are determined from 
!  "niib" automatically.
!
!  Since constant values are often used and are easy to input via 
!  namelists, a single scalar variable is read in for nflo on each 
!  boundary.  These are read in as:
!
!      niis(1),nois(1),nijs(1),nojs(1),niks(1),noks(1)
!
!  For more complicated boundary conditions, values of nflo at each 
!  boundary zone are stored in the 2-D arrays:
!
!      niib,noib,nijb,nojb,nikb,nokb
!
!  These arrays are automatically set to the constant input values.  For
!  more complicated problems, these arrays may be altered in the Problem
!  Generator.
!
!  2) Since the radiation boundary types may differ from the fluid 
!  boundary types, we may specify the former independently of the 
!  latter.  Thus, the radiation boundary types are specified by lflo,
!  where
!
!      lflo = 0  =>  interior boundary (get data from neighboring tile)
!           = 1  =>  reflecting
!           = 2  =>  flow out
!           = 3  =>  flow in
!           = 4  =>  periodic
!
!  Constant values for lflo are read in as:
!
!      niis(2),nois(2),nijs(2),nojs(2),niks(2),noks(2)
!
!  while the 2-D arrays are liib,loib,lijb,lojb,likb,lokb.
!
!  3) Boundary conditions on the gravitational potential are
!     specified by igr, where
!
!       igr = 0  =>  interior boundary (get data from neighboring tile)
!           = 1  =>  reflecting (dgp/d(normal) = 0 "von Neumann")
!           = 2  =>  outflow (equivalent to reflecting)
!           = 3  =>  gp specified (Dirichlet)
!           = 4  =>  periodic
!
!     The flags are defined by analogy with the hydro nflo flag.  This
!     is quite different from ZEUS-2D.  The flags igr are read in as:
!
!      niis(3),nois(3),nijs(3),nojs(3),niks(3),noks(3)
!
!     No 2-D integer arrays for gravity BCs, I guess.
!
!  4) For flow-in boundaries, boundary values of d,e,v1,v2,[v3],...
!  must be input.  Since constant values are often used and are easy to 
!  input via namelists, a single scalar variable is read in for each 
!  field variable on each boundary, for example, the quantities
!
!    fiis(1),fois(1),fijs(1),fojs(1),fiks(1),foks(1) 
!
!  give the constant boundary values for d.  There is a set of f's
!  for e, v, b, er, and gp.
!
!  For more complicated inflow boundary conditions, arrays are used to 
!  store specified values at each boundary zone for each function; for 
!  example:
!
!    diib(j,k,2) is inner i boundary density for sweep j,k at ism2
!    diib(j,k,1) is inner i boundary density for sweep j,k at ism1
!    doib(j,k,1) is outer i boundary density for sweep j,k at iep1
!    doib(j,k,2) is outer i boundary density for sweep j,k at iep2
!
!    dijb(k,i,2) is inner j boundary density for sweep k,i at jsm2
!    dijb(k,i,1) is inner j boundary density for sweep k,i at jsm1
!    dojb(k,i,1) is outer j boundary density for sweep k,i at jep1
!    dojb(k,i,2) is outer j boundary density for sweep k,i at jep2
!
!    dikb(i,j,2) is inner k boundary density for sweep i,j at ksm2
!    dikb(i,j,1) is inner k boundary density for sweep i,j at ksm1
!    dokb(i,j,1) is outer k boundary density for sweep i,j at kep1
!    dokb(i,j,2) is outer k boundary density for sweep i,j at kep2
!
!    These arrays are initialized automatically to constant values;
!    the user's Problem Generator may alter these arrays, e.g. to
!    specify an inlet for a jet or for certain advection tests.
!
! 5) For corner zones [(ii-1,ji-1),(io+1,ji-1),etc] there is a
!    "pecking order" of precedence which attempts to pick the BC
!    that provides the most stable solution.
!    THE USER MAY WANT TO OVERIDE THIS CHOICE in the Problem Generator.
!
!  Set defaults and read in boundary values.  buf_in still contains
!  the "floor" or default values for each field variable.
!
      if (myid_w .eq. 0) then
        do 10 i=1,3
          niis(i) = 0
          nijs(i) = 0
          niks(i) = 0
          nois(i) = 0
          nojs(i) = 0
          noks(i) = 0
10      continue
        do 20 i=1,nbvar
          fiis(i) = buf_in(i+4)
          fois(i) = buf_in(i+4)
          fijs(i) = buf_in(i+4)
          fojs(i) = buf_in(i+4)
          fiks(i) = buf_in(i+4)
          foks(i) = buf_in(i+4)
20      continue
!
        read (1,iib)
        write(2,iib)
!
        read (1,oib)
        write(2,oib)
!
        read (1,ijb)
        write(2,ijb)
!
        read (1,ojb)
        write(2,ojb)
!
        read (1,ikb)
        write(2,ikb)
!
        read (1,okb)
        write(2,okb)
!
      endif
!
!
! Tell the others what the master has read.
!
      call MPI_BCAST(niis, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nois, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nijs, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nojs, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(niks, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(noks, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
      call MPI_BCAST(fiis,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fois,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fijs,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fojs,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fiks,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(foks,nbvar,MPI_DOUBLE_PRECISION  ,0 &
                                          ,MPI_COMM_WORLD,ierr)
!
!-------------------------  SET BOUNDARY TYPES  ------------------------
!
      if (coords(1) .gt. 0) then
        niis(1) = 0
        niis(2) = 0
        niis(3) = 0
      endif
      if (coords(2) .gt. 0) then
        nijs(1) = 0
        nijs(2) = 0
        nijs(3) = 0
      endif
      if (coords(3) .gt. 0) then
        niks(1) = 0
        niks(2) = 0
        niks(3) = 0
      endif
      if (coords(1) .lt. ntiles(1) - 1) then
        nois(1) = 0
        nois(2) = 0
        nois(3) = 0
      endif
      if (coords(2) .lt. ntiles(2) - 1) then
        nojs(1) = 0
        nojs(2) = 0
        nojs(3) = 0
      endif
      if (coords(3) .lt. ntiles(3) - 1) then
        noks(1) = 0
        noks(2) = 0
        noks(3) = 0
      endif
!
!-------------------------  INITIALIZE GRID  ---------------------------
!
!  Routine ggen reads ggen1, ggen2, ggen3; computes the grid for each 
!  tile.
!
      call MPI_BARRIER(comm3d, ierr)
      call ggen
!
!-------------------------  SET MPI DATATYPES  -------------------------
!
! Define MPI Derived Datatypes for passing 2-D slices of 3-D arrays.
!
      call MPI_TYPE_VECTOR (  jn*kn,1 ,in &
      , MPI_DOUBLE_PRECISION, i_slice,ierr)
      call MPI_TYPE_COMMIT (i_slice,ierr)
!
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm ,in*neqm &
      , MPI_DOUBLE_PRECISION, ils_slice,ierr)
      call MPI_TYPE_COMMIT (ils_slice,ierr)
!
      call MPI_TYPE_VECTOR (  jn*kn*nspec,1 ,in &
      , MPI_DOUBLE_PRECISION, iab_slice,ierr)
      call MPI_TYPE_COMMIT (iab_slice,ierr)
!
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm*neqm ,in*neqm*neqm &
      , MPI_DOUBLE_PRECISION, ilsm_slice,ierr)
      call MPI_TYPE_COMMIT (ilsm_slice,ierr)
!
      call MPI_TYPE_VECTOR (    kn  ,in,in*jn,MPI_DOUBLE_PRECISION &
                                       ,j_slice,ierr)
      call MPI_TYPE_COMMIT (j_slice,ierr)
!
      call MPI_TYPE_VECTOR (    kn  ,in*neqm,in*jn*neqm, &
                            MPI_DOUBLE_PRECISION &
                                       ,jls_slice,ierr)
      call MPI_TYPE_COMMIT (jls_slice,ierr)
!
      call MPI_TYPE_VECTOR (    kn*nspec  ,in,in*jn, &
                            MPI_DOUBLE_PRECISION &
                                       ,jab_slice,ierr)
      call MPI_TYPE_COMMIT (jab_slice,ierr)
!
      call MPI_TYPE_VECTOR (    kn  ,in*neqm*neqm,in*jn*neqm*neqm, &
                            MPI_DOUBLE_PRECISION &
                                       ,jlsm_slice,ierr)
      call MPI_TYPE_COMMIT (jlsm_slice,ierr)
!
      call MPI_TYPE_VECTOR (1,in*jn ,1    ,MPI_DOUBLE_PRECISION &
                                       ,k_slice,ierr)
      call MPI_TYPE_COMMIT (k_slice,ierr)
!
      call MPI_TYPE_VECTOR (1,in*jn*neqm ,1    , &
                            MPI_DOUBLE_PRECISION &
                                       ,kls_slice,ierr)
      call MPI_TYPE_COMMIT (kls_slice,ierr)
!
      call MPI_TYPE_VECTOR (nspec,in*jn ,in*jn*kn, &
                            MPI_DOUBLE_PRECISION &
                                       ,kab_slice,ierr)
      call MPI_TYPE_COMMIT (kab_slice,ierr)
!
      call MPI_TYPE_VECTOR (1,in*jn*neqm*neqm ,1    , &
                            MPI_DOUBLE_PRECISION &
                                       ,klsm_slice,ierr)
      call MPI_TYPE_COMMIT (klsm_slice,ierr)
!
!------------------------  GRAVITY CONTROL  ----------------------------
!
!  Self-gravity can be included
!  by solving the Poisson equation in the GRAVITY module.
!  Point mass potentials are included directly in the momentum eqn
!  by using a non-zero value for the variable ptmass.  Point mass
!  potentials do not require defining GRAV, do not call the GRAVITY
!  module, and are not included in the array phi but are explicitly
!  added to the momentum eqn terms in FORCES.
!     g              gravitational constant
!     tgrav          time when gravitation is switched on
!     ptmass         mass of fixed point mass 
!     x1ptm          x1 coordinate of the point mass
!     x2ptm
!     x3ptm
!
!  Boundary flags more general than the ones below are input in the
!  BC namelists described above.  Keep these for now.  Strange, though.
!
!     igrijb = ijb flag (0 for     symmetric (Dirichlet) boundary      )
!                       (1 for non-symmetric             boundary whose
!                        value is calculated using multipole expansion )
!     igrojb = ojb flag ("  "      "          "         "              )
!
!  These are for the multigrid solver.  R: real, I: integer
!
!     gsup       I   Cycle number to calculate gravitational potential
      if (myid_w .eq. 0) then
        guniv     = 6.667d-8
        tgrav     = 0.0
        gsup      = 1
        xwedge    = .false.
!
        read (1,grvcon, iostat=iost)
        write(2,grvcon)
      endif
      igcall = 0
!
!------------------------  RADIATION CONTROL  --------------------------
!
      if(myid_w .eq. 0) then
       ifld  = 1
       epsme = 1.0D-8
       demax = 0.2D0
       dermax = 0.2D0
       nmeiter = 20
       radth = 1.0D0
       epsrad = 1.0D-8
       cnvcrit = 1
       ernom   = 1.0D0
       ennom   = 1.0D0
       epsmaxd = 0.05D0
       cgerrcrit = 666
       ipcflag   = 666
       rmfp0     = huge
       xnu       = huge
       powr      = huge
       rho0      = huge
       t_0       = huge
       read (1,radcon)
       write(2,radcon)
      endif ! myid_w
!
!------------------------  EQUATION OF STATE  --------------------------
!
!      gamma = ratio of specific heats
!      ciso  = isothermal sound speed
!
      if (myid_w .eq. 0) then
        gamma = 0.0
        ciso  = 0.0
        mmw   = 1.0
!
        read (1,eqos)
        write(2,eqos)
      endif
!
! Broadcast grvcon, radcon, and eqos to the others.
!
      if (myid_w .eq. 0) then
         buf_in( 1) = guniv
         buf_in( 2) = tgrav
         buf_in( 3) = ptmass
         buf_in( 4) = x1ptm
         buf_in( 5) = x2ptm
         buf_in( 6) = x3ptm
         buf_in( 7) = epsme
         buf_in( 8) = demax
         buf_in( 9) = dermax
         buf_in(10) = radth
         buf_in(11) = epsrad
         buf_in(12) = epsmaxd
         buf_in(13) = epsmaxd ! there used to be an obsolete parameter here
         buf_in(14) = rmfp0
         buf_in(15) = xnu
         buf_in(16) = powr
         buf_in(17) = rho0
         buf_in(18) = t_0
         buf_in(19) =  gamma
         buf_in(20) =  ciso
         buf_in(21) =  mmw
         buf_in(22) =  ernom
         buf_in(23) =  ennom
        ibuf_in( 1) = gsup
        ibuf_in( 2) = nmeiter
        ibuf_in( 3) = cnvcrit
        ibuf_in( 4) = cgerrcrit
        ibuf_in( 5) = ipcflag
        ibuf_in( 6) = ifld
      endif
      call MPI_BCAST( buf_in,23,MPI_DOUBLE_PRECISION  ,0,comm3d,ierr)
      call MPI_BCAST(ibuf_in,6,MPI_INTEGER,0,comm3d,ierr)
      call MPI_BCAST(xwedge ,1,MPI_LOGICAL,0,comm3d,ierr)
      if (myid_w .ne. 0) then
        guniv      =  buf_in( 1)
        tgrav      =  buf_in( 2)
        ptmass     =  buf_in( 3)
        x1ptm      =  buf_in( 4)
        x2ptm      =  buf_in( 5)
        x3ptm      =  buf_in( 6)
        epsme      =  buf_in( 7)
        demax      =  buf_in( 8)
        dermax     =  buf_in( 9)
        radth      =  buf_in(10)
        epsrad     =  buf_in(11)
        epsmaxd    =  buf_in(12)
        epsmaxd    =  buf_in(13) ! there used to be an obsolete parameter here
        rmfp0      =  buf_in(14)
        xnu        =  buf_in(15)
        powr       =  buf_in(16)
        rho0       =  buf_in(17)
        t_0        =  buf_in(18)
        gamma      =  buf_in(19)
        ciso       =  buf_in(20)
        mmw        =  buf_in(21)
        ernom      =  buf_in(22)
        ennom      =  buf_in(23)
        gsup       = ibuf_in( 1)
        nmeiter    = ibuf_in( 2)
        cnvcrit    = ibuf_in( 3)
        cgerrcrit  = ibuf_in( 4)
        ipcflag    = ibuf_in( 5)
        ifld       = ibuf_in( 6)
      endif
      gamm1 = gamma - 1.0
!
!-----------------------  SET CONSTANT VALUES IN BOUNDARY ARRAYS -------
!
!  For the loops below, set the boundary flags to zero for internal
!  boundaries.  Copy the constant values into the 2-D arrays.
!
!-----------------------  Constant values for IIB  ---------------------
!
       do 40 k=ks-2,ke+3
         do 30 j=js-2,je+3
           if (coords(1) .eq. 0) then
             niib    (j,k) = niis(1)
            if(lrad .ne. 0) then
             liib    (j,k) = niis(2)
            endif
!             if(xmhd) then
              niib2   (j,k) = niis(1)
              niib3   (j,k) = niis(1)
              niib23  (j,k) = niis(1)
!             endif
           else
             niib    (j,k) = 0
            if(lrad .ne. 0) then
             liib    (j,k) = 0
            endif
!             if(xmhd) then
              niib2   (j,k) = 0
              niib3   (j,k) = 0
              niib23  (j,k) = 0
!             endif
           endif
           diib    (j,k,1) = fiis(1)
           diib    (j,k,2) = huge
           diib    (j,k,3) = huge
          if(xiso .eqv. .false.) then
           eiib    (j,k,1) = fiis(2)
           eiib    (j,k,2) = huge
          endif
           v1iib   (j,k,1) = fiis(3)
           v1iib   (j,k,2) = huge
           v2iib   (j,k,1) = fiis(4)
           v2iib   (j,k,2) = huge
           v3iib   (j,k,1) = fiis(5)
           v3iib   (j,k,2) = huge
          if(xmhd) then
           b1iib   (j,k,1) = fiis(6)
           b1iib   (j,k,2) = huge
           b2iib   (j,k,1) = fiis(7)
           b2iib   (j,k,2) = huge
           b3iib   (j,k,1) = fiis(8)
           b3iib   (j,k,2) = huge
           emf1iib (j,k,1) = fiis(9)
           emf1iib (j,k,2) = huge
           emf1iib (j,k,3) = huge
           emf2iib (j,k,1) = fiis(10)
           emf2iib (j,k,2) = huge
           emf2iib (j,k,3) = huge
           emf3iib (j,k,1) = fiis(11)
           emf3iib (j,k,2) = huge
           emf3iib (j,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           eriib   (j,k,1) = fiis(12)
           eriib   (j,k,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpiib   (j,k,1) = fiis(13)
30       continue
40     continue
!
!-----------------------  Constant values for OIB  ---------------------
!
       do 60 k=ks-2,ke+3
         do 50 j=js-2,je+3
           if (coords(1) .eq. ntiles(1) - 1) then
             noib    (j,k) = nois(1)
            if(lrad .ne. 0) then
             loib    (j,k) = nois(2)
            endif
!             if(xmhd) then
              noib2   (j,k) = nois(1)
              noib3   (j,k) = nois(1)
              noib23  (j,k) = nois(1)
!             endif
           else
             noib    (j,k) = 0
            if(lrad .ne. 0) then
             loib    (j,k) = 0
            endif
!             if(xmhd) then
              noib2   (j,k) = 0
              noib3   (j,k) = 0
              noib23  (j,k) = 0
!             endif
           endif
           doib    (j,k,1) = fois(1)
           doib    (j,k,2) = huge
           doib    (j,k,3) = huge
          if(xiso .eqv. .false.) then
           eoib    (j,k,1) = fois(2)
           eoib    (j,k,2) = huge
          endif ! xiso
           v1oib   (j,k,1) = fois(3)
           v1oib   (j,k,2) = huge
           v2oib   (j,k,1) = fois(4)
           v2oib   (j,k,2) = huge
           v3oib   (j,k,1) = fois(5)
           v3oib   (j,k,2) = huge
          if(xmhd) then
           b1oib   (j,k,1) = fois(6)
           b1oib   (j,k,2) = huge
           b2oib   (j,k,1) = fois(7)
           b2oib   (j,k,2) = huge
           b3oib   (j,k,1) = fois(8)
           b3oib   (j,k,2) = huge
           emf1oib (j,k,1) = fois(9)
           emf1oib (j,k,2) = huge
           emf1oib (j,k,3) = huge
           emf2oib (j,k,1) = fois(10)
           emf2oib (j,k,2) = huge
           emf2oib (j,k,3) = huge
           emf3oib (j,k,1) = fois(11)
           emf3oib (j,k,2) = huge
           emf3oib (j,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           eroib   (j,k,1) = fois(12)
           eroib   (j,k,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpoib   (j,k,1) = fois(13)
50       continue
60     continue
!
!-----------------------  Constant values for IJB  ---------------------
!
       do 80 k=ks-2,ke+3
         do 70 i=is-2,ie+3
           if (coords(2) .eq. 0) then
             nijb    (i,k) = nijs(1)
            if(lrad .ne. 0) then
             lijb    (i,k) = nijs(2)
            endif
!             if(xmhd) then
              nijb3   (i,k) = nijs(1)
              nijb1   (i,k) = nijs(1)
              nijb31  (i,k) = nijs(1)
!             endif
           else
             nijb    (i,k) = 0
            if(lrad .ne. 0) then
             lijb    (i,k) = 0
            endif
!             if(xmhd) then
              nijb3   (i,k) = 0
              nijb1   (i,k) = 0
              nijb31  (i,k) = 0
!             endif
           endif
           dijb    (i,k,1) = fijs(1)
           dijb    (i,k,2) = huge
           dijb    (i,k,3) = huge
          if(xiso .eqv. .false.) then
           eijb    (i,k,1) = fijs(2)
           eijb    (i,k,2) = huge
          endif ! xiso
           v1ijb   (i,k,1) = fijs(3)
           v1ijb   (i,k,2) = huge
           v2ijb   (i,k,1) = fijs(4)
           v2ijb   (i,k,2) = huge
           v3ijb   (i,k,1) = fijs(5)
           v3ijb   (i,k,2) = huge
          if(xmhd) then
           b1ijb   (i,k,1) = fijs(6)
           b1ijb   (i,k,2) = huge
           b2ijb   (i,k,1) = fijs(7)
           b2ijb   (i,k,2) = huge
           b3ijb   (i,k,1) = fijs(8)
           b3ijb   (i,k,2) = huge
           emf1ijb (i,k,1) = fijs(9)
           emf1ijb (i,k,2) = huge
           emf1ijb (i,k,3) = huge
           emf2ijb (i,k,1) = fijs(10)
           emf2ijb (i,k,2) = huge
           emf2ijb (i,k,3) = huge
           emf3ijb (i,k,1) = fijs(11)
           emf3ijb (i,k,2) = huge
           emf3ijb (i,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erijb   (i,k,1) = fijs(12)
           erijb   (i,k,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpijb   (i,k,1) = fijs(13)
70       continue
80     continue
!
!-----------------------  Constant values for OJB  ---------------------
!
       do 100 k=ks-2,ke+3
         do 90 i=is-2,ie+3
           if (coords(2) .eq. ntiles(2) - 1) then
             nojb    (i,k) = nojs(1)
            if(lrad .ne. 0) then
             lojb    (i,k) = nojs(2)
            endif
!             if(xmhd) then
              nojb3   (i,k) = nojs(1)
              nojb1   (i,k) = nojs(1)
              nojb31  (i,k) = nojs(1)
!             endif
           else
             nojb    (i,k) = 0
            if(lrad .ne. 0) then
             lojb    (i,k) = 0
            endif
!             if(xmhd) then
              nojb3   (i,k) = 0
              nojb1   (i,k) = 0
              nojb31  (i,k) = 0
!             endif
           endif
           dojb    (i,k,1) = fojs(1)
           dojb    (i,k,2) = huge
           dojb    (i,k,3) = huge
          if(xiso .eqv. .false.) then
           eojb    (i,k,1) = fojs(2)
           eojb    (i,k,2) = huge
          endif ! xiso
           v1ojb   (i,k,1) = fojs(3)
           v1ojb   (i,k,2) = huge
           v2ojb   (i,k,1) = fojs(4)
           v2ojb   (i,k,2) = huge
           v3ojb   (i,k,1) = fojs(5)
           v3ojb   (i,k,2) = huge
          if(xmhd) then
           b1ojb   (i,k,1) = fojs(6)
           b1ojb   (i,k,2) = huge
           b2ojb   (i,k,1) = fojs(7)
           b2ojb   (i,k,2) = huge
           b3ojb   (i,k,1) = fojs(8)
           b3ojb   (i,k,2) = huge
           emf1ojb (i,k,1) = fojs(9)
           emf1ojb (i,k,2) = huge
           emf1ojb (i,k,3) = huge
           emf2ojb (i,k,1) = fojs(10)
           emf2ojb (i,k,2) = huge
           emf2ojb (i,k,3) = huge
           emf3ojb (i,k,1) = fojs(11)
           emf3ojb (i,k,2) = huge
           emf3ojb (i,k,3) = huge
          endif
          if(lrad .ne. 0) then
           erojb   (i,k,1) = fojs(12)
           erojb   (i,k,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpojb   (i,k,1) = fojs(13)
90       continue
100    continue
!
!-----------------------  Constant values for IKB  ---------------------
!
       do 120 j=js-2,je+3
         do 110 i=is-2,ie+3
           if (coords(3) .eq. 0) then
             nikb    (i,j) = niks(1)
            if(lrad .ne. 0) then
             likb    (i,j) = niks(2)
            endif
!             if(xmhd) then
              nikb1   (i,j) = niks(1)
              nikb2   (i,j) = niks(1)
              nikb12  (i,j) = niks(1)
!             endif
           else
             nikb    (i,j) = 0
            if(lrad .ne. 0) then
             likb    (i,j) = 0
            endif
!            if(xmhd) then
             nikb1   (i,j) = 0
             nikb2   (i,j) = 0
             nikb12  (i,j) = 0
!            endif
           endif
           dikb    (i,j,1) = fiks(1)
           dikb    (i,j,2) = huge
           dikb    (i,j,3) = huge
          if(xiso .eqv. .false.) then
           eikb    (i,j,1) = fiks(2)
           eikb    (i,j,2) = huge
          endif
           v1ikb   (i,j,1) = fiks(3)
           v1ikb   (i,j,2) = huge
           v2ikb   (i,j,1) = fiks(4)
           v2ikb   (i,j,2) = huge
           v3ikb   (i,j,1) = fiks(5)
           v3ikb   (i,j,2) = huge
          if(xmhd) then
           b1ikb   (i,j,1) = fiks(6)
           b1ikb   (i,j,2) = huge
           b2ikb   (i,j,1) = fiks(7)
           b2ikb   (i,j,2) = huge
           b3ikb   (i,j,1) = fiks(8)
           b3ikb   (i,j,2) = huge
           emf1ikb (i,j,1) = fiks(9)
           emf1ikb (i,j,2) = huge
           emf1ikb (i,j,3) = huge
           emf2ikb (i,j,1) = fiks(10)
           emf2ikb (i,j,2) = huge
           emf2ikb (i,j,3) = huge
           emf3ikb (i,j,1) = fiks(11)
           emf3ikb (i,j,2) = huge
           emf3ikb (i,j,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erikb   (i,j,1) = fiks(12)
           erikb   (i,j,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpikb   (i,j,1) = fiks(13)
110      continue
120    continue
!
!-----------------------  Constant values for OKB  ---------------------
!
       do 140 j=js-2,je+3
         do 130 i=is-2,ie+3
           if (coords(3) .eq. ntiles(3) - 1) then
            nokb    (i,j) = noks(1)
            if(lrad .ne. 0) then
             lokb    (i,j) = noks(2)
            endif
!            if(xmhd) then
             nokb1   (i,j) = noks(1)
             nokb2   (i,j) = noks(1)
             nokb12  (i,j) = noks(1)
!            endif
           else
            nokb    (i,j) = 0
            if(lrad .ne. 0) then
             lokb    (i,j) = 0
            endif
!            if(xmhd) then
             nokb1   (i,j) = 0
             nokb2   (i,j) = 0
             nokb12  (i,j) = 0
!            endif
           endif
           dokb    (i,j,1) = foks(1)
           dokb    (i,j,2) = huge
           dokb    (i,j,3) = huge
          if(xiso .eqv. .false.) then
           eokb    (i,j,1) = foks(2)
           eokb    (i,j,2) = huge
          endif
           v1okb   (i,j,1) = foks(3)
           v1okb   (i,j,2) = huge
           v2okb   (i,j,1) = foks(4)
           v2okb   (i,j,2) = huge
           v3okb   (i,j,1) = foks(5)
           v3okb   (i,j,2) = huge
          if(xmhd) then
           b1okb   (i,j,1) = foks(6)
           b1okb   (i,j,2) = huge
           b2okb   (i,j,1) = foks(7)
           b2okb   (i,j,2) = huge
           b3okb   (i,j,1) = foks(8)
           b3okb   (i,j,2) = huge
           emf1okb (i,j,1) = foks(9)
           emf1okb (i,j,2) = huge
           emf1okb (i,j,3) = huge
           emf2okb (i,j,1) = foks(10)
           emf2okb (i,j,2) = huge
           emf2okb (i,j,3) = huge
           emf3okb (i,j,1) = foks(11)
           emf3okb (i,j,2) = huge
           emf3okb (i,j,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erokb   (i,j,1) = foks(12)
           erokb   (i,j,2) = huge
          endif
          if(xgrav .or. xgrvfft) &
           gpokb   (i,j,1) = foks(13)
130      continue
140    continue
!
!-------------------------  benchmark GENERATOR  -------------------------
!
!  benchmark is a user-defined cpp macro name representing a subroutine
!  which intializes all field variables for the particular problem to
!  be studied.  benchmark should initialize the field variable arrays 
!  for both active zones and at least the first layer of boundary zones,
!  unless the default or input constant values already specify the 
!  desired problem.
!
!  For non-uniform initial magnetic field configurations, to satisfy 
!  the constraint DIV(B)=0 benchmark should initialize
!  b1, b2, b3 by differencing a vector potential.
!
!  First initialize all field variables to default (input) values.
!
      do 170 k=ks-2,ke+3
        do 160 j=js-2,je+3
          do 150 i=is-2,ie+3
            d  (i,j,k) = dfloor
            if(xiso .eqv. .false.) e  (i,j,k) = efloor
            v1 (i,j,k) = v1floor
            v2 (i,j,k) = v2floor
            v3 (i,j,k) = v3floor
           if(xmhd) then
            b1 (i,j,k) = b1floor * g2bi(i) * g31bi(i)
            b2 (i,j,k) = b2floor * g32bi(j)
            b3 (i,j,k) = b3floor * g2bi(i)
           endif ! xmhd
           if(lrad .ne. 0) then
            er (i,j,k) = erfloor
           endif
150       continue
160     continue
170   continue
!
! Set the bvstat array to 0, indicating that the boundary values
! will need to be updated.
!
      do 175 j=1,nbvar
        do 175 i=1,6
          bvstat(i,j) = 0
175   continue
      dt     = huge
!
      call benchmark
!
! Note that the 2-D boundary arrays were set to constant values above,
! which should be overwritten with any desired spatially-dependent 
! functions in the Problem Generator.
!
! Set the scalar boundary value flags niis, etc. to nonzero values if 
! the corresponding boundaries are physical ones.
!
       if (niib(js,ks).ne.0 .and. coords(1).eq.0        ) &
                                    niis(1) = niib(js,ks)
       if (noib(js,ks).ne.0 .and. coords(1).eq.ntiles(1) - 1) &
                                    nois(1) = noib(js,ks)
       if (nijb(is,ks).ne.0 .and. coords(2).eq.0        ) &
                                    nijs(1) = nijb(is,ks)
       if (nojb(is,ks).ne.0 .and. coords(2).eq.ntiles(2) - 1) &
                                    nojs(1) = nojb(is,ks)
       if (nikb(is,js).ne.0 .and. coords(3).eq.0        ) &
                                    niks(1) = nikb(is,js)
       if (nokb(is,js).ne.0 .and. coords(3).eq.ntiles(3) - 1) &
                                    noks(1) = nokb(is,js)
!
!-----------------------  INITIALIZE ADDITIONAL BOUNDARY ARRAYS --------
!
!      Set the secondary boundary flags ("niib2", "niib3", "niib23",
!  etc.) from the primary boundary flags ("niib", etc.).
!
       call bndyflgs
!
!      If the second (and third) layer boundaries were not set by the
!  user, set them equal to the first layer boundary values.
!
       do 190 k=ks-2,ke+3
         do 180 j=js-2,je+3
           if (   diib(j,k,2) .eq. huge)    diib(j,k,2) =    diib(j,k,1)
           if (   diib(j,k,3) .eq. huge)    diib(j,k,3) =    diib(j,k,2)
          if(xiso .eqv. .false.) then
           if (   eiib(j,k,2) .eq. huge)    eiib(j,k,2) =    eiib(j,k,1)
          endif ! xiso
           if (  v1iib(j,k,2) .eq. huge)   v1iib(j,k,2) =   v1iib(j,k,1)
           if (  v2iib(j,k,2) .eq. huge)   v2iib(j,k,2) =   v2iib(j,k,1)
           if (  v3iib(j,k,2) .eq. huge)   v3iib(j,k,2) =   v3iib(j,k,1)
          if(xmhd) then
           if (  b1iib(j,k,2) .eq. huge)   b1iib(j,k,2) =   b1iib(j,k,1)
           if (  b2iib(j,k,2) .eq. huge)   b2iib(j,k,2) =   b2iib(j,k,1)
           if (  b3iib(j,k,2) .eq. huge)   b3iib(j,k,2) =   b3iib(j,k,1)
           if (emf1iib(j,k,2) .eq. huge) emf1iib(j,k,2) = emf1iib(j,k,1)
           if (emf1iib(j,k,3) .eq. huge) emf1iib(j,k,3) = emf1iib(j,k,2)
           if (emf2iib(j,k,2) .eq. huge) emf2iib(j,k,2) = emf2iib(j,k,1)
           if (emf2iib(j,k,3) .eq. huge) emf2iib(j,k,3) = emf2iib(j,k,2)
           if (emf3iib(j,k,2) .eq. huge) emf3iib(j,k,2) = emf3iib(j,k,1)
           if (emf3iib(j,k,3) .eq. huge) emf3iib(j,k,3) = emf3iib(j,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  eriib(j,k,2) .eq. huge)   eriib(j,k,2) =   eriib(j,k,1)
          endif
           if (   doib(j,k,2) .eq. huge)    doib(j,k,2) =    doib(j,k,1)
           if (   doib(j,k,3) .eq. huge)    doib(j,k,3) =    doib(j,k,2)
          if(xiso .eqv. .false.) then
           if (   eoib(j,k,2) .eq. huge)    eoib(j,k,2) =    eoib(j,k,1)
          endif
           if (  v1oib(j,k,2) .eq. huge)   v1oib(j,k,2) =   v1oib(j,k,1)
           if (  v2oib(j,k,2) .eq. huge)   v2oib(j,k,2) =   v2oib(j,k,1)
           if (  v3oib(j,k,2) .eq. huge)   v3oib(j,k,2) =   v3oib(j,k,1)
          if(xmhd) then
           if (  b1oib(j,k,2) .eq. huge)   b1oib(j,k,2) =   b1oib(j,k,1)
           if (  b2oib(j,k,2) .eq. huge)   b2oib(j,k,2) =   b2oib(j,k,1)
           if (  b3oib(j,k,2) .eq. huge)   b3oib(j,k,2) =   b3oib(j,k,1)
           if (emf1oib(j,k,2) .eq. huge) emf1oib(j,k,2) = emf1oib(j,k,1)
           if (emf1oib(j,k,3) .eq. huge) emf1oib(j,k,3) = emf1oib(j,k,2)
           if (emf2oib(j,k,2) .eq. huge) emf2oib(j,k,2) = emf2oib(j,k,1)
           if (emf2oib(j,k,3) .eq. huge) emf2oib(j,k,3) = emf2oib(j,k,2)
           if (emf3oib(j,k,2) .eq. huge) emf3oib(j,k,2) = emf3oib(j,k,1)
           if (emf3oib(j,k,3) .eq. huge) emf3oib(j,k,3) = emf3oib(j,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  eroib(j,k,2) .eq. huge)   eroib(j,k,2) =   eroib(j,k,1)
          endif
180      continue
190    continue
!
       do 210 k=ks-2,ke+3
         do 200 i=is-2,ie+3
           if (   dijb(i,k,2) .eq. huge)    dijb(i,k,2) =    dijb(i,k,1)
           if (   dijb(i,k,3) .eq. huge)    dijb(i,k,3) =    dijb(i,k,2)
          if(xiso .eqv. .false.) then
           if (   eijb(i,k,2) .eq. huge)    eijb(i,k,2) =    eijb(i,k,1)
          endif
           if (  v1ijb(i,k,2) .eq. huge)   v1ijb(i,k,2) =   v1ijb(i,k,1)
           if (  v2ijb(i,k,2) .eq. huge)   v2ijb(i,k,2) =   v2ijb(i,k,1)
           if (  v3ijb(i,k,2) .eq. huge)   v3ijb(i,k,2) =   v3ijb(i,k,1)
          if(xmhd) then
           if (  b1ijb(i,k,2) .eq. huge)   b1ijb(i,k,2) =   b1ijb(i,k,1)
           if (  b2ijb(i,k,2) .eq. huge)   b2ijb(i,k,2) =   b2ijb(i,k,1)
           if (  b3ijb(i,k,2) .eq. huge)   b3ijb(i,k,2) =   b3ijb(i,k,1)
           if (emf1ijb(i,k,2) .eq. huge) emf1ijb(i,k,2) = emf1ijb(i,k,1)
           if (emf1ijb(i,k,3) .eq. huge) emf1ijb(i,k,3) = emf1ijb(i,k,2)
           if (emf2ijb(i,k,2) .eq. huge) emf2ijb(i,k,2) = emf2ijb(i,k,1)
           if (emf2ijb(i,k,3) .eq. huge) emf2ijb(i,k,3) = emf2ijb(i,k,2)
           if (emf3ijb(i,k,2) .eq. huge) emf3ijb(i,k,2) = emf3ijb(i,k,1)
           if (emf3ijb(i,k,3) .eq. huge) emf3ijb(i,k,3) = emf3ijb(i,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erijb(i,k,2) .eq. huge)   erijb(i,k,2) =   erijb(i,k,1)
          endif
           if (   dojb(i,k,2) .eq. huge)    dojb(i,k,2) =    dojb(i,k,1)
           if (   dojb(i,k,3) .eq. huge)    dojb(i,k,3) =    dojb(i,k,2)
          if(xiso .eqv. .false.) then
           if (   eojb(i,k,2) .eq. huge)    eojb(i,k,2) =    eojb(i,k,1)
          endif
           if (  v1ojb(i,k,2) .eq. huge)   v1ojb(i,k,2) =   v1ojb(i,k,1)
           if (  v2ojb(i,k,2) .eq. huge)   v2ojb(i,k,2) =   v2ojb(i,k,1)
           if (  v3ojb(i,k,2) .eq. huge)   v3ojb(i,k,2) =   v3ojb(i,k,1)
          if(xmhd) then
           if (  b1ojb(i,k,2) .eq. huge)   b1ojb(i,k,2) =   b1ojb(i,k,1)
           if (  b2ojb(i,k,2) .eq. huge)   b2ojb(i,k,2) =   b2ojb(i,k,1)
           if (  b3ojb(i,k,2) .eq. huge)   b3ojb(i,k,2) =   b3ojb(i,k,1)
           if (emf1ojb(i,k,2) .eq. huge) emf1ojb(i,k,2) = emf1ojb(i,k,1)
           if (emf1ojb(i,k,3) .eq. huge) emf1ojb(i,k,3) = emf1ojb(i,k,2)
           if (emf2ojb(i,k,2) .eq. huge) emf2ojb(i,k,2) = emf2ojb(i,k,1)
           if (emf2ojb(i,k,3) .eq. huge) emf2ojb(i,k,3) = emf2ojb(i,k,2)
           if (emf3ojb(i,k,2) .eq. huge) emf3ojb(i,k,2) = emf3ojb(i,k,1)
           if (emf3ojb(i,k,3) .eq. huge) emf3ojb(i,k,3) = emf3ojb(i,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erojb(i,k,2) .eq. huge)   erojb(i,k,2) =   erojb(i,k,1)
          endif
200      continue
210    continue
!
       do 230 j=js-2,je+3
         do 220 i=is-2,ie+3
           if (   dikb(i,j,2) .eq. huge)    dikb(i,j,2) =    dikb(i,j,1)
           if (   dikb(i,j,3) .eq. huge)    dikb(i,j,3) =    dikb(i,j,2)
          if(xiso .eqv. .false.) then
           if (   eikb(i,j,2) .eq. huge)    eikb(i,j,2) =    eikb(i,j,1)
          endif
           if (  v1ikb(i,j,2) .eq. huge)   v1ikb(i,j,2) =   v1ikb(i,j,1)
           if (  v2ikb(i,j,2) .eq. huge)   v2ikb(i,j,2) =   v2ikb(i,j,1)
           if (  v3ikb(i,j,2) .eq. huge)   v3ikb(i,j,2) =   v3ikb(i,j,1)
          if(xmhd) then
           if (  b1ikb(i,j,2) .eq. huge)   b1ikb(i,j,2) =   b1ikb(i,j,1)
           if (  b2ikb(i,j,2) .eq. huge)   b2ikb(i,j,2) =   b2ikb(i,j,1)
           if (  b3ikb(i,j,2) .eq. huge)   b3ikb(i,j,2) =   b3ikb(i,j,1)
           if (emf1ikb(i,j,2) .eq. huge) emf1ikb(i,j,2) = emf1ikb(i,j,1)
           if (emf1ikb(i,j,3) .eq. huge) emf1ikb(i,j,3) = emf1ikb(i,j,2)
           if (emf2ikb(i,j,2) .eq. huge) emf2ikb(i,j,2) = emf2ikb(i,j,1)
           if (emf2ikb(i,j,3) .eq. huge) emf2ikb(i,j,3) = emf2ikb(i,j,2)
           if (emf3ikb(i,j,2) .eq. huge) emf3ikb(i,j,2) = emf3ikb(i,j,1)
           if (emf3ikb(i,j,3) .eq. huge) emf3ikb(i,j,3) = emf3ikb(i,j,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erikb(i,j,2) .eq. huge)   erikb(i,j,2) =   erikb(i,j,1)
          endif
           if (   dokb(i,j,2) .eq. huge)    dokb(i,j,2) =    dokb(i,j,1)
           if (   dokb(i,j,3) .eq. huge)    dokb(i,j,3) =    dokb(i,j,2)
          if(xiso .eqv. .false.) then
           if (   eokb(i,j,2) .eq. huge)    eokb(i,j,2) =    eokb(i,j,1)
          endif
           if (  v1okb(i,j,2) .eq. huge)   v1okb(i,j,2) =   v1okb(i,j,1)
           if (  v2okb(i,j,2) .eq. huge)   v2okb(i,j,2) =   v2okb(i,j,1)
           if (  v3okb(i,j,2) .eq. huge)   v3okb(i,j,2) =   v3okb(i,j,1)
          if(xmhd) then
           if (  b1okb(i,j,2) .eq. huge)   b1okb(i,j,2) =   b1okb(i,j,1)
           if (  b2okb(i,j,2) .eq. huge)   b2okb(i,j,2) =   b2okb(i,j,1)
           if (  b3okb(i,j,2) .eq. huge)   b3okb(i,j,2) =   b3okb(i,j,1)
           if (emf1okb(i,j,2) .eq. huge) emf1okb(i,j,2) = emf1okb(i,j,1)
           if (emf1okb(i,j,3) .eq. huge) emf1okb(i,j,3) = emf1okb(i,j,2)
           if (emf2okb(i,j,2) .eq. huge) emf2okb(i,j,2) = emf2okb(i,j,1)
           if (emf2okb(i,j,3) .eq. huge) emf2okb(i,j,3) = emf2okb(i,j,2)
           if (emf3okb(i,j,2) .eq. huge) emf3okb(i,j,2) = emf3okb(i,j,1)
           if (emf3okb(i,j,3) .eq. huge) emf3okb(i,j,3) = emf3okb(i,j,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erokb(i,j,2) .eq. huge)   erokb(i,j,2) =   erokb(i,j,1)
          endif
220      continue
230    continue
!
!-------------------------  GRID MOTION CONTROL  -----------------------
!
!  x1fac     x1 motion factor
!            < 0 gives "Lagrangian" tracking in x1 lines
!  x2fac     x2 motion factor
!            < 0 gives "Lagrangian" tracking in x2 lines
!  x3fac     x3 motion factor
!            < 0 gives "Lagrangian" tracking in x3 lines
!  iga       i<ia => zone ratio is preserved in x1 lines
!  jga       j<ja => zone ratio is preserved in x2 lines
!  kga       k<ka => zone ratio is preserved in x3 lines
!  igcon     selects grid treatment:
!            =0 => separate motion
!            =1 => averaged motion
!            =2 => tracking x1, x2, and x3 boundaries
!            =3 => averaged boundary tracking
!            =4 => input grid boundary speeds
!                  vg1(io) = x1fac * central sound speed
!                  vg2(jo) = x2fac * central sound speed
!                  vg3(ko) = x3fac * central sound speed
!
      if (myid_w .eq. 0) then
        x1fac = -1.0
        x2fac = 0.0
        x3fac = 0.0
        iga    = is
        jga    = js
        kga    = ks
        igcon = 5
!
!        read (1,gcon)
!        write(2,gcon)
!
        iga =  min (  max ( iga, is ), ie )
        jga =  min (  max ( jga, js ), je )
        kga =  min (  max ( kga, ks ), ke )
!
      endif
      if (myid_w .eq. 0) then
         buf_in(1) = x1fac
         buf_in(2) = x2fac
         buf_in(3) = x3fac
        ibuf_in(1) = iga
        ibuf_in(2) = jga
        ibuf_in(3) = kga
        ibuf_in(4) = igcon
      endif
      call MPI_BCAST( buf_in,3,MPI_DOUBLE_PRECISION  ,0,comm3d,ierr)
      call MPI_BCAST(ibuf_in,4,MPI_INTEGER,0,comm3d,ierr)
      if (myid_w .ne. 0) then
        x1fac =  buf_in(1)
        x2fac =  buf_in(2)
        x3fac =  buf_in(3)
        iga    = ibuf_in(1)
        jga    = ibuf_in(2)
        kga    = ibuf_in(3)
        igcon = ibuf_in(4)
      endif
!      if(xvgrid) then
       do 320 i=is-2,ie+3
         vg1(i) = 0.0
320    continue
       do 330 j=js-2,je+3
         vg2(j) = 0.0
330    continue
       do 340 k=ks-2,ke+3
         vg3(k) = 0.0
340    continue  
!
!----------------------  INITIAL TIMESTEP  -----------------------------
!
       nreq = 0
!
! Compute the momentum densities from the initial conditions.  We
! need some density boundary values, but we won't try to be efficient
! here.  Be sure to call bvald in sequence so that the values on edges
! get set properly.  
!
       call MPI_BARRIER(comm3d, ierr)
       nreq = 0
       nsub = nsub + 1
       call bvald(1,0,0,0,0,0,d)
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
       nreq = 0
       nsub = nsub + 1
       call bvald(0,0,1,0,0,0,d)
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
       nreq = 0
       nsub = nsub + 1
       call bvald(0,0,0,0,1,0,d)
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
       if(ldimen .eq. 3) then
        kone = 1
       else
        kone = 0
       endif
       if(ldimen .gt. 1) then
        jone = 1
       else
        jone = 0
       endif
!
       do 370 k=ks,ke
        do 360 j=js,je
         do 350 i=is,ie
          w3da(i,j,k) = v1(i,j,k)*0.5*(d(i-1,j     ,k     )+d(i,j,k))
          w3db(i,j,k) = v2(i,j,k)*0.5*(d(i  ,j-jone,k     )+d(i,j,k)) &
                      * g2b(i)
          w3dc(i,j,k) = v3(i,j,k)*0.5*(d(i  ,j     ,k-kone)+d(i,j,k)) &
                      * g31b(i) * g32b(j)
350      continue
360     continue
370    continue
!
!      initialize radiation scratch array
!
      if(lrad .ne. 0) then
       do k = ks, ke
        do j = js, je
         do i = is, ie
          w3dh(i,j,k) = er(i,j,k)/d(i,j,k)
         enddo
        enddo
       enddo
      endif ! lrad
      if(nspec .gt. 1) then
       do n = 1, nspec
        do k = ks, ke
         do j = js, je
          do i = is, ie
           w4da(i,j,k,n) = abun(i,j,k,n)
          enddo
         enddo
        enddo
       enddo
      endif ! nspec
!
! Set the artificial viscosity time step dtqqi2 to zero.
!
       dtqqi2 = 0.0
       dtmin  = tiny
!
!
      if(lrad .ne. 0) then
       dtimrdi2 = 1.0/(dtrat * tlim)**2
       dtnri2   = 1.0 / dt**2
      endif
       call nudt
       dtmin = dtrat * dt
       dt    = 10.0 * dtmin
      if(lrad .ne. 0) then
       dtnri2   = 1.0 / dt**2
       dtrdi2   = 1.0 / dt**2
       dtimrdi2   = 1.0 / dt**2
      endif
!
!-----------------------  INITIALIZE NEW GRID  -------------------------
!
!  Compute grid velocities, and new grid positions [in routine newgrid].
!  Note newgrid will recompute "n+1/2" and "n+1" grid lines, but only
!  if x1fac or x2fac .ne. 0.  The "n+1/2" and "n+1" grid lines are
!  initialized in ggen to the old values in case the grid never moves.  Note
!  newgrid must be called after nudt since it needs the timestep.  Thus,
!  the initial timestep computed above did not account for grid motion.
!
      if(xvgrid) call newgrid
!
!------------------------  Initialize everything else  -----------------
!
      ix1x2x3 = 1
      jx1x2x3 = 1  !  For rad_pf* routines
!
      return
      end
