c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine restart
c
c  PURPOSE:  Sets up restarted run
c
c  EXTERNALS:
c     PROBRES -- macroname which is defined to be the user-supplied
c                subroutine name which re-initializes user-defined
c                variables for the problem to be restarted.
c                PROBRES is undefined by default for backward 
c                compatiability with earlier versions of ZEUS-2D.
c     nudt    -- computes timestep
c
c  LOCALS:
c
c  LAST MODIFIED: Jan. 7, 1999 by JCH for ZEUS-MP
c  LAST MODIFIED: Mar     1999 by efh, including debug
c
c               : DJW 11-28-06 for 9-species primordial chemistry/RT
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use bndry
      use scratch
      use gravmod
      use radiation
      use impsoln
      use opac
      use opac_law
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use chem
c
      implicit NONE
c
      integer  :: i   , j, k, n, jone, kone
      integer  :: iord, istp
      real(rl) :: dtrat
c
      integer  :: nbl, igrid
      real(rl) :: x1min, x1max, x1rat, dx1min,
     .            x2min, x2max, x2rat, dx2min,
     .            x3min, x3max, x3rat, dx3min
c
      logical  lgrid
c
      namelist /pcon/ nlim,tlim,cpulim,tsave,mbatch
      namelist /hycon/
     . qcon,qlin,courno,dtrat,iord,istp,iCMA,
     . iordd,iorde,iords1,iords2,iords3,iordb1,iordb2,iordb3,iorder,
     . istpd,istpe,istps1,istps2,istps3,istpb1,istpb2,istpb3,istper,
     . dfloor,efloor,v1floor,v2floor,v3floor,b1floor,b2floor,b3floor,
     . emf1floor,emf2floor,emf3floor,erfloor,gpfloor
c
      namelist /ggen1/ nbl,x1min,x1max,igrid,x1rat,dx1min,lgrid
      namelist /ggen2/ nbl,x2min,x2max,igrid,x2rat,dx2min,lgrid
      namelist /ggen3/ nbl,x3min,x3max,igrid,x3rat,lgrid
c
      namelist /rtcon/ ifreq, ispct, iOTS, iPWA, iLW, IHM, ibkgnd,
     .                 iextinct, nphdot, nsrc, ephot, t_on, t_off, 
     .                 r_sep, t_star, hnumax, alpha, nLWdot, L_star,
     .                 J_21
c
      namelist /chemcon/ iceco, icico, ireco, ibrco, ipiht, ih2co,
     .                   icmpt, iDMco, iHDco, icycle, t_cutoff, rdshft,
     .                   fh   , idust, z_sol
c
      namelist /iib/     niis, fiis
      namelist /oib/     nois, fois
      namelist /ijb/     nijs, fijs
      namelist /ojb/     nojs, fojs
      namelist /ikb/     niks, fiks
      namelist /okb/     noks, foks
c
      namelist /mpitop/ ntiles, periodic
      namelist /grvcon/tgrav,ptmass,x1ptm,x2ptm,x3ptm,gsup,guniv
      namelist /radcon/ ifld,epsme,demax,dermax,nmeiter,radth,epsrad,
     .     cnvcrit,ernom,ennom,epsmaxd,cgerrcrit,
     .     ipcflag,rmfp0,xnu,powr,rho0,t_0,kpfrac
      namelist /eqos/ gamma, ciso, mmw
      namelist /gcon/ x1fac,x2fac,x3fac,iga,jga,kga,igcon
c
c=======================================================================
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c
c------------------------  PROB CONTROL  -------------------------------
c
c   nlim   = cycles to run
c   tlim   = physical (problem) time to stop calculation
c cpulim   = CPU time in seconds to stop the calculation (default  3.6M)
c  tsave   = CPU time to reserve for terminating the run (default 30.0s)
c mbatch   = 0 interactive mode                          (default     0)
c          = 1 batch mode (does not scan for keyboard input)
c
      nlim   = 1 000 000
      tlim   = 0.0
      cpulim = 3 600 000.0
      tsave  = 30.0
      mbatch = 0
c
      nred   = 0
      if (myid_w .eq. 0) then
       read (1,pcon)
       write(2,pcon)
      endif
c
c------------------------  HYDRO CONTROL  ------------------------------
c
c  qcon   = quadratic artificial viscosity (q) constant
c  qlin   = linear    artificial viscosity (q) constant
c  courno = courant number
c  dtrat  = ratio of initial dt to dtmin (used to compute dtmin below)
c  iord   = default order of advection scheme for all variables
c  iCMA   = use consistent advection scheme of Plewa & Mueller, A&A, 342
c           179 to correct abundance conservation problems not addressed
c           by the MCONSRV module in some 9-species I-front problems. 
c           iCMA=0 is no consistent advection; iCMA=1 means H and He are
c           separately enforced to be fh and 1-fh; iCMA=2 means the sum
c           of the H, He, and metal abundances are required to be 1.
c  istp   = default steepening flag for all variables
c  iord** = order of advection scheme to be used for variable **
c  iostp**= steepening flag for 3rd order advection.  When istp**=1,
c           use the discontinuity detection to steepen shocks during
c           interpolation for variable ** in X1INT,X1INTFC,X2INT,X2INTFC
c  **floor = smallest value desired for variable ** on grid
c            Can also be used to set a default value for initialization.
c            Note that no attempt is made to ensure that actual values
c            stay above the floor values.
c
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
        d floor = tiny
        e floor = tiny
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
c
      read (1,hycon)
      write(2,hycon)
c
c Set flags to default values unless they were set in the input deck.
c
        if(iordd  .eq. 0) iordd  = iord
        if(iorde  .eq. 0) iorde  = iord
        if(iords1 .eq. 0) iords1 = iord
        if(iords2 .eq. 0) iords2 = iord
        if(iords3 .eq. 0) iords3 = iord
        if(iordb1 .eq. 0) iordb1 = iord
        if(iordb2 .eq. 0) iordb2 = iord
        if(iordb3 .eq. 0) iordb3 = iord
        if(iorder .eq. 0) iorder = iord
c
        if(istpd  .eq. 2) istpd  = istp
        if(istpe  .eq. 2) istpe  = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istper .eq. 2) istper = istp
c
c copy input flags to a buffer for later use and broadcasting.
c
#ifdef MPI_USED
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
c
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
        ibuf_in(21) = idebug
        ibuf_in(22) = iCMA
#endif
      endif
c
c Broadcast pcon and hycon to the others (use arrays).
c
#ifdef MPI_USED
      call MPI_BCAST( buf_in,20,MPI_FLOAT  ,0
     &                                       ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibuf_in,22,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (myid_w .ne. 0) then
        qcon    =  buf_in( 1)
        qlin    =  buf_in( 2)
        courno  =  buf_in( 3)
        dtrat   =  buf_in( 4)
        d floor =  buf_in( 5)
        e floor =  buf_in( 6)
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
c
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
        idebug  = ibuf_in(21)
        iCMA    = ibuf_in(22)
      endif
#endif
c------------------------------ RT CONTROL -----------------------------
c
c  ichem     = activate 9-species chemistry solver (off=0, on=1)--read
c              in by configure.F in the rchmconf namelist, not by setup 
c              in rtcon
c  iRT       = specify which RT scheme to apply (none=0, raytracing=1,
c              photon packet=2)--this is read in by configure.F in the 
c              rchmconf namelist, not by setup in rtcon
c  nnu1,nnu2 = number of photon bins below and above 13.6 eV, respectively.
c              These are also read into the rchmconf namelist rather than
c              in rtcon because they are used to allocate arrays at runtime,
c              as is iRT
c  ifreq     = specify monochromatic (=1), multifrequency (=2), or multigroup
c              (=3) transport 
c  ispct     = type of source spectrum (monochromatic=1, black body=2, QSO=3
c              time-dependent=4 --see spectrum.F). ispct=4 should only be 
c              used if ifreq=2 and for nspec > 6 (full H_2 chemistry); this
c              option requires that photon rates for the freqency bins be
c              read in from external tables (can be spectra for individual 
c              stars or spectral energy distributions (SED's) for stellar
c              populations or galaxies (but the data must be photons/sec for
c              each energy bin)
c  iOTS      = use the on-the-spot approximation (in other words, do not 
c              explicitly transport diffuse (recombination) radiation) (no=0,
c              yes=1)
c  iLW       = include Lyman-Werner radiation (no=0, single band=1, multiline
c              =2)
c  ibkgnd    = include uniform (metagalactic) Lyman-Werner background (no=0, 
c              yes=1).  Only a frequency-averaged flux is currently imple-
c              mented.
c  iHM       = include H- photodetachment radiation (no=0, yes=1)
c  iPWA      = in XYZ or ZRP geometries, attenuate the intensity of the plane 
c              wave according to its distance from its source by 1/r**2 (no=0,
c              yes=1)
c  iextinct  = include dust extinction in photon-conservative transfer (need
c	       user-supplied tables of tau vs hnu)(no=0, yes=1)
c  nsrc      = the number of stars comprising the UV source
c  nphdot    = emission rate of ionizing source photons (s^-1)
c  ephot     = energy (in eV) of source photons, if monochromatic
c  t_on      = time in the problem at which the radiation source should be
c              turned on (in sec)
c  t_off     = time in the problem at which the radiation source should be
c              turned off (in sec)
c  r_sep     = if we are assuming a plane wave UV source, specify r_sep (in
c              pc) to compute photon number flux at far left yz-face (by 
c              nphdot / (4 pi r_sep^2) if XYZ or ZRP are defined together
c              with RT in zeusmp.def)
c  t_star    = effective temperature of black body source (in K)
c  hnumax    = cutoff energy (in eV) for the normalization integral we 
c              compute to obtain an n(nu)dnu that is consistent with nphdot 
c  alpha     = specify QSO power-law spectrum (according to F(nu) being 
c              proportional to nu^(-alpha)
c  nLWdot    = emission rate of photons in the Lyman-Werner band
c  L_star    = log total luminosity (erg/s) of the star (needed to normalize  
c              the H- photodetachment spectrum, obtainable from Schaerer, 2002)
c  J_21      = frequency-averaged flux assumed for the uniform LW background
c              (in units of 10^-21 erg cm^-2 Hz^-1 str^-1, eg Haiman et al 2001)
c              Recall that 4*pi*J_nu = F_nu.
c

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
#ifdef MPI_USED
c
c Broadcast rtcon to the others.
c
      call MPI_BCAST(ibuf_in, 9,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( buf_in,11,MPI_FLOAT  ,
     .                0,MPI_COMM_WORLD,ierr)

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
#endif /* MPI_USED */
c
c------------------------ CHEMISTRY CONTROL  ---------------------------
c
c  nratec = number of bins in which temperatures are interpolated in
c           the cooling and rate tables computed by CTABLE.  Read in
c           by configure in the rchmconf namelist rather than here in
c           setup because it is used to allocate array sizes at run
c           time
c
c  iceco  = include collisional excitation cooling (yes=1, no=0)
c  icico  = include collisional ionization cooling (yes=1, no=0)
c  ireco  = include recombination cooling          
c  ibrco  = include bremmstrahlung cooling
c  ipiht  = include photoionization heating
c  ih2co  = include ro-vibrational cooling by H2 (0=none, 1=Lepp & Shull, 
c           2=Galli & Palla)
c  icmpt  = include Compton cooling          
c  iDMco  = include Dalgarno McCray (1972) cooling 
c  iHDco  = include HD cooling
c  idust  = include H_2 production by dust grains in reaction network (0=no,
c	    1=yes) (Jura, M. ApJ 1975, 197, 575.)  Gas cooling due to dust
c  icycle = choice of timestep hierarchy (1 for original, which is slightly
c           more accurate but much slower at transporting R-type I-fronts or
c           2 for optimized/compressed subcycling, with speedups of 2 - 15
c           for R-type fronts (and is a little faster with D-type fronts too))
c  rdshft = redshift (needed to specify Compton cooling)          
c  t_cutoff = temperature below which cooling isn't included in edot
c  fh     = mass fraction of hydrogen in the gas (fh = 1.0 for pure hydrogen
c           and 0.76 if primordial abundances for H and He are assumed)
c  z_sol  = fraction of solar metallicity assumed for gas (set to tiny for
c           zero-metallicity gas to avoid floating point exceptions)
c
c  Note:  electron collisional excitation of H and He is contained in the
c         Dalgarno & McCray cooling curves, so iceco should be set to 0 if
c         iDMco = 1.  Apparently, electron ionizational cooling of H and He
c         are not contained in DM so icico should remain on.

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
#ifdef MPI_USED
c
c Broadcast chemcon to the others.
c
      call MPI_BCAST(ibuf_in,11,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( buf_in,4,MPI_FLOAT  
     .,0,MPI_COMM_WORLD,ierr)

      if (myid_w .ne. 0) then
c
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
#endif /* MPI_USED */
c
c  Call ctable to set up chemistry coefficient tables
c
      call ctable

c
c-------------------------  SET MPI DATATYPES  -------------------------c
c Define MPI Derived Datatypes for passing 2-D slices of 3-D arrays.
c
#ifdef MPI_USED
      call MPI_TYPE_VECTOR (  jn*kn,1 ,in
     &, MPI_FLOAT, i_slice,ierr)
      call MPI_TYPE_COMMIT (i_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm ,in*neqm
     &, MPI_FLOAT, ils_slice,ierr)
      call MPI_TYPE_COMMIT (ils_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn*nspec,1 ,in
     &, MPI_FLOAT, iab_slice,ierr)
      call MPI_TYPE_COMMIT (iab_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm*neqm ,in*neqm*neqm
     &, MPI_FLOAT, ilsm_slice,ierr)
      call MPI_TYPE_COMMIT (ilsm_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in,in*jn,MPI_FLOAT
     &                                 ,j_slice,ierr)
      call MPI_TYPE_COMMIT (j_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in*neqm,in*jn*neqm,
     .                      MPI_FLOAT
     &                                 ,jls_slice,ierr)
      call MPI_TYPE_COMMIT (jls_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in*neqm*neqm,in*jn*neqm*neqm,
     .                      MPI_FLOAT
     &                                 ,jlsm_slice,ierr)
      call MPI_TYPE_COMMIT (jlsm_slice,ierr)
c
      call MPI_TYPE_VECTOR (    kn*nspec  ,in,in*jn,
     .                      MPI_FLOAT
     &                                 ,jab_slice,ierr)
      call MPI_TYPE_COMMIT (jab_slice,ierr)
c
      call MPI_TYPE_VECTOR (1,in*jn ,1    ,MPI_FLOAT
     &                                 ,k_slice,ierr)
      call MPI_TYPE_COMMIT (k_slice,ierr)
      call MPI_TYPE_VECTOR (1,in*jn*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,kls_slice,ierr)
      call MPI_TYPE_COMMIT (kls_slice,ierr)
      call MPI_TYPE_VECTOR (1,in*jn*neqm*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,klsm_slice,ierr)
      call MPI_TYPE_COMMIT (klsm_slice,ierr)
c
      call MPI_TYPE_VECTOR (nspec,in*jn ,in*jn*kn,
     .                      MPI_FLOAT
     &                                 ,kab_slice,ierr)
      call MPI_TYPE_COMMIT (kab_slice,ierr)
#endif /* MPI */
c
c-------------------------  READ GRID NAMELISTS  -----------------------
c
CJH   --- GRID data handled in RESTART ---
c
c------------------------  GRAVITY CONTROL  ----------------------------
c
c  Gravitational self-potentials can be included in both 1-D and 2-D
c  problems by solving the Poisson equation in the GRAVITY module.
c  Point mass potentials are included directly in the momentum eqn
c  by using a non-zero value for the variable ptmass.  Point mass
c  potentials do not require defining GRAV, do not call the GRAVITY
c  module, and are not included in the array phi but are explicitely
c  added to the momentum eqn terms in the routines STV1 and STV2.
c     tgrav  = time when gravitation is switched on
c     ptmass = fixed central point mass object
c     izero  = i index of Z=0 (x1=0) for odd symmetry case
c              (cylindrical geometry only)
c     igrijb = ijb flag (0 for     symmetric (Dirichlet) boundary      )
c                       (1 for non-symmetric             boundary whose
c                        value is calculated using multipole expansion )
c     igrojb = ojb flag ("  "      "          "         "              )
c     epsgrv = error limit              for ICCGAF
c     maxgrv = maximum iteration count  for ICCGAF
c     ks0grv = level of cyclic reduction in ICCGAF
c     phibverr = error criteria for multipole moments in PHIBV
c     phibvnmx = max number of moments to be taken in PHIBV
c     graverr = error criteria for computing a soln in GRAVITY
c
      if(myid_w .eq. 0) then
        tgrav     = 0.0        ! efh 99/05/14
        ptmass    = 0.0
        x1ptm     = x1a(is)
        x2ptm     = x2a(js)
        x3ptm     = x3a(ks)
c
       read (1,grvcon)
       write(2,grvcon)
      endif
c
c------------------------  RADIATION CONTROL  --------------------------
c
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
c
c------------------------  EQUATION OF STATE  --------------------------
c
c      gamma = ratio of specific heats
c      ciso  = isothermal sound speed
c
      if(myid_w .eq. 0) then
        gamma = 0.0
        ciso  = 0.0
c
        read (1,eqos)
        write(2,eqos)
      endif
#ifdef MPI_USED
c
c Broadcast grvcon, radcon, and eqos to the others.
c
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
      call MPI_BCAST( buf_in,23,MPI_FLOAT  ,0,comm3d,ierr)
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
#endif /* MPI */
      gamm1 = gamma - 1.0
c
c--------------RADIATION TILE BOUNDARY ARRAYS -------------------------
c
c  For the loops below, set the boundary flags to zero for internal
c  boundaries.  Copy the constant values into the 2-D arrays.
c
c----------------------  Constant values for IIB  ---------------------
c
      if(lrad .ne. 0) then
       do 40 k=ks-2,ke+3
         do 30 j=js-2,je+3
           if (coords(1) .eq. 0) then
             liib    (j,k) = niis(2)
           else
             liib    (j,k) = 0
           endif
30       continue
40     continue
c
c-----------------------  Constant values for OIB  ---------------------
c
       do 60 k=ks-2,ke+3
         do 50 j=js-2,je+3
           if (coords(1) .eq. ntiles(1) - 1) then
             loib    (j,k) = nois(2)
           else
             loib    (j,k) = 0
           endif
50       continue
60     continue
c
c-----------------------  Constant values for IJB  ---------------------
c
       do 80 k=ks-2,ke+3
         do 70 i=is-2,ie+3
           if (coords(2) .eq. 0) then
             lijb    (i,k) = nijs(2)
           else
             lijb    (i,k) = 0
           endif
70       continue
80     continue
c
c-----------------------  Constant values for OJB  ---------------------
c
       do 100 k=ks-2,ke+3
         do 90 i=is-2,ie+3
           if (coords(2) .eq. ntiles(2) - 1) then
             lojb    (i,k) = nojs(2)
           else
             lojb    (i,k) = 0
           endif
90       continue
100    continue
c
c-----------------------  Constant values for IKB  ---------------------
c
       do 120 j=js-2,je+3
         do 110 i=is-2,ie+3
           if (coords(3) .eq. 0) then
             likb    (i,j) = niks(2)
           else
             likb    (i,j) = 0
           endif
110      continue
120    continue
c
c-----------------------  Constant values for OKB  ---------------------
c
       do 140 j=js-2,je+3
         do 130 i=is-2,ie+3
           if (coords(3) .eq. ntiles(3) - 1) then
             lokb    (i,j) = noks(2)
           else
             lokb    (i,j) = 0
           endif
130      continue
140    continue
      endif ! lrad
c
c-------------------------  PROBLEM RESTART  ---------------------------
c
c  PROBRES is a user defined cpp macroname representing a subroutine
c  which re-intializes all user-defined variables for the particular 
c  problem to be restarted.  PROBRES should re-initialize only the
c  variables defined by the user in PROBLEM which are not in any
c  ZEUS-2D common blocks and are therefore not saved in the restart
c  dump.  PROBRES must read the same namelist(s) as PROBLEM.
c
      call PROBRES
c
c-------------------------  GRID MOTION CONTROL  -----------------------
c
c        igcon = 1 and x1fac < 0 gives "lagrangian" tracking in x1 lines
c        igcon = 1 and x2fac < 0 gives "lagrangian" tracking in x2 lines
c        igcon = 2 for input grid boundary speeds
c                vg1(io) = x1fac * central soundspeed
c                vg2(jo) = x2fac * central soundspeed
c        igcon = 3 for constant motion at x1[2]fac
c
      if(myid_w .eq. 0) then
        read (1,gcon)
        write(2,gcon)
      endif
c
c----------------------  COMPUTE TIMESTEP  -----------------------------
      if(ldimen .eq. 3) then
       kone = 1
      else
       kone = 0
      endif 
      if(ldimen .eq. 1) then
       jone = 0 
      else
       jone = 1
      endif
c
      do 370 k=ks,ke
       do 360 j=js,je
        do 350 i=is,ie
         w3da(i,j,k) = v1(i,j,k)*0.5*(d(i-1,j     ,k     )+d(i,j,k))
         w3db(i,j,k) = v2(i,j,k)*0.5*(d(i  ,j-jone,k     )+d(i,j,k))
     1                  * g2b(i)
         w3dc(i,j,k) = v3(i,j,k)*0.5*(d(i  ,j     ,k-kone)+d(i,j,k))
     1                  * g31b(i) * g32b(j)
         if(lrad .ne. 0) then
          w3dh(i,j,k) = er(i,j,k) / d(i,j,k)
         endif ! lrad
350     continue
360    continue
370   continue
c
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
      endif
c
      if(.not. xvgrid) then
       do i = 1, in
        vg1(i) = 0.0D0
       enddo
       do j = 1, jn
        vg2(j) = 0.0D0
       enddo
       do k = 1, kn
        vg3(k) = 0.0D0
       enddo
      endif
c
c      call nudt
c
      return
      end
