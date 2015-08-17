!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////             C O N F I G U R E           \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!     PURPOSE:  Reads the code configuration namelists that debuted
!               with ZEUS-MP Version 2.  These namelists serve the
!               functions formerly performed by CPP macros in the
!               discarded "zeusmp.def" file used to configure
!               ZEUS-MP Version 1.
!
!     Written by:  John Hayes; way back in '03.
!     Modified by: John Hayes; repeatedly since then.
!
!=======================================================================
!
      subroutine configure
!
      use real_prec
      use param
      use config
      use mpiyes
      use mpipar
      use chem
!
      implicit NONE
!
      integer ::  confi_buf(15)
      logical ::  confl_buf(20)
      real(rl)::  confr_buf(2)
!
      namelist /geomconf/  lgeom, ldimen
      namelist /physconf/  lrad   , xhydro , xgrav, xmhd  , xgrvfft, &
                           xptmass, xtotnrg, xiso , xvgrid, xsubav , &
                           xforce , xsphgrv, xhse , leos  , nspec  , &
                           lopac , xchem
      namelist /rchmconf/  ichem  , iRT    , nnu1 , nnu2  , nratec
      namelist /ioconf/    xascii , xhdf, xrestart, xtsl, xhst
      namelist /preconf/   small_no, large_no
      namelist /arrayconf/ izones, jzones, kzones, maxijk
!
!----------------------------------------------------------------------
!     If parallel execution, start up MPI
!----------------------------------------------------------------------
!
       call MPI_INIT( ierr )
       call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
!
!----------------------------------------------------------------------
!     Open zmp_conf run configuration file
!----------------------------------------------------------------------
!
      if (myid_w .eq. 0) then 
        open(unit=1,file='zmp_inp',status='old')
        open(unit=2,file='zmp_log',status='unknown')
      endif
!
!
!----------------------------------------------------------------------
!     Initialize all parameters to default values before continuing
!     with configuration READ
!----------------------------------------------------------------------
!
      lgeom    = 1
      ldimen   = 3
      lrad     = 0
      leos     = 1
      lopac    = 0
      nspec    = 1
      izones   = 32
      jzones   = 32
      kzones   = 32
      maxijk   = 32
      xvgrid   = .false.
      xhydro   = .false.
      xforce   = .false.
      xgrav    = .false.
      xgrvfft  = .false.
      xsphgrv  = .false.
      xptmass  = .false.
      xmhd     = .false.
      xhse     = .false.
      xtotnrg  = .false.
      xiso     = .false.
      xascii   = .true.
      xhdf     = .false.
      xhst     = .false.
      xrestart = .false.
      xtsl     = .false.
      ichem    = 0
      iRT      = 0
      nratec   = 600
      nnu1     = 0
      nnu2     = 1
      large_no = 1.0D99
      small_no = 1.0D-99
!
!----------------------------------------------------------------------
!     Read remaining namelists
!----------------------------------------------------------------------
!
      if(myid_w .eq. 0) then
       read(1,geomconf)
       write(2,geomconf)
       confi_buf(1) = lgeom
       confi_buf(2) = ldimen
!
       read(1,physconf)
       write(2,physconf)
       confi_buf(3) = lrad
       confi_buf(4) = lopac
       confi_buf(5) = leos
       confi_buf(6) = nspec
       confl_buf(1) = xhydro
       confl_buf(2) = xforce
       confl_buf(3) = xmhd
       confl_buf(4) = xgrav
       confl_buf(5) = xgrvfft
       confl_buf(6) = xptmass
       confl_buf(7) = xtotnrg
       confl_buf(8) = xiso
       confl_buf(9) = xvgrid
       confl_buf(10) = xsubav
       confl_buf(11) = xsphgrv
       confl_buf(12) = xhse
!
       read(1,ioconf)
       confl_buf(13) = xascii
       confl_buf(14) = xhdf
       confl_buf(15) = xrestart
       confl_buf(16) = xtsl
       confl_buf(17) = xhst
!
       read(1,preconf)
       confr_buf(1) = small_no
       confr_buf(2) = large_no
!
       read(1,arrayconf)
       confi_buf(12) = izones
       confi_buf(13) = jzones
       confi_buf(14) = kzones
       confi_buf(15) = maxijk
!
       read(1,rchmconf)
       confi_buf(7) = nratec
       confi_buf(8) = nnu1
       confi_buf(9) = nnu2
       confi_buf(10)= iRT
       confi_buf(11)= ichem
!      close zmp_conf run configuration file
       close(1)
!
      endif ! myid_w
!
       call MPI_BCAST(confi_buf,15,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(confl_buf,17,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(confr_buf,2,MPI_DOUBLE_PRECISION,0, &
                                  MPI_COMM_WORLD,ierr)
       if(myid_w .ne. 0) then
        lgeom  = confi_buf(1)
        ldimen = confi_buf(2)
        lrad   = confi_buf(3)
        lopac  = confi_buf(4)
        leos   = confi_buf(5)
        nspec  = confi_buf(6)
        nratec = confi_buf(7)
        nnu1   = confi_buf(8)
        nnu2   = confi_buf(9)
        iRT    = confi_buf(10)
        ichem  = confi_buf(11)
        izones = confi_buf(12)
        jzones = confi_buf(13)
        kzones = confi_buf(14)
        maxijk = confi_buf(15)
        xhydro  = confl_buf(1)
        xforce  = confl_buf(2)
        xmhd    = confl_buf(3)
        xgrav   = confl_buf(4)
        xgrvfft = confl_buf(5)
        xptmass = confl_buf(6)
        xtotnrg = confl_buf(7)
        xiso    = confl_buf(8)
        xvgrid  = confl_buf(9)
        xsubav  = confl_buf(10)
        xsphgrv = confl_buf(11)
        xhse    = confl_buf(12)
        xascii   = confl_buf(13)
        xhdf     = confl_buf(14)
        xrestart = confl_buf(15)
        xtsl     = confl_buf(16)
        xhst     = confl_buf(17)
        small_no = confr_buf(1)
        large_no = confr_buf(2)
       endif ! myid_w
!
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
      call alloc_arrays
!
      return
      end
