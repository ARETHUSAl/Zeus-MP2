!-----------------------------------------------------------------------
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
       subroutine auto_h4
!
! ZEUS-MP Post-processor: AUTO_H4.
!
!      Purpose:  this is like h4splice, except that the DO loop
!      over dump number is automated such that all dumps beginning
!      with Nbeg and ending with Nend are processed, where Nbeg
!      and Nend are input by the user.
!
!.......................................................................
!
      use real_prec
      use config
      use param
      use mpino
      use mpipar
      use grid
      use root
!
      implicit NONE
!
      !integer   :: irestart
!
      integer   :: nbl,igrid,imin,imax,jmin,jmax,kmin,kmax
      logical   :: lgrid
!
      integer   :: iskip,jskip,kskip,izone,jzone,kzone
      integer   :: n, indx, icmb, ndump, nbeg, nend
      integer   :: it,jt,kt,i,j,k,itil
      integer   :: nfunc,lfunc,ret,shape(3),rank,rcmb,shcmb(3)
!
      real(rl4), dimension(:), allocatable :: xscale,yscale,zscale
!
      integer, parameter :: max_pe = 256
      integer, parameter :: max_fc = 257
!
      real(rl4), dimension(:), allocatable :: xscmb,yscmb, zscmb
      real(rl4), dimension(:), allocatable :: data ,dcmb
!
      real(rl) :: x1min,x1max,x1rat,dx1min,dfndx1r,x1r,deltx1r,errx1r, &
                  x2min,x2max,x2rat,dx2min,dfndx2r,x2r,deltx2r,errx2r, &
                  x3min,x3max,x3rat,dx3min,dfndx3r,x3r,deltx3r,errx3r
!
      character*8 phrase
      character*23 hdf(max_pe)
      character*15 cmb
      character*120 line
      character strng(max_fc)*32,qunit*1,qfrmt*1,qcoord*16
!
      integer  dsgdims,dsgdisc,dsgdata,dsgdast
      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
      external dsgdims,dsgdisc,dsgdata,dsgdast
      external dssdims,dssdast,dssdisc,dsadata,dspdata
!PS
      integer ntotal,inz,jnz,knz
!
      namelist /geomconf/  lgeom, ldimen
      namelist /physconf/  lrad   , xhydro , xgrav, xmhd    , xgrvfft, &
                           xptmass, xtotnrg, xiso , xvgrid  , xsubav, &
                           xforce , xsphgrv, leos , nspec, xchem
      namelist /ioconf/    xascii , xhdf, xrestart, xtsl
      namelist /preconf/   small_no, large_no
      namelist /arrayconf/ izones, jzones, kzones, maxijk
      namelist /rescon/ irestart,tdump,dtdump,id,resfile
      namelist /mpitop/ ntiles, periodic
      namelist /ggen1/  nbl,x1min,x1max,igrid,x1rat,dx1min,lgrid
      namelist /ggen2/  nbl,x2min,x2max,igrid,x2rat,dx2min,lgrid
      namelist /ggen3/  nbl,x3min,x3max,igrid,x3rat,dx3min,lgrid
!
!.......................................................................
!
! Read ZEUS-MP input file "zmp_inp" to determine the number of zones
! in each direction in the physical mesh (nx1z,nx2z,nx3z), and the 
! number of tiles in each direction "ntiles".
!
      open(1,file='zmp_inp',status='old')
!
!     read the code configuration namelists
!
      read(1,geomconf)
      read(1,physconf)
      read(1,ioconf)
      read(1,preconf)
      read(1,arrayconf)
!
      in = izones + 5
      jn = jzones + 5
      kn = kzones + 5
!
      allocate(xscmb(max_pe*in))
      allocate(yscmb(max_pe*jn))
      allocate(zscmb(max_pe*kn))
      allocate(data (in*jn*kn))
      allocate(dcmb(max_pe*in*jn*kn))
      allocate(xscale(in))
      allocate(yscale(jn))
      allocate(zscale(kn))
!
!
!   Look for and read the "mpitop" and "ggen*" namelists.
!
!       do 10 n=1,100
!         read(1,"(a8)") phrase
!         if (phrase(3:8) .eq. 'MPITOP' .or. &
!             phrase(3:8) .eq. 'mpitop') goto 20
!   10  continue
!       write(*,"(/'HSPLICE: ERROR -- Did not find mpitop namelist.')")
!       close (1)
!       return
!!
!   20  continue
!       backspace (1)
       read(1,mpitop)
       write(*,"('HSPLICE: ntiles is ',3i6)") ntiles(1),ntiles(2) &
       , ntiles(3)
!
!   Read the "rescon" namelist to find the run's "id".
!
       irestart = 0
        tdump = 0.0
       dtdump = 0.0
       id     = 'aa'
       resfile= 'resaa000000.000'
       read (1,rescon)
       write(*,"(/'HSPLICE: id is ',a)") id
!
       imax  = 3
       imin = imax
       nbl   = 1
       lgrid = .false.
!
   30  continue
       read(1,ggen1)
       imax = imax + nbl
       if (.not. lgrid) go to 30
!
       nx1z = ( imax - imin ) / ntiles(1)
!PS
!
       jmax  = 3
       jmin  = jmax
       nbl   = 1
       lgrid = .false.
!
   40  continue
       read(1,ggen2)
       jmax = jmax + nbl
       if (.not. lgrid) go to 40
!
       nx2z = ( jmax - jmin ) / ntiles(2)
!
!
       kmax  = 3
       kmin = kmax
       nbl   = 1
       lgrid = .false.
!
   50  continue
       read(1,ggen3)
       kmax = kmax + nbl
       if (.not. lgrid) go to 40
!
       nx3z = ( kmax - kmin ) / ntiles(3)
!PS
!
       write(*,"('HSPLICE: nx1z,nx2z,nx3z are ',3i6)") nx1z,nx2z,nx3z
       close(1)
!.......................................................................
!
! Ask for the dump number.  Loop back to this point when finished
! with this dump number.  Negative input dump number terminates
! execution.
!
   60  continue
!
      write(*,"(/'AUTO_HSPLICE: Enter beginning and ending dump', &
                  ' numbers'/)")
!
      read(*,*) nbeg, nend
!
! Generate names of hdf files to read.
!
      DO NDUMP = nbeg, nend
       indx = 0
       do 90 kt=0,ntiles(3)-1
         do 80 jt=0,ntiles(2)-1
           do 70 it=0,ntiles(1)-1
             indx = indx + 1
             write(hdf(indx),"(a2,i4.4,a1,a3,a2,3i2.2,'.',i3.3)") & 
                   'DD',ndump,'/','hdf',id, it,jt,kt,ndump
   70      continue
   80    continue
   90  continue
       nprocs = indx
!.......................................................................
!
! Read each function from each tile's dump and write it into
! the combined file 'cmb' ("hdf<id>.<ndump>").
!
       write(cmb,"(a3,a2,'.',i3.3)") 'hdf',id,ndump
!
!   Loop over the number of functions to be read -- determine number.
!
       nfunc = 5
       if(xmhd) nfunc = nfunc + 3
       if(lrad .ne. 0) nfunc = nfunc + 1
!
       do lfunc=1,nfunc
!
!   For each tile, read one function.  On the first pass, read
!   and save the mesh points for each tile.
!
!PS
         indx = 0
         inz=0
         jnz=0
         knz=0
         do kt=0,ntiles(3)-1
           do jt=0,ntiles(2)-1
             do it=0,ntiles(1)-1
               indx = indx + 1
               if (lfunc.eq.1) then
!
                 ret = dsgdims(hdf(indx),rank,shape,3)
                 ret = dsgdisc(1,shape(1),xscale)
                 ret = dsgdisc(2,shape(2),yscale)
                 ret = dsgdisc(3,shape(3),zscale)
!
!         Copy the mesh points for the tile into the cmb mesh arrays.
!
!PS
                 do i=1,nx1z
                   xscmb(i+it*nx1z) = xscale(i)
                 enddo ! i
                 do j=1,nx2z
                   yscmb(j+jt*nx2z) = yscale(j)
                 enddo ! j
                 do k=1,nx3z
                   zscmb(k+kt*nx3z) = zscale(k)
                 enddo ! k
               endif ! lfunc = 1
!
!
!         Read the tile HDF file for the lfunc-th function.
!
               do n=1,lfunc
                 ret = dsgdata(hdf(indx),rank,shape,data)
                 ret = dsgdast(strng(lfunc),qunit,qfrmt,qcoord)
               enddo ! n
!
!         Copy data values for this tile into the cmb data array.
!
               do k=1,nx3z
                 do j=1,nx2z
                   do i=1,nx1z
                     itil       = i + (j-1)*nx1z + (k-1)*nx1z*nx2z
!PS
                     icmb       = i + it*nx1z &
                                + (j+jt*nx2z-1)*(ntiles(1)*nx1z) &
                                + (k+kt*nx3z-1)*(ntiles(1)*nx1z) &
                                * (ntiles(2)*nx2z)
!
                     dcmb(icmb) = data(itil)
                   enddo ! i
                 enddo ! j
               enddo ! k
!
!PS
             enddo ! it
           enddo ! jt
         enddo ! kt
!
!
!   For the first function only, name the cmb mesh HDF file and
!   prepare it to receive data (after reading all tiles mesh points).
!
         rcmb     = 3
!PS
         shcmb(1) = nx1z*ntiles(1)
         shcmb(2) = nx2z*ntiles(2)
         shcmb(3) = nx3z*ntiles(3)
!
         ret = dssdims(rcmb,shcmb)
         ret = dssdisc(1,shcmb(1),xscmb)
         ret = dssdisc(2,shcmb(2),yscmb)
         ret = dssdisc(3,shcmb(3),zscmb)
         ret = dssdast(strng(lfunc),' ',' ',qcoord)
         if (lfunc.eq.1) then
!
!   Write the cmb mesh data into the combined file for this function.
!
           ret = dspdata(cmb,rcmb,shcmb,dcmb)
         else ! lfunc = 1
!
!   Append the cmb mesh data into the combined file for this function.
!
           ret = dsadata(cmb,rcmb,shcmb,dcmb)
         endif ! lfunc = 1
       enddo ! lfunc
!
      ENDDO ! NDUMP
!
      return
!.......................................................................
!
! Error messages
!
  999  continue
       write(*,"(/'HSPLICE: ERROR -- could not open file ',a15 &
       , ' for tile ',i4)") hdf(indx),indx
!
       return
       end
