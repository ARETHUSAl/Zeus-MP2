!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine mstart
!
!  PURPOSE:  Starts a run.
!
!  EXTERNALS: SETUP, MGET, RESTART
!
!  LOCALS:
!
!  LAST MODIFIED: 7/20/01 by PSLi
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use mpiyes
      use mpipar
!
      implicit NONE
!
!      integer  :: irestart
!
      integer incr,strtoi, iost
!
      namelist /mpitop/ ntiles, periodic
      namelist /rescon/ irestart,tdump,dtdump,id,resfile
      namelist /iocon/ thdf,dthdf,thist,dthist,tusr,dtusr,t_out &
                      ,ttsl,dttsl
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================
!
!
! Open input and log files -- master thread only.
!
      if (myid_w .eq. 0) then
        open(1,file='zmp_inp',status='old', iostat=iost)
!        write(*,*) 'IOSTAT = ',iost
        !open(2,file='zmp_log',status='unknown')
      endif
!JH
!------------------------  MPI TOPOLOGY  -------------------------------c
!  ntiles:   elements equal the number of tiles in each direction.
!  periodic: elements are true if grid is periodic in that direction --
!            we check hydro BC flags for defaults, but can override.
!
      ntiles(1) = 1
      ntiles(2) = 1
      ntiles(3) = 1
      periodic(1) = .false.
      periodic(2) = .false.
      periodic(3) = .false.
      if (myid_w .eq. 0) then 
        read(1,mpitop)
        write(2,mpitop)
      endif
!
!
! Tell the others what the master has read.
!
      call MPI_BCAST(ntiles  , 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(periodic, 3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!
! Quit if the number of processors indicated on the command line
! differs from the number of tiles specified in the input file.
!
      if (nprocs_w .ne. ntiles(1)*ntiles(2)*ntiles(3)) then
        if (myid_w .eq. 0) &
        write(*,"(/'MSTART: The number of threads ',i3,' does not match', &
                  /'MSTART: the number input via ntiles ',i3,' in', &
                  /'MSTART: input file zmp_inp; aborting the run...')") &
          nprocs_w, ntiles(1)*ntiles(2)*ntiles(3)
        call MPI_FINALIZE( ierr )
        stop
      endif
!
! Create a virtual Cartesian topology for the domain decomposition.
!
      call MPI_CART_CREATE( MPI_COMM_WORLD, 3, ntiles, periodic &
                          , reorder, comm3d, ierr )
      call MPI_COMM_RANK( comm3d, myid,     ierr )
      call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
!
! Find the ranks of my neighbors; find my virtual Cartesian coords.
!
      call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
      call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
      call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
!
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )
!
!
!------------------------  RESTART CONTROL  ----------------------------
!
!  irestart: set to one for calculation restarted from restart dump
!  tdump: time of last restart dump
!  dtdump: time between restart dumps
!  id: character*2 tag attended to filenames to identify run
!  resfile: name of restart file to restart from
!
      if (myid_w .eq. 0) then
        irestart = 0
        tdump = 0.0
        dtdump = 0.0
        id     = 'aa'
        resfile= 'resaa000000.0000'
        read (1,rescon)
        write(2,rescon)
      endif
!
! Tell the others what the master has read (do tdump, dtdump later). 
!
      call MPI_BCAST(irestart, 1,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(id      , 2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(resfile ,16,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!
! Read the remaining namelists (except for iocon) in setup/restart.
!
      if (myid_w .eq. 0 ) close(unit=1)
      if (irestart .eq. 0) then
        call setup
        time = 0.0D0
      else
        incr = strtoi(resfile,13,16)
        write(resfile,"(a3,a2,3i2.2,'.',i4.4)") 'res',id,coords(1) &
                                            ,coords(2),coords(3),incr
        call mget(resfile)
        nwarn = 0
        ifsen(1) = 0
        call restart
      endif
!
!------------------------  I/O CONTROL ---------------------------------
!
!  thdf: time of last HDF dump
!  dthdf: time between HDF dumps
!  ttsl: time of last tslice dump
!  dttsl: time between tslice dumps
!  thist: time of last history dump
!  dthist: time between history dumps
!  tusr: time of last user dump
!  dtusr: time between user dumps
!
      if (myid_w .eq. 0) then
        open(1, file='zmp_inp', status='old', iostat=iost)
        if (irestart .eq. 0) then
           thdf  = 0.0
          dthdf  = 0.0
           ttsl  = 0.0
          dttsl  = 0.0
           thist = 0.0
          dthist = 0.0
           tusr  = 0.0
          dtusr  = 3.1536e07
          do incr=1,nbuff-8
            t_out(incr) = 0.0
          enddo ! incr
        endif
        read (1,iocon)
        write(2,iocon)
        !dtusr  = 3.1536e07
      endif
      if (myid_w .eq. 0) then
        buf_in(1) = thdf
        buf_in(2) = dthdf
        buf_in(3) = thist
        buf_in(4) = dthist
        buf_in(5) = tusr
        buf_in(6) = dtusr
        buf_in(7) = tdump
        buf_in(8) = dtdump
        buf_in(9) = 0.0
        buf_in(10) = 0.0
        buf_in(11) = ttsl
        buf_in(12) = dttsl
        do incr=13,nbuff
          buf_in(incr) = t_out(incr-12)
        enddo ! incr
      endif
!
! Tell the others what the master has read.  An array
! is used to pass just one message.
!
      call MPI_BCAST(buf_in,nbuff,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD &
                    ,ierr)
      if (myid_w .ne. 0) then
          thdf  = buf_in(1)
         dthdf  = buf_in(2)
          thist = buf_in(3)
         dthist = buf_in(4)
          tusr  = buf_in(5)
         dtusr  = buf_in(6)
         tdump  = buf_in(7)
        dtdump  = buf_in(8)
        ttsl   = buf_in(11)
        dttsl  = buf_in(12)
        do incr=1,nbuff-12
          t_out(incr) = buf_in(incr+12)
        enddo ! incr
      endif
!
! Output file names are of the form "hdfidccc.n" for ease of
! use with graphics packages that can process a series of files.
!
      if (irestart .eq. 0) then
        incr = 0
        ifsen(2) = 1
        ifsen(3) = 1
        ifsen(4) = 1
        ifsen(5) = 1
        ifsen(6) = 1
        write(tslfile,"(a3,i4.4,a2)") 'tsl',incr,id
        write(resfile,"(a3,a2,3i2.2,'.',i4.4)") 'res',id,coords(1) &
                                            ,coords(2),coords(3),incr
        write(hdffile,"(a3,a2,3i2.2,'.',i4.4)") 'hdf',id,coords(1) &
                                            ,coords(2),coords(3),incr
        write(hstfile,"(a3,a2,3i2.2,'.',i4.4)") 'hst',id,coords(1) &
                                            ,coords(2),coords(3),incr
        write(usrfile,"(a3,a2,3i2.2,'.',i4.4)") 'usr',id,coords(1) &
                                            ,coords(2),coords(3),incr
      else
        ifsen(2) = 0
        ifsen(3) = 0
        ifsen(4) = 0
        ifsen(5) = 0
        ifsen(6) = 0
        incr = strtoi(tslfile,4,7) + 1
        write(tslfile,"(a3,i4.4,a2)") 'tsl',incr,id 
        incr = strtoi(resfile,13,16) + 1
        write(resfile,"(a3,a2,3i2.2,'.',i4.4)") 'res',id,coords(1) &
                                            ,coords(2),coords(3),incr
        incr = strtoi(hdffile,13,16) + 1
        write(hdffile,"(a3,a2,3i2.2,'.',i4.4)") 'hdf',id,coords(1) &
                                            ,coords(2),coords(3),incr
        incr = strtoi(usrfile,13,16) + 1
        write(usrfile,"(a3,a2,3i2.2,'.',i4.4)") 'usr',id,coords(1) &
                                            ,coords(2),coords(3),incr
      endif
      if(xtsl) then
       if (myid_w .eq. 0) then
         open (unit=31,file=tslfile,status='unknown')
       endif
      endif ! xtsl
!
      if(xhst) then
       open (unit=3, file=hstfile, status='unknown')
      endif ! xhst
      if (incr+1 .le. nbuff-10) then
!
! Adjust dump times to hit the next t_out, if defined non-zero.
!
        if (t_out(incr+1) .gt. 0.0) then
          tusr = t_out(incr+1) - dtusr
          thdf = t_out(incr+1) - dthdf
        endif
      endif
!
!  Close unit=1 (input deck).  Unit=2 (zeuslp) is kept open throughout
!  entire run to accept warning messages.  It is closed in zeusmp.
!
      if (myid_w .eq. 0) close(1)
!
! Write out initial data dumps.
!
      if(irestart .eq. 0) then
!      call dataio( iswres, iswhdf, iswhst, iswusr
!    .            )
       call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
       nhy = 0
      else
       nhy = 0
      endif
!
      return
      end
