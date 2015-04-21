c=======================================================================
c
c    \\\\\\\\\\        B E G I N   P R O G R A M          //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
      program zeusmp
c
c PURPOSE
c   Main program for 3-D MPI version of ZEUS.
c
c AUTHOR
c   Robert A. Fiedler
c
c LAST MODIFIED by JCH, for F90
c   6/26/02.
c
c updated 12.05.06 by Daniel Whalen to accommodate multispecies
c reactive flow and radiative transfer
c.......................................................................
c
c DECLARATIONS
c
      use real_prec
      use config
      use param
      use field
      use grid
      use root
      use cons
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use gravmod
      use clockmod
      use impsoln
      use chem
c
      implicit NONE
c
      real(rl4) :: cputime, wclock
c
      real(rl) :: zcs, etot, etot_glb
      real(rl) :: cpuall, rshock, vshkmax, rfront
c
      integer :: i , j , k, nwrite, ie_orig, ie_new
      integer :: nx, ny, nz, snz, maxn, ie_prev, ie_old
c
      ifsen(1) = 0
      ifsen(2) = 0
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
c
      myid_w    = 0
      myid      = 0
      nprocs_w  = 1
      nprocs    = 1
      coords(1) = 0
      coords(2) = 0
      coords(3) = 0
      reorder   = .false.
      totcgit = 0
      ncgcall = 0
c
c  call CONFIGURE
c
      call configure
c
c  Master writes greeting.
c
      if (myid_w .eq. 0) then
       call options
      endif
c
c Set up the problem: read input deck and possibly a restart file.
c
      call mstart
c
c Write out initial data dumps.
c
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
c Create best calculating plans for FFTW to perform forward
c and backward FFT.
c 
#ifdef FFT
      if(xgrvfft) then
        nx=ie-is+1
        ny=je-js+1
        nz=ke-ks+1
        call create_plan(ntiles(1)*nx,ntiles(2)*ny,ntiles(3)*nz)
      endif
#endif
c
c  Initialize cpu and wall clocks.  The quantities "cputime" and
c  "wclock" are the CPU and wall clock times (in seconds) in the main 
c  loop.
c
        wclock0 = 0.0
        cputime0 = 0.0
        call clocks (cputime, wclock)
        wclock0 = wclock
        cputime0 = cputime
c
      if (myid .eq. 0)
     &  write(6,"(/,' Set-up complete with ',i2,' warning(s):'
     &             ,' entering main loop...')") nwarn
c
      if(ldimen .eq. 1) nwrite = 1
      if(ldimen .eq. 2) nwrite = 1
      if(ldimen .eq. 3) nwrite = 1

 888  format(1p3e12.5)

c      open(unit=70,file='ifront.pos',status='unknown')


c
c--------------------------  start of main loop  -----------------------
c
c Execution ends when INTCHK returns a value of 1 for ifsen(1).
c
      ie_orig = ie
      ie_prev = is

1000  continue
      do i=is,ie
         if (abun(i,js,ks,1) .le. 0.5) rfront = x1b(i)
         if (v1(i,js,ks) .gt. v1(i-1,js,ks) .and.
     .       v1(i,js,ks) .gt. vshkmax)       then
           rshock  = x1b(i)
           vshkmax = v1(i,js,ks)
         endif

      enddo

c      write(70,888) time, rshock/cmkpc, rfront/cmkpc

      nsub = 1

c      if (time .lt. 4.7e12) then
c        ie = is+20
c      else 
c        ie = ie_orig
c      endif

c      ie_old = is
c      ie_new = is

c      do k=ks,ke
c      do j=js,je
c      do i=is,ie_orig

c         if (abun(i,j,k,2) .lt. 0.01 .and. dabs(v1(i,j,k)/cmkm) .lt. 
c     .       1.0d-5) then 
c           ie_new = i + 5
c           goto 70
c         endif

c      enddo

c      if (irestart .eq. 1 .and. ie_new .eq. is) ie_new = ie_orig

 70   continue
      
c      if (ie_old .ge. ie_new) ie_new = ie_old
c      ie_old = ie_new

c      enddo
c      enddo

cc      stop

c      if (ie_prev .ge. ie_new) then
c         ie_new  = ie_prev
c      else
c         ie_prev = ie_new
c      endif
c      ie = ie_new
c      if (ie .gt. ie_orig) ie = ie_orig
#ifdef MPI_USED
c
c Find the largest ie among all the tiles, and send the result
c to all in buf_out.  This preserves solution concurrency among
c all the tiles
c
c      ibuf_in(1) = ie
c      call MPI_ALLREDUCE( ibuf_in(1), ibuf_out(1), 1
c     .                   , MPI_2INTEGER
c     .                   , MPI_MAXLOC, comm3d, ierr)
c      ie = ibuf_out(1)
#endif /* MPI_USED */
cc      print*,"ie is: ",ie


c
c Call module that mimics the GSF 1996 paper approximate RT.
c

c      call gsfr
c
c      the lower line is the original one
c      if (ispct .ge. 4 .or. nhy .eq. 0) call spectrum
       if (ispct .ge. 4) call spectrum

c
c Call primordial chemistry solver
c Disabled to get the code compiling
c      if(ichem .eq. 1) then
c        call coolchem3D
c      endif 
c
c Solve Poisson's equation for the gravitational potential.
c
      if(xgrav .or. xgrvfft) then
       call gravity
      endif ! xgrav

c
c Call time-dependent boundary updates, if any
c
c      call SPECIALSRC
c      call SPECIALSRC1
c      call SPECIALSRC2

c
c Evaluate all non-advective terms in the evolution equations.
c
      if(lrad .eq. 0) then
       if(myid_w .eq. 0) then
        if(mod(nhy,nwrite) .eq. 0) then
         write(*,"('nhy, time, dt = ',i6,1p2d12.4)")
     .              nhy, time, dt
        endif ! mod
       endif ! myid_w
      endif ! lrad
c
      call srcstep

c      call vzero
c
c Compute the advection of all field variables.
c
      if(ldimen .eq. 3) call transprt
      if(ldimen .eq. 2) call transprt_2D
      if(ldimen .eq. 1) call transprt_1D

c      if(myid_w .eq. 0) write(*,*) 'transport done'
c      call vzero
c
c Update the step counter and current time value.
c
      nhy   = nhy   + 1
      time  = time  + dt
c
c Check the CPU time, number of steps, output times to determine if
c a stopping criterion has been met or output is desired.
c Also check for keyboard input, depending on the value of mbatch.
c
      call intchk( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6) )
c
c Compute new timestep
c
      call nudt

c      call vzero
c
c Update the grid and related quantites.
c
      if(xvgrid) call newgrid
c
c Write out any desired output files now that everything has been
c updated.
c Skip dataio if the run is being terminated to avoid duplicate output.
c
      ie = ie_orig
      if (ifsen(1) .eq. 1) goto 2000
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
      goto 1000  !  Loop back to begin the next time step.
c
c--------------------------  end of main loop  -------------------------
c
c Terminate the run by making final dumps, write goodbyes
c
2000  continue

c      close(unit=70)

      call clocks (cputime, wclock)
#ifndef ARCH_CRAY
      tused = real(cputime)
#else
      tused = wclock
#endif
      ifsen(2) = 1
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
      ifsen(6) = 1
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5) , ifsen(6))
c
#ifdef MPI_USED
c
c Sum everyone's cputime (stored in tused) to get CPU time for all 
c processes.
c
      call MPI_REDUCE(tused, cpuall, 1, MPI_FLOAT,
     &                MPI_SUM, 0, comm3d, ierr )
#else /* MPI */
       cpuall = tused
#endif /* MPI */
      if (myid .eq. 0) then      
c
c Let's assume tused is user + system time on the master thread.
c One would also like to keep track of wall-clock time and the sum
c of CPU times used by each processor.
c
        zcs = real(nprocs_w*nhy*nx1z*nx2z*nx3z)/(tused+tiny)
        write(6,"(/' Execution terminated with ',i4,' warning(s)')") 
     &     nwarn
        write(6,"(/' Performance summary:')")
        write(6,"('  zone-cycles per cpu second =',1pe12.5)") zcs
        write(6,"('  Master CPU                 =',1pe12.5, ' sec')") 
     &     tused
        write(6,"('  Average CPU/node           =',1pe12.5, ' sec')") 
     &     cpuall/real(nprocs_w)
        write(6,"('  Wall Clock                 =',1pe12.5, ' sec')") 
     &     wclock
        write(6,"()")
c
        if(xgrav .and. (.not. xsphgrv)) then
         if(lgeom .ne. 1) then
          write(6,"(/' GRAVITY SUMMARY:')")
          write(6,"('  Percentage of cycles with Phi updates: ',i3,
     .              '%')")int(100*float(ncgcall)/float(nhy))
          write(6,"('  Average number of iterations/update  : ',
     .              1pd12.4)")float(totcgit)/float(ncgcall)
         endif ! lgeom
        endif ! xgrav
        if(lrad .ne. 0) then
         write(6,"(/' RADIATION SUMMARY:')")
         write(6,"('  Average number of N-R iterations/cycle: ',
     .              1pd12.4)")totnrit/float(nhy)
         write(6,"('  Average number of CG  iterations/cycle: ',
     .              1pd12.4)")totlsit/float(nhy)
        endif ! lrad
c
        close(unit=2)
C        close(unit=3)
        close(unit=30)
       if(xtsl) then
        close(unit=31)
       endif ! xtsl
      endif
c
c Turn off MPI
c
#ifdef MPI_USED
      call MPI_FINALIZE ( ierr )
#endif
c
c=======================================================================
c
c    \\\\\\\\\\          E N D  P R O G R A M             //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c=======================================================================
c
      end
