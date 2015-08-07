!=======================================================================
!
!    \\\\\\\\\\        B E G I N   P R O G R A M          //////////
!    //////////               Z E U S M P                 \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
program zeusmp
!
! PURPOSE
!   Main program for 3-D MPI version of ZEUS.
!
! AUTHOR
!   Robert A. Fiedler
!
! LAST MODIFIED by JCH, for F90
!   6/26/02.
!
! updated 12.05.06 by Daniel Whalen to accommodate multispecies
! reactive flow and radiative transfer
!.......................................................................
!
! DECLARATIONS
!
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
!
      implicit NONE
!
      real(rl4) :: cputime, wclock
!
      real(rl) :: zcs, etot, etot_glb
      real(rl) :: cpuall, rshock, vshkmax=0, rfront
!
      integer :: i , j , k, nwrite, ie_orig, ie_new
      integer :: nx, ny, nz, snz, maxn, ie_prev, ie_old
!
      ifsen(1) = 0
      ifsen(2) = 0
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
!
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
!
!  call CONFIGURE
!
      call configure
!
!  Master writes greeting.
!
      if (myid_w .eq. 0) then
       call options
      endif
!
! Set up the problem: read input deck and possibly a restart file.
!
      call mstart
!
! Write out initial data dumps.
!
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
!
! Create best calculating plans for FFTW to perform forward
! and backward FFT.
! 
#ifdef FFT
      if(xgrvfft) then
        nx=ie-is+1
        ny=je-js+1
        nz=ke-ks+1
        call create_plan(ntiles(1)*nx,ntiles(2)*ny,ntiles(3)*nz)
      endif
#endif
!
!  Initialize cpu and wall clocks.  The quantities "cputime" and
!  "wclock" are the CPU and wall clock times (in seconds) in the main 
!  loop.
!
        wclock0 = 0.0
        cputime0 = 0.0
        call clocks (cputime, wclock)
        wclock0 = wclock
        cputime0 = cputime
!
      if (myid .eq. 0) &
        write(6,"(/,' Set-up complete with ',i2,' warning(s):' &
                   ,' entering main loop...')") nwarn
!
      if(ldimen .eq. 1) nwrite = 1
      if(ldimen .eq. 2) nwrite = 1
      if(ldimen .eq. 3) nwrite = 1
 888  format(1p,3e12.5)
!      open(unit=70,file='ifront.pos',status='unknown')
!
!--------------------------  start of main loop  -----------------------
!
! Execution ends when INTCHK returns a value of 1 for ifsen(1).
!
      ie_orig = ie
      ie_prev = is
1000  continue
      do i=is,ie
         if (abun(i,js,ks,1) .le. 0.5) rfront = x1b(i)
         if (v1(i,js,ks) .gt. v1(i-1,js,ks) .and. &
             v1(i,js,ks) .gt. vshkmax)       then
           rshock  = x1b(i)
           vshkmax = v1(i,js,ks)
         endif
      enddo
!      write(70,888) time, rshock/cmkpc, rfront/cmkpc
      nsub = 1
!      if (time .lt. 4.7e12) then
!        ie = is+20
!      else 
!        ie = ie_orig
!      endif
!      ie_old = is
!      ie_new = is
!      do k=ks,ke
!      do j=js,je
!      do i=is,ie_orig
!         if (abun(i,j,k,2) .lt. 0.01 .and. dabs(v1(i,j,k)/cmkm) .lt. 
!     .       1.0d-5) then 
!           ie_new = i + 5
!           goto 70
!         endif
!      enddo
!      if (irestart .eq. 1 .and. ie_new .eq. is) ie_new = ie_orig
 70   continue
!      if (ie_old .ge. ie_new) ie_new = ie_old
!      ie_old = ie_new
!      enddo
!      enddo
!c      stop
!      if (ie_prev .ge. ie_new) then
!         ie_new  = ie_prev
!      else
!         ie_prev = ie_new
!      endif
!      ie = ie_new
!      if (ie .gt. ie_orig) ie = ie_orig
#ifdef MPI_USED
!
! Find the largest ie among all the tiles, and send the result
! to all in buf_out.  This preserves solution concurrency among
! all the tiles
!
      ibuf_in(1) = ie
      call MPI_ALLREDUCE( ibuf_in(1), ibuf_out(1), 1 &
                         , MPI_2INTEGER &
                         , MPI_MAXLOC, comm3d, ierr)
      ie = ibuf_out(1)
#endif /* MPI_USED */
!c      print*,"ie is: ",ie
!
! Call module that mimics the GSF 1996 paper approximate RT.
!
!      call gsfr
!
!      the lower line is the original one
!      if (ispct .ge. 4 .or. nhy .eq. 0) call spectrum
       if (ispct .ge. 4) call spectrum
!
! Call primordial chemistry solver
! Disabled to get the code compiling
!      if(ichem .eq. 1) then
!        call coolchem3D
!      endif 
!
! Solve Poisson's equation for the gravitational potential.
!
      if(xgrav .or. xgrvfft) then
       call gravity
      endif ! xgrav
!
! Call time-dependent boundary updates, if any
!
!      call SPECIALSRC
!      call SPECIALSRC1
!      call SPECIALSRC2
!
! Evaluate all non-advective terms in the evolution equations.
!
      if(lrad .eq. 0) then
       if(myid_w .eq. 0) then
        if(mod(nhy,nwrite) .eq. 0) then
         write(*,"('nhy, time, dt = ',i6,1p,2d12.4)") &
                    nhy, time, dt
        endif ! mod
       endif ! myid_w
      endif ! lrad
!
      call srcstep
!      call vzero
!
! Compute the advection of all field variables.
!
      if(ldimen .eq. 3) call transprt
      if(ldimen .eq. 2) call transprt_2D
      if(ldimen .eq. 1) call transprt_1D
!      if(myid_w .eq. 0) write(*,*) 'transport done'
!      call vzero
!
! Update the step counter and current time value.
!
      nhy   = nhy   + 1
      time  = time  + dt
!
! Check the CPU time, number of steps, output times to determine if
! a stopping criterion has been met or output is desired.
! Also check for keyboard input, depending on the value of mbatch.
!
      call intchk( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6) )
!
! Compute new timestep
!
      call nudt
!      call vzero
!
! Update the grid and related quantites.
!
      if(xvgrid) call newgrid
!
! Write out any desired output files now that everything has been
! updated.
! Skip dataio if the run is being terminated to avoid duplicate output.
!
      ie = ie_orig
      if (ifsen(1) .eq. 1) goto 2000
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
!
      goto 1000  !  Loop back to begin the next time step.
!
!--------------------------  end of main loop  -------------------------
!
! Terminate the run by making final dumps, write goodbyes
!
2000  continue
!      close(unit=70)
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
!
#ifdef MPI_USED
!
! Sum everyone's cputime (stored in tused) to get CPU time for all 
! processes.
!
      call MPI_REDUCE(tused, cpuall, 1, MPI_FLOAT, &
                      MPI_SUM, 0, comm3d, ierr )
#else /* MPI */
       cpuall = tused
#endif /* MPI */
      if (myid .eq. 0) then      
!
! Let's assume tused is user + system time on the master thread.
! One would also like to keep track of wall-clock time and the sum
! of CPU times used by each processor.
!
        zcs = real(nprocs_w*nhy*nx1z*nx2z*nx3z)/(tused+tiny)
        write(6,"(/' Execution terminated with ',i4,' warning(s)')") &
           nwarn
        write(6,"(/' Performance summary:')")
        write(6,"('  zone-cycles per cpu second =',1pe12.5)") zcs
        write(6,"('  Master CPU                 =',1pe12.5, ' sec')") &
           tused
        write(6,"('  Average CPU/node           =',1pe12.5, ' sec')") &
           cpuall/real(nprocs_w)
        write(6,"('  Wall Clock                 =',1pe12.5, ' sec')") &
           wclock
        write(6,"()")
!
        if(xgrav .and. (.not. xsphgrv)) then
         if(lgeom .ne. 1) then
          write(6,"(/' GRAVITY SUMMARY:')")
          write(6,"('  Percentage of cycles with Phi updates: ',i3, &
                    '%')")int(100*float(ncgcall)/float(nhy))
          write(6,"('  Average number of iterations/update  : ', &
                    1pd12.4)")float(totcgit)/float(ncgcall)
         endif ! lgeom
        endif ! xgrav
        if(lrad .ne. 0) then
         write(6,"(/' RADIATION SUMMARY:')")
         write(6,"('  Average number of N-R iterations/cycle: ', &
                    1pd12.4)")totnrit/float(nhy)
         write(6,"('  Average number of CG  iterations/cycle: ', &
                    1pd12.4)")totlsit/float(nhy)
        endif ! lrad
!
        close(unit=2)
!        close(unit=3)
        close(unit=30)
       if(xtsl) then
        close(unit=31)
       endif ! xtsl
      endif
!
! Turn off MPI
!
#ifdef MPI_USED
      call MPI_FINALIZE ( ierr )
#endif
!
!=======================================================================
!
!    \\\\\\\\\\          E N D  P R O G R A M             //////////
!    //////////               Z E U S M P                 \\\\\\\\\\
!
!=======================================================================
!
end program zeusmp
