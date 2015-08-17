!=======================================================================
!
       subroutine clocks (cputime, wclock)
!
! This routine obtains the CPU and wall-clock times in seconds for the 
! calling process since the times (cputime0, wclock0) at which the 
! clocks were initalized (passed through common block /clocks/). 
!
! With MPI, MPI_WALL time is used except on the EXEMPLAR.
!
! For systems without a CPU timer, cputime is set to zero.
!
! Written by RAF, last modified 3/25/96.
!......................................................................
!
      use real_prec
      use param
      use config
      use clockmod
!
      implicit NONE
!
      real(rl4) :: cputime, wclock
!
!
      real(rl) :: etime_
!
!
      real(rl) :: MPI_WTIME, wall
!
!
! Wall clock time is easy to get with MPI:
!
      wall   = MPI_WTIME()
      wclock = real(wall) - wclock0
!
! Get the CPU time for this process/thread.
!
! For systems with the standard UNIX etime/itime routines.  Note that
! the output from itime is an integer array with (hours,minutes,seconds)
! and therefor is accurate only to the nearest second.  The output from
! function etime itself is the sum of user plus system times, which can
! be significantly longer than the user time alone (stored in 
! tarray(1)). 
!
!
      cputime = etime_ ( tarray ) - cputime0
!
      end
