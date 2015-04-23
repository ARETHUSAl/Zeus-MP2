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
#ifndef ARCH_CRAY
      real(rl4) :: cputime, wclock
#else
      real(rl ) :: cputime, wclock
#endif
!
#ifdef ARCH_IBM
#define ETIME   etime_
#else
#define ETIME   etime
#endif
!
#ifndef ARCH_TERAGRID
      real(rl) :: ETIME
#endif
!
!
#ifdef MPI_USED
      real(rl) :: MPI_WTIME, wall
#endif
!
#ifdef MPI_USED
!
! Wall clock time is easy to get with MPI:
!
      wall   = MPI_WTIME()
#ifndef ARCH_CRAY
      wclock = real(wall) - wclock0
#else
      wclock = wall - wclock0
#endif
#endif /* MPI_USED */
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
#ifndef MPI_USED
      cputime = 0.0
#endif
!
#ifndef ARCH_TERAGRID
#ifndef ARCH_CRAY
      cputime = ETIME ( tarray ) - cputime0
#else
      cputime = wclock
#endif
#endif
!
      end
