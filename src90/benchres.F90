      subroutine benchres
      use real_prec
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use root
      implicit NONE
      real(rl) :: ampl0
      namelist /pgen/ ampl0

      ampl0=1.0e-4
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
        buf_in(1) = ampl0
      endif
!      call MPI_BCAST( ibuf_in, 8, MPI_INTEGER &
!                     , 0, comm3d, ierr )
      call MPI_BCAST(buf_in, 1, MPI_FLOAT &
                     , 0, comm3d, ierr )
      if (myid .ne. 0) then
        ampl0 = buf_in(1)
      endif
      return
      end
