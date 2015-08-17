!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                S T R E A M                \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
!
      subroutine stream
!
!  PURPOSE:  Problem generator for optically thin propagation of a
!  square pulse of radiation energy density.
!
!  Converted for ZEUS-MP from Jim Stone's ZEUS-2D version.
!
!  RAF, 2/18/97.
!
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use bndry
      use mpiyes
      use mpipar
      use cons
      use opac_law
!
      implicit NONE
!
      integer  :: i, j, k, idirect
      real(rl) :: qa,amp,x0
      real(rl) :: ros_mfp, flx_lim, dmc_max, dx_min
!
      namelist /pgen/ idirect,amp,x0
      idirect = 1
      amp     = 1.0
      x0      = 0.1
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
         buf_in( 1) = amp
         buf_in( 2) = x0
        ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 2, MPI_DOUBLE_PRECISION &
                     , 0, comm3d, ierr )
       call MPI_BCAST(ibuf_in, 1, MPI_INTEGER &
                     , 0, comm3d, ierr )
       if (myid .ne. 0) then
         amp     =  buf_in( 1)
         x0      =  buf_in( 2)
         idirect = ibuf_in( 1)
       endif
!
      qa = dfloor*8.625e7/(gamma-1.0)*(clight/(4.0*pi*1.8044e-5))**0.25
      do 1 k=1,kn
      do 1 j=1,jn
      do 1 i=1,in
        er(i,j,k) = erfloor
        if (idirect .eq. 1 .and. x1b(i) .le. x0) er(i,j,k) = er(i,j,k) &
          + amp
        if (idirect .eq. 2 .and. x2b(j) .le. x0) er(i,j,k) = er(i,j,k) &
          + amp
        if (idirect .eq. 3 .and. x3b(k) .le. x0) er(i,j,k) = er(i,j,k) &
          + amp
        e (i,j,k) = qa*er(i,j,k)**0.25
1     continue
!
      if (idirect .eq. 1) then
        do 10 k=1,kn
        do 10 j=1,jn
          niib(j,k) =  3
          liib(j,k) =  3
          noib(j,k) =  2
          loib(j,k) =  2
          eriib(j,k,1) = fiis(12)
          eriib(j,k,2) = fiis(12)
          eiib(j,k,1) = qa*eriib(j,k,1)**0.25
          eiib(j,k,2) = qa*eriib(j,k,2)**0.25
10      continue
      endif
!
      if (idirect .eq. 2) then
        do 20 k=1,kn
        do 20 i=1,in
          nijb(i,k) = 3
          lijb(i,k) = 3
          nojb(i,k) = 2
          lojb(i,k) = 2
          erijb(i,k,1) = fijs(12)
          erijb(i,k,2) = fijs(12)
          eijb(i,k,1) = qa*erijb(i,k,1)**0.25
          eijb(i,k,2) = qa*erijb(i,k,2)**0.25
20      continue
      endif
!
      if (idirect .eq. 3) then
        do 30 j=1,jn
        do 30 i=1,in
          nikb(i,j) = 3
          likb(i,j) = 3
          nokb(i,j) = 2
          lokb(i,j) = 2
          erikb(i,j,1) = fiks(12)
          erikb(i,j,2) = fiks(12)
          eikb(i,j,1) = qa*erikb(i,j,1)**0.25
          eikb(i,j,2) = qa*erikb(i,j,2)**0.25
30      continue
      endif
!
      return
      end
