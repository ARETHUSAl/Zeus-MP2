!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          M S E N D R E C _ B N D          \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!     PURPOSE: update matrix elements in ghost zones for CG solver
!
!     Written by: F. Douglas Swesty and John Hayes
!
!=======================================================================
      subroutine msendrec_bnd( &
                        isx,iex,isy,iey,isz,iez, &
                        x)
!
      use real_prec
      use config
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif /* MPI_USED */
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez, i, j, k, n, nelm
!
      real(rl) :: x(neqm,neqm,in,jn,kn)
!
!-----------------------------------------------------------------------
!
#ifdef MPI_USED
      nelm = 1
!
!                               ! send south; receive north
      if(ntiles(1) .gt. 1) then
       if((periodic(1) .eqv. .true.) .or. (coords(1) .ne. 0)) then
        nreq = nreq + 1
        call mpi_isend( &
              x(1,1,isx,1,1),nelm,ilsm_slice,n1m,16300+nsub, &
              comm3d,req(nreq),ierr)
       endif
       if((periodic(1) .eqv. .true.) .or. (coords(1) .ne. ntiles(1)-1)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
              x(1,1,iex+1,1,1),nelm,ilsm_slice,n1p,16300+nsub, &
              comm3d,req(nreq),ierr)
       endif
!
!                               ! send north; receive south
       if((periodic(1) .eqv. .true.) .or. (coords(1) .ne. ntiles(1)-1)) &
        then
        nreq = nreq + 1
        call mpi_isend( &
              x(1,1,iex,1,1),nelm,ilsm_slice,n1p,16400+nsub, &
              comm3d,req(nreq),ierr)
       endif
       if((periodic(1) .eqv. .true.) .or. (coords(1) .ne. 0)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
              x(1,1,isx-1,1,1),nelm,ilsm_slice,n1m,16400+nsub, &
              comm3d,req(nreq),ierr)
       endif
      else ! ntiles(1)
       if((periodic(1) .eqv. .true.)) then
        do k = 1, kn
         do j = 1, jn
          x(1,1,isx-1,j,k) = x(1,1,iex  ,j,k)
          x(1,1,isx-2,j,k) = x(1,1,iex-1,j,k)
          x(1,1,iex+1,j,k) = x(1,1,isx  ,j,k)
          x(1,1,iex+2,j,k) = x(1,1,isx+1,j,k)
         enddo
        enddo
       endif
      endif ! ntiles(1)
!
      if(ldimen .gt. 1) then
!                            ! send east; receive west
      if(ntiles(2) .gt. 1) then
       if((periodic(2) .eqv. .true.) .or. (coords(2) .ne. 0)) &
        then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,isy,1),nelm,jlsm_slice,n2m,16500+nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if((periodic(2) .eqv. .true.) .or. (coords(2) .ne. ntiles(2)-1)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,iey+1,1),nelm,jlsm_slice,n2p,16500+nsub, &
                          comm3d,req(nreq),ierr)
       endif
!
!                            ! send west; receive east
       if((periodic(2) .eqv. .true.) .or. (coords(2) .ne. ntiles(2)-1)) &
        then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,iey,1),nelm,jlsm_slice,n2p,16600+nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if((periodic(2) .eqv. .true.) .or. (coords(2) .ne. 0)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,isy-1,1),nelm,jlsm_slice,n2m,16600+nsub, &
                          comm3d,req(nreq),ierr)
       endif
      else ! ntiles(2)
       if((periodic(2) .eqv. .true.)) then
        do k = 1, kn
         do i = 1, in
          x(1,1,i,isy-1,k) = x(1,1,i,iey  ,k)
          x(1,1,i,isy-2,k) = x(1,1,i,iey-1,k)
          x(1,1,i,iey+1,k) = x(1,1,i,isy  ,k)
          x(1,1,i,iey+2,k) = x(1,1,i,isy+1,k)
         enddo
        enddo
       endif
      endif ! ntiles(2)
!
      if(ldimen .gt. 2) then
!                            ! send above; receive below
      if(ntiles(3) .gt. 1) then
       if((periodic(3) .eqv. .true.) .or. (coords(3) .ne. ntiles(3)-1)) &
        then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,1,iez),nelm,klsm_slice,n3p,16700+nsub, &
               comm3d,req(nreq),ierr)
       endif
       if((periodic(3) .eqv. .true.) .or. (coords(3) .ne. 0)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,1,isz-1),nelm,klsm_slice,n3m,16700+nsub, &
               comm3d,req(nreq),ierr)
       endif
!
!                            ! send below; receive above
       if((periodic(3) .eqv. .true.) .or. (coords(3) .ne. 0)) &
        then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,1,isz),nelm,klsm_slice,n3m,16800+nsub, &
               comm3d,req(nreq),ierr)
       endif
       if((periodic(3) .eqv. .true.) .or. (coords(3) .ne. ntiles(3)-1)) &
        then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,1,iez+1),nelm,klsm_slice,n3p,16800+nsub, &
               comm3d,req(nreq),ierr)
       endif
      else ! ntiles(3)
       if((periodic(3) .eqv. .true.)) then
        do j = 1, jn
         do i = 1, in
          x(1,1,i,j,isz-1) = x(1,1,i,j,iez  )
          x(1,1,i,j,isz-2) = x(1,1,i,j,iez-1)
          x(1,1,i,j,iez+1) = x(1,1,i,j,isz  )
          x(1,1,i,j,iez+2) = x(1,1,i,j,isz+1)
         enddo
        enddo
       endif
      endif ! ntiles(3)
!
      endif ! ldimen > 2
      endif ! ldimen > 1
!
#endif /* MPI_USED */
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////          M S E N D R E C _ B N D          \\\\\\\\\\
!
!=======================================================================
