!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////           S E N D R E C _ B N D           \\\\\\\\\\
!
!=======================================================================
      subroutine updt_vec_bnd_1(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(neqm,in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                               ! send south; receive north
!-----------------------------------------------------------------------
!
      if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
              x(1,isx,1,1),nelm,ils_slice,n1m,15700+nsub, &
              comm3d,req(nreq),ierr)
      endif
      if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
              x(1,iex+1,1,1),nelm,ils_slice,n1p,15700+nsub, &
              comm3d,req(nreq),ierr)
      endif
!
!-----------------------------------------------------------------------
!                               ! send north; receive south
!-----------------------------------------------------------------------
!
      if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
        nreq = nreq + 1
         call mpi_isend( &
              x(1,iex,1,1),nelm,ils_slice,n1p,15800+nsub, &
              comm3d,req(nreq),ierr)
      endif
      if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
              x(1,isx-1,1,1),nelm,ils_slice,n1m,15800+nsub, &
              comm3d,req(nreq),ierr)
      endif
!
#endif /* MPI_USED */
!
 999  return
      end
!
      subroutine updt_vec_bnd_2(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(neqm,in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send east; receive west
!-----------------------------------------------------------------------
!
      if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,isy,1),nelm,jls_slice,n2m,15900+nsub, &
                          comm3d,req(nreq),ierr)
      endif
      if(nojs(3) .eq. 0 .or. nojs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,iey+1,1),nelm,jls_slice,n2p,15900+nsub, &
                          comm3d,req(nreq),ierr)
      endif
!
!-----------------------------------------------------------------------
!                            ! send west; receive east
!-----------------------------------------------------------------------
!
      if(nojs(3) .eq. 0 .or. nojs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,iey,1),nelm,jls_slice,n2p,16000+nsub, &
                          comm3d,req(nreq),ierr)
      endif
      if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,isy-1,1),nelm,jls_slice,n2m,16000+nsub, &
                          comm3d,req(nreq),ierr)
      endif
#endif /* MPI_USED */
!
 999  return
      end
!
      subroutine updt_vec_bnd_3(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(neqm,in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send above; receive below
!-----------------------------------------------------------------------
!
      if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,iez),nelm,kls_slice,n3p,16100+nsub, &
               comm3d,req(nreq),ierr)
      endif
      if(niks(3) .eq. 0 .or. niks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,isz-1),nelm,kls_slice,n3m,16100+nsub, &
               comm3d,req(nreq),ierr)
      endif
!
!-----------------------------------------------------------------------
!                            ! send below; receive above
!-----------------------------------------------------------------------
!
      if(niks(3) .eq. 0 .or. niks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,isz),nelm,kls_slice,n3m,16200+nsub, &
               comm3d,req(nreq),ierr)
      endif
      if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
         call mpi_irecv( &
               x(1,1,1,iez+1),nelm,kls_slice,n3p,16200+nsub, &
               comm3d,req(nreq),ierr)
      endif
#endif /* MPI_USED */
!
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////           S E N D R E C _ B N D           \\\\\\\\\\
!
!=======================================================================
      subroutine updt_gp_bnd_1(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                               ! send south; receive north
!-----------------------------------------------------------------------
!
       if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
              x(isx,1,1),nelm,i_slice,n1m,14000+25*nsub, &
              comm3d,req(nreq),ierr)
       endif
       if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
              x(iex+1,1,1),nelm,i_slice,n1p,14000+25*nsub, &
              comm3d,req(nreq),ierr)
       endif
!
!-----------------------------------------------------------------------
!                               ! send north; receive south
!-----------------------------------------------------------------------
!
       if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
        nreq = nreq + 1
         call mpi_isend( &
              x(iex,1,1),nelm,i_slice,n1p,14100+25*nsub, &
              comm3d,req(nreq),ierr)
       endif
       if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
              x(isx-1,1,1),nelm,i_slice,n1m,14100+25*nsub, &
              comm3d,req(nreq),ierr)
       endif
!
#endif /* MPI_USED */
!
 999  return
      end
!
      subroutine updt_gp_bnd_2(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send east; receive west
!-----------------------------------------------------------------------
!
       if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,isy,1),nelm,j_slice,n2m,14200+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if(nojs(3) .eq. 0 .or. nojs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,iey+1,1),nelm,j_slice,n2p,14200+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
!
!-----------------------------------------------------------------------
!                            ! send west; receive east
!-----------------------------------------------------------------------
!
       if(nojs(3) .eq. 0 .or. nojs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,iey,1),nelm,j_slice,n2p,14300+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,isy-1,1),nelm,j_slice,n2m,14300+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
#endif /* MPI_USED */
!
 999  return
      end
!
      subroutine updt_gp_bnd_3(isx,iex,isy,iey,isz,iez,x)
!
      use real_prec
      use config
      use param
      use bndry
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
      real(rl) :: x(in,jn,kn)
!
#ifdef MPI_USED 
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send above; receive below
!-----------------------------------------------------------------------
!
       if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,iez),nelm,k_slice,n3p,14400+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
       if(niks(3) .eq. 0 .or. niks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,isz-1),nelm,k_slice,n3m,14400+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
!
!-----------------------------------------------------------------------
!                            ! send below; receive above
!-----------------------------------------------------------------------
!
       if(niks(3) .eq. 0 .or. niks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,isz),nelm,k_slice,n3m,14500+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
       if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
         call mpi_irecv( &
               x(1,1,iez+1),nelm,k_slice,n3p,14500+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
#endif /* MPI_USED */
!
 999  return
      end
