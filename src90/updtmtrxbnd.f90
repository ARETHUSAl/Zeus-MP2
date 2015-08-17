!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          M S E N D R E C _ B N D          \\\\\\\\\!
!=======================================================================
      subroutine updt_mtrx_bnd_1(isx,iex,isy, iey, isz, iez, x)
!
      use real_prec
      use config
      use param
      use bndry
      use mpiyes
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
      nelm = 1
!
!                               ! send south; receive north
      if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
       nreq = nreq + 1
       call mpi_isend( &
             x(1,1,isx,1,1),nelm,ilsm_slice,n1m,16300+25*nsub, &
             comm3d,req(nreq),ierr)
      endif
      if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
       nreq = nreq + 1
       call mpi_irecv( &
             x(1,1,iex+1,1,1),nelm,ilsm_slice,n1p,16300+25*nsub, &
             comm3d,req(nreq),ierr)
      endif
!
!                               ! send north; receive south
!
      if(nois(3) .eq. 0 .or. nois(3) .eq. 4) then
       nreq = nreq + 1
       call mpi_isend( &
             x(1,1,iex,1,1),nelm,ilsm_slice,n1p,16400+25*nsub, &
             comm3d,req(nreq),ierr)
      endif
      if(niis(3) .eq. 0 .or. niis(3) .eq. 4) then
       nreq = nreq + 1
       call mpi_irecv( &
             x(1,1,isx-1,1,1),nelm,ilsm_slice,n1m,16400+25*nsub, &
             comm3d,req(nreq),ierr)
      endif
!
 999  return
      end
!
      subroutine updt_mtrx_bnd_2(isx, iex, isy, iey, isz, iez, x)
!
      use real_prec
      use config
      use param
      use bndry
      use mpiyes
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez, i, j, k, n, nelm
!
      real(rl) :: x(neqm,neqm,in,jn,kn)
!
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send east; receive west
!-----------------------------------------------------------------------
!
       if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,isy,1),nelm,jlsm_slice,n2m,16500+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if(nojs(3) .eq. 0 .or. nojs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,iey+1,1),nelm,jlsm_slice,n2p,16500+25*nsub, &
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
               x(1,1,1,iey,1),nelm,jlsm_slice,n2p,16600+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
       if(nijs(3) .eq. 0 .or. nijs(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,isy-1,1),nelm,jlsm_slice,n2m,16600+25*nsub, &
                          comm3d,req(nreq),ierr)
       endif
!
 999  return
      end
!
      subroutine updt_mtrx_bnd_3(isx, iex, isy, iey, isz, iez, x)
!
      use real_prec
      use config
      use param
      use bndry
      use mpiyes
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez, i, j, k, n, nelm
!
      real(rl) :: x(neqm,neqm,in,jn,kn)
!
      nelm = 1
!
!-----------------------------------------------------------------------
!                            ! send above; receive below
!-----------------------------------------------------------------------
!
       if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_isend( &
               x(1,1,1,1,iez),nelm,klsm_slice,n3p,16700+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
       if(niks(3) .eq. 0 .or. niks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,1,isz-1),nelm,klsm_slice,n3m,16700+25*nsub, &
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
               x(1,1,1,1,isz),nelm,klsm_slice,n3m,16800+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
       if(noks(3) .eq. 0 .or. noks(3) .eq. 4) then
        nreq = nreq + 1
        call mpi_irecv( &
               x(1,1,1,1,iez+1),nelm,klsm_slice,n3p,16800+25*nsub, &
               comm3d,req(nreq),ierr)
       endif
!
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////          M S E N D R E C _ B N D          \\\\\\\\\!
!=======================================================================
