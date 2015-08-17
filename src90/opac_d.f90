!======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               O P A C _ D                 \\\\\\\\\!======================================================================
!
!
      subroutine opac_d
!
!     Driver routine for computing opacities when radiation is included
!
!     Written by: John Hayes
!     modified: 4-13-99; added entry for "kem" to OPACITY calls
!     modified: 5-26-2003; rewritten for F90 code
!
!
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use bndry
      use scratch
      use mpiyes
      use mpipar
      use radiation
      use opac
      use cons
!
      implicit NONE
!
      integer  :: i,j,k
!
      real(rl) :: temp
!
!----------------------------------------------------------------------
!
!     update ghost zones along 1-coordinate
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,3,0,0,0,0,d)
       call bvale  (3,3,0,0,0,0,e)
!
!     compute opacities in all interior zones
!
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                    is, ie, js, je, ks, ke)
!
!     wait for communications to complete
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
      if(ldimen .eq. 1) then
!
!     compute opacities in 1-coordinate ghost zones and return
!
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                    is-2, is-1, js, je, ks, ke)
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                    ie+1, ie+2, js, je, ks, ke)
!
       go to 999
      endif ! ldimen = 1
      if(ldimen .gt. 1) then
!
!     update ghost zones along 2-coordinate
!
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,3,3,0,0,d)
       call bvale  (0,0,3,3,0,0,e)
!
!     compute opacities in 1-coordinate ghost zones
!
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                    is-2, is-1, js, je, ks, ke)
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                    ie+1, ie+2, js, je, ks, ke)
!
!     wait for communications to complete
!
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
       if(ldimen .eq. 2) then
!
!     compute opacities in 2-coordinate ghost zones and return
!
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, js-2, js-1, ks, ke)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, je+1, je+2, ks, ke)
!
        go to 999
       else ! ldimen > 2
        nreq = 0
        nsub = nsub + 1
        call bvald  (0,0,0,0,3,3,d)
        call bvale  (0,0,0,0,3,3,e)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, js-2, js-1, ks, ke)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, je+1, je+2, ks, ke)
        if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, js, je, ks-2, ks-1)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt, &
                     is, ie, js, je, ke+1, ke+2)
!
        go to 999
       endif ! ldimen > 2
      endif ! ldimen > 1
!
999   return
      end
