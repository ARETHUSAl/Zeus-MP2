!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////           L S _ P R C N _ B N D           \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!     PURPOSE:  compute diagonal preconditioning
!
!     Written by: F. Douglas Swesty and John Hayes
!
!=======================================================================
      subroutine sym_prcn_bnd(isx,iex,isy,iey,isz,iez, &
                             ipcflag,itrn, &
                             d, dp1, &
                                dp2, &
                                dp3, &
                             x,      rhs)
!
      use real_prec
      use config
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
!
      integer :: isx, iex, isy, iey, isz, iez, &
                 ipcflag, itrn
!
      real(rl) ::   d(neqm,neqm,in,jn,kn)
      real(rl) :: dp1(neqm,neqm,in,jn,kn)
      real(rl) :: dp2(neqm,neqm,in,jn,kn)
      real(rl) :: dp3(neqm,neqm,in,jn,kn)
      real(rl) :: x  (neqm,     in,jn,kn)
      real(rl) :: rhs(neqm,     in,jn,kn)
!
!                            loop indices
!
      integer :: i, j, k
!
!                            diagonal scale factor
!
      real(rl) :: atil
!
      if(ipcflag.eq.2 .or. ipcflag .eq. 0) then
        call sym_diag_bnd(isx,iex,isy,iey,isz,iez,ipcflag,itrn, &
                          d, dp1,dp2,dp3,x,rhs)
      else
        write(*,*)  ' ls_prcn: precon. flag ',ipcflag, 'not supported'
        write(*,*)  ' ls_prcn: i am stopping '
        stop
      endif
!
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////           L S _ P R C N _ B N D           \\\\\\\\\\
!
!=======================================================================
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////           L S _ P R C N _ I N T           \\\\\\\\\\
!
!=======================================================================
      subroutine sym_prcn_int(isx,iex,isy,iey,isz,iez, &
                             ipcflag,itrn, &
                             d, dp1, &
                                dp2, &
                                dp3, &
                             x,      rhs)
!
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit none
!
      integer  :: isx, iex, isy, iey, isz, iez, &
                  ipcflag, itrn
!
      real(rl) ::   d(neqm,neqm,in,jn,kn)
      real(rl) :: dp1(neqm,neqm,in,jn,kn)
      real(rl) :: dp2(neqm,neqm,in,jn,kn)
      real(rl) :: dp3(neqm,neqm,in,jn,kn)
      real(rl) :: x  (neqm,     in,jn,kn)
      real(rl) :: rhs(neqm,     in,jn,kn)
!
!                            loop indices
!
      integer  :: i, j, k
!
!                            diagonal scale factor
!
      real(rl) :: atil
!
      if(ipcflag .eq. 0 .or. ipcflag.eq.2) then
        call sym_diag_int(isx,iex,isy,iey,isz,iez,ipcflag,itrn, &
                         d, dp1,dp2,dp3,x,rhs)
      else
        write(*,*)  ' ls_prcn: precon. flag ',ipcflag, 'not supported'
        write(*,*)  ' ls_prcn: i am stopping '
        stop
      endif
!
 999  return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////           L S _ P R C N _ I N T           \\\\\\\\\\
!
!=======================================================================
