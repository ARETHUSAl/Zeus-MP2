
!=======================================================================
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                  O P A C                  \\\\\\\\\!
!=======================================================================
!
      subroutine opacity(e, d, gam, kr, kp, sg, dkpdt, km, dkedt, &
                         ibeg, iend, jbeg, jend, kbeg, kend)
!
      use real_prec
      use config
      use param
      use cons
      use grid
      use opac_law
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, ibeg, iend, j, jbeg, jend, k, kbeg, kend
!
      real(rl) :: e    (in,jn,kn), d   (in,jn,kn), t    (in,jn,kn), &
                  p    (in,jn,kn), kp   (in,jn,kn), &
                  dtde (in,jn,kn), dpde(in,jn,kn), dbbdt(in,jn,kn), &
                  dkpdt(in,jn,kn), sg  (in,jn,kn), er   (in,jn,kn), &
                  kr   (in,jn,kn), km  (in,jn,kn), dkedt(in,jn,kn)
!
      real(rl) :: gam, coef, temp, binv, gam1, stef
      real(rl) :: rmfp0i
!
      real(rl) :: kap_floor, kap_ceil, sp_xscn, hundred
!
      real(rl) :: so_eps , tau    , eriibn, timarg
!
      common /soui/ so_eps, tau, eriibn
!
      hundred = 100.0d0
      stef = clight * rad_con / (4.D0 * pi)
      gam1 = gam - 1.0D0
!
      binv = 1.D0 / boltz
      do 20 k = kbeg, kend
      do 20 j = jbeg, jend
      do 20 i = ibeg, iend
       t(i,j,k) = (gam1)*mmw*mh*binv*e(i,j,k)/d(i,j,k)
20    continue
!
      rmfp0i  = 1.D0 / rmfp0
      do 50 k = kbeg, kend
      do 50 j = jbeg, jend
      do 50 i = ibeg, iend
         kr  (i,j,k) = rmfp0i * (d(i,j,k)/rho0)**xnu * &
                               (t(i,j,k)/ t_0)**(-powr)
         kr  (i,j,k) = max(min_coef, kr(i,j,k))
         kr  (i,j,k) = min(max_coef, kr(i,j,k))
         sg   (i,j,k) = 0.D0
         kp   (i,j,k) = min(hundred,kpfrac*kr(i,j,k))
         dkpdt(i,j,k) = -powr * kp(i,j,k) / t(i,j,k)
50    continue
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         km   (i,j,k) = kp   (i,j,k)
         dkedt(i,j,k) = dkpdt(i,j,k)
        enddo
       enddo
      enddo
!
      return
      end
!=======================================================================
!
!    \\\\\\\\\\        E N D  S U B R O U T I N E      //////////
!    //////////                  O P A C               \\\\\\\\\!
!=======================================================================
