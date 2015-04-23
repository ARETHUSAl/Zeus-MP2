!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                E J P R O F                \\\\\\\\\\
!
!           Subroutine to calculate the supernova ejecta profile
!             at a fixed time.
!
!
!           The profile of the supernova ejecta is as follows:
!
!                  rho0 * t^(-3)                   v < v0 
!           rho = 
!                  rho0 * t^(-3) * (v/v0)^(-k)     v > v0  
!
!                     INPUT
!           E      :Ejecta energy                               [erg]
!           Mej    :Ejecta mass                                 [g]
!           k      :Slope of outer density drop                 [*] 
!           rmax   :Radius of last grid cell                    [cm]
!           rmin   :Radius of first grid cell                   [cm]
!           vmax   :Maximum velocity of supernova ejecta        [cm/s]
!           ncells :Number of grid cells
!
!
!                     OUTPUT
!           rho0   :Scaling constant                            [g/cm^3]
!           v0     :Scaling constant                            [cm/s]
!
!
!
!
!           Developed by Bob van Veelen
!             Last updated on: 20-03-2007
!    
!=======================================================================
!
      subroutine ejprof(Eej, Mej, k, rmax, rmin, vmax, ncells, rho0, v0)
      use real_prec
      use param
      use config
      use root
      use field
      use bndry
      use grid
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      implicit none
      integer i, iter, ncells
!-----INPUT VARIABLES-------------
      real(rl) :: Eej, Mej, k, rmax, rmin, vmax
!-----OUTPUT VARIABLES------------
      real(rl) :: rho0, v0, tset
!-----PHYSICS CONSTANTS-----------
      real(rl) :: parsec
!-----LOOP VARIABLES
      real(rl) :: r(ncells), v(ncells), dr(ncells), rho(ncells)
      real(rl) :: mass, energy, vincr, vmin
      parsec     = cmpc
!------Grid parameters---------------------------------
      iter       = 0
      vincr      = vmax/100.
      tset = rmax / vmax
      vmin = rmin / tset
!     rmin = rmax * (vmin/vmax)
!      print*,"tset,vmin,rmin are: ",tset,vmin,rmin
      do i=1,ncells
      r(i)  = rmin + (real(i-1)/real(ncells-1)) * (rmax - rmin)
      v(i)  = r(i) / tset
       if(i.eq.1)then
        dr(i) = r(i+1) - r(i)
       else
        dr(i) = r(i) - r(i-1)
       endif
      enddo
!---First guess for v0------------------------------------
      v0         = v(1)
!---------------------------------------------------------
200   continue
      rho0       = 1.
!------Loop to calculate v and rho------------------------
      do i=1,ncells
        if (v(i).lt.v0) then
         rho(i)  = rho0
        else
         rho(i)  = rho0*(v(i)/v0)**(-k)
        endif
      enddo
      mass   = 0.
      energy = 0.
      do i=1,ncells
       mass  = mass+4.*pi*r(i)*r(i)*rho(i)*dr(i)
       energy= energy+0.5*rho(i)*v(i)*v(i)*4.*pi*r(i)*r(i)*dr(i)
      enddo
      rho0   = Mej/mass
      mass   = 0.
      energy = 0.
!------Loop to calculate mass and energy using new rho0
      do i=1,ncells
        if (v(i).lt.v0) then
         rho(i)  = rho0
        else
         rho(i)  = rho0*(v(i)/v0)**(-k)
        endif
      enddo
      do i=1,ncells
       mass  = mass+4.*pi*r(i)*r(i)*rho(i)*dr(i)
       energy= energy+0.5*rho(i)*v(i)*v(i)*4.*pi*r(i)*r(i)*dr(i)
      enddo
!------------------------------------------------------
!------Checking if the solution is found---------------
!      if (iter.eq.0.or.iter.eq.1.e3) then
!      print*,'E/E_ej = ',energy/E
!      endif
      if (energy/Eej.gt.1.00005.or.energy/Eej.lt.0.99995) then
!      if (mass/Mej.gt.1.00005.or.mass/Mej.lt.0.99995) then
       if (iter.lt.1e7) then
        v0= v0 + vincr
        if (energy/Eej.gt.1.00005) then
!          print*,'Exceeded convergence limit',energy/Eej
!          print*,'Choosing new increment for v0'
          v0         = 0.95*v0
          vincr      = vincr/5.
!          iterations = 0
          goto 200
        endif
        if (v0.gt.v(ncells).or.v0.lt.v(1)) then
         print*,'Velocity falls out of range!'
         goto 400
        else
         iter=iter+1
        goto 200
        endif
       else
        print*,'Maximum number of iterations reached'
        goto 400
       endif
      else
!       print*,'Solution found for supernova ejecta profile'
       goto 300
      endif
!------------------------------------------------------
300   continue
!      print*,'Vmax = ',v(ncells)/1.e5,' [km/s]'
!      print*,'Vo   = ',v0/1.e5,' [km/s]'
!      print*,'Vmin = ',v(1)/1.e5,' [km/s]'
!      print*,''
!      print*,'Number of time steps = ',ncells
!      print*,'Number of iterations = ',iter
!      print*,'Output written to supernova.dat'
!------------------------------------------------------
!      open(1,file='supernova.dat')
!      do i=1,ncells
!      write(1,2000)r(i),v(i),rho(i),99.
!      enddo
!      close(1)
      goto 500
400   continue
      print*,'No solution found!!!'
500   continue
2000  format(4e15.7)
      return
      end
