c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                E J P R O F                \\\\\\\\\\
c
c           Subroutine to calculate the supernova ejecta profile
c             at a fixed time.
c
c
c           The profile of the supernova ejecta is as follows:
c
c                  rho0 * t^(-3)                   v < v0 
c           rho = 
c                  rho0 * t^(-3) * (v/v0)^(-k)     v > v0  
c
c                     INPUT
c           E      :Ejecta energy                               [erg]
c           Mej    :Ejecta mass                                 [g]
c           k      :Slope of outer density drop                 [*] 
c           rmax   :Radius of last grid cell                    [cm]
c           rmin   :Radius of first grid cell                   [cm]
c           vmax   :Maximum velocity of supernova ejecta        [cm/s]
c           ncells :Number of grid cells
c
c
c                     OUTPUT
c           rho0   :Scaling constant                            [g/cm^3]
c           v0     :Scaling constant                            [cm/s]
c
c
c
c
c           Developed by Bob van Veelen
c             Last updated on: 20-03-2007
c    
c=======================================================================
c
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
c-----INPUT VARIABLES-------------
      real(rl) :: Eej, Mej, k, rmax, rmin, vmax
c-----OUTPUT VARIABLES------------
      real(rl) :: rho0, v0, tset
c-----PHYSICS CONSTANTS-----------
      real(rl) :: parsec
c-----LOOP VARIABLES
      real(rl) :: r(ncells), v(ncells), dr(ncells), rho(ncells)
      real(rl) :: mass, energy, vincr, vmin



      parsec     = cmpc
c------Grid parameters---------------------------------
      iter       = 0
      vincr      = vmax/100.

      tset = rmax / vmax
      vmin = rmin / tset
c     rmin = rmax * (vmin/vmax)

c      print*,"tset,vmin,rmin are: ",tset,vmin,rmin

      do i=1,ncells
      r(i)  = rmin + (real(i-1)/real(ncells-1)) * (rmax - rmin)
      v(i)  = r(i) / tset
       if(i.eq.1)then
        dr(i) = r(i+1) - r(i)
       else
        dr(i) = r(i) - r(i-1)
       endif
      enddo


c---First guess for v0------------------------------------
      v0         = v(1)
c---------------------------------------------------------

200   continue

      rho0       = 1.

c------Loop to calculate v and rho------------------------
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

c------Loop to calculate mass and energy using new rho0
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

c------------------------------------------------------


c------Checking if the solution is found---------------
c      if (iter.eq.0.or.iter.eq.1.e3) then
c      print*,'E/E_ej = ',energy/E
c      endif

      if (energy/Eej.gt.1.00005.or.energy/Eej.lt.0.99995) then
c      if (mass/Mej.gt.1.00005.or.mass/Mej.lt.0.99995) then
       if (iter.lt.1e7) then
        v0= v0 + vincr
        if (energy/Eej.gt.1.00005) then
c          print*,'Exceeded convergence limit',energy/Eej
c          print*,'Choosing new increment for v0'
          v0         = 0.95*v0
          vincr      = vincr/5.
c          iterations = 0
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
c       print*,'Solution found for supernova ejecta profile'
       goto 300
      endif
c------------------------------------------------------

300   continue
c      print*,'Vmax = ',v(ncells)/1.e5,' [km/s]'
c      print*,'Vo   = ',v0/1.e5,' [km/s]'
c      print*,'Vmin = ',v(1)/1.e5,' [km/s]'
c      print*,''
c      print*,'Number of time steps = ',ncells
c      print*,'Number of iterations = ',iter
c      print*,'Output written to supernova.dat'
c------------------------------------------------------

c      open(1,file='supernova.dat')

c      do i=1,ncells
c      write(1,2000)r(i),v(i),rho(i),99.
c      enddo
c      close(1)
      goto 500


400   continue
      print*,'No solution found!!!'
500   continue


2000  format(4e15.7)
      return
      end
