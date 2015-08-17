!=======================================================================
!
!                            Developed by
!                Laboratory for Computational Astrophysics
!                  University of California at San Diego
!
      subroutine gravity
!
!-----------------------------------------------------------------------
!     driver routine which selects an algorithm for computing the 
!     gravitational potential, gp(i,j,k), based upon user-specified
!     values of XGRAV, XSPHGRV, XGRVFFT, LDIMEN, and LGEOM:
!
!     General potentials:
!
!     LGEOM = 1 (XYZ)
!      (1) Dirichlet/Neuman BC's --> use grav3D_MG
!      (2) Triply-periodic  BC's --> use fftwgrav
!
!     LGEOM = 2 (ZRP):
!      (1) LDIMEN = 2 --> use grav2D_CG
!      (2) LDIMEN = 3 --> use grav3D_CG
!
!     LGEOM = 3 (RTP):
!      (1) LDIMEN = 2 --> use grav2D_CG
!      (2) LDIMEN = 3 --> use grav3D_CG
!
!     Spherically symmetric potential (GM/r):
!
!     LGEOM = 3 (RTP)
!      (1) LDIMEN = 1 --> use spherical_gravity
!      (2) LDIMEN = 2 --> use spherical_gravity
!-----------------------------------------------------------------------
!
      use real_prec
      use config
      use param
      use cons
      use domain
      use root
      use field
      use grid
      use bndry
      use impsoln
      use mpiyes
      use mpipar
!
      implicit NONE
!
!-----------------------------------------------------------------------
!     spherical potential (GM/r); supported if LGEOM = 3 and 
!     LDIMEN = 1 or 2
!-----------------------------------------------------------------------
!
      if(xsphgrv) then
!
!     check for 1D or 2D...
!
       if(ldimen .eq. 3) then
        if(myid .eq. 0) then
         write(*,"('XSPHGRV: LDIMEN must equal 1 or 2!')")
         write(*,"('aborting run...')")
        endif
        go to 1
       endif
!
!     check for RTP geometry...
!
       if(lgeom .ne. 3) then
        if(myid .eq. 0) then
         write(*,"('XSPHGRV: LGEOM must equal 3!')")
         write(*,"('aborting run...')")
        endif
        go to 1
       endif
!
       call spherical_gravity
       return
!
1      continue
       call mpi_finalize(ierr)
       stop
      endif
!
!-----------------------------------------------------------------------
!     potential for triply-periodic 3D problems (LGEOM =1)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     Solve Poisson's equation for 2D/3D problems which are not
!     triply-periodic
!-----------------------------------------------------------------------
!
      if(ldimen .eq. 2) then
!
!     2D, so must be spherical or cylindrical
!
       if(lgeom .eq. 1) then
        if(myid .eq. 0) then
         write(*,"('XGRAV: Can not do cartesian gravity in 2D!')")
         write(*,"('aborting run...')")
        endif
        go to 3
       endif
!
       call grav2D_CG
       return
!
3      continue
       call mpi_finalize(ierr)
       stop
      endif
!
!-----------------------------------------------------------------------
!     In theory, any geometry is legal here
!-----------------------------------------------------------------------
!
      if(ldimen .eq. 3) then
       if(lgeom .eq. 1) then
        call grav3D_CG
       else
        call grav3D_CG
       endif
      endif
!
      return
      end
!
!=======================================================================
      subroutine spherical_gravity
!
!     computes M(r) for use in forces.F in the case of 1D or 2D 
!     spherical problems, in which spherical gravity is assumed
!     in place of computing a potential
!
!     Written by: John Hayes
!
!     Modified 1: 12/12/05, fixed typo in "dml" expression
!     for 2D case when xwedge=false.  (deleted a "dcos2" reference).
!
      use real_prec
      use config
      use param
      use cons
      use field
      use grid
      use bndry
      use gravmod
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, index
!
      real(rl), dimension(:,:), allocatable :: dml, dmu
      real(rl), dimension(:  ), allocatable :: intm_half
!
!
      real(rl) :: local_mass, theta1, theta2, alpha1, alpha2, &
                  dcos1, dcos2
!
      real(rl) :: third, halfpi, onep5pi
!
      allocate(dml (in,jn))
      allocate(dmu (in,jn))
!
      do i = 1, in
       intm(i) = 0.0D0
      enddo
!
      third   = 1.0D0/3.0D0
      halfpi  = 0.5D0*pi
      onep5pi = 1.5D0*pi
!
      allocate(intm_half(2*in))
!
!-----------------------------------------------------------------------
!     differential mass elements (lower half cell and upper half cell)
!-----------------------------------------------------------------------
!
      if(ldimen .eq. 2) then
       if(xwedge) then
        do j = js, je
         theta1 = x2a(j)
         theta2 = x2a(j+1)
         if(x2a(j+1) .le. halfpi) then
          alpha1 = halfpi - x2a(j)
          alpha2 = halfpi - x2a(j+1)
          dcos1    = cos(theta1) - cos(theta2)
          dcos2    = cos(alpha2) - cos(alpha1)
         else
          alpha1 = onep5pi - x2a(j)
          alpha2 = onep5pi - x2a(j+1)
          dcos1    = cos(theta1) - cos(theta2)
          dcos2    = cos(alpha2) - cos(alpha1)
         endif
         do i = is, ie
          dmu(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1a(i+1)**3 - x1b(i)**3) &
                   * (dcos1 + dcos2)
         enddo
         if(coords(1) .eq. 0) then
          do i = is+1, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3) &
                    * (dcos1 + dcos2)
          enddo
          dml(is,j) = 2.0D0*third*pi*d(is,j,ks)*x1b(is)**3 &
                   * (dcos1 + dcos2)
         else ! coords(1)
          do i = is, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3) &
                    * (dcos1 + dcos2)
          enddo
         endif ! coords(1)
        enddo
       else ! xwedge
        do j = js, je
         theta1 = x2a(j)
         theta2 = x2a(j+1)
         dcos1    = cos(theta1) - cos(theta2)
         do i = is, ie
          dmu(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1a(i+1)**3 - x1b(i)**3) &
                   * dcos1
         enddo
         if(coords(1) .eq. 0) then
          do i = is+1, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3) &
                    * dcos1
          enddo
          dml(is,j) = 2.0D0*third*pi*d(is,j,ks)*x1b(is)**3 &
                    * dcos1
         else ! coords(1)
          do i = is, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3) &
                    * dcos1
          enddo
         endif ! coords(1)
        enddo
       endif ! xwedge
      else ! ldimen
       do i = is, ie
        dmu(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1a(i+1)**3 - x1b(i)**3)
       enddo
       if(coords(1) .eq. 0) then
        do i = is+1, ie
         dml(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1b(i)**3 - x1a(i)**3)
        enddo
        dml(is,js) = 4.0D0*third*pi*x1b(is)**3*d(is,js,ks)
       else
        do i = is, ie
         dml(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1b(i)**3 - x1a(i)**3)
        enddo
       endif ! coords(1)
      endif ! ldimen
!
!-----------------------------------------------------------------------
!     total interior mass (local processor) contained between x1a(is)
!     and x1a(ie+1), computed at every half cell interface
!-----------------------------------------------------------------------
!
      do i = 1, 2*in
       intm_half(i) = 0.0D0
      enddo
      local_mass = 0.0D0
!
      index = 0
      do i = is, ie
       index = index + 2
       do j = js, je
        local_mass = local_mass + dml(i,j)
       enddo
       intm_half(index  ) = local_mass
       do j = js, je
        local_mass = local_mass + dmu(i,j)
       enddo
       intm_half(index+1) = local_mass
      enddo
!
!
!-----------------------------------------------------------------------
!     sum mass arrays across rows of processors (2-coord)
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(comm3d, ierr)
!
      if(ntiles(2) .gt. 1) call sum_over_two_coord(intm_half,local_mass)
!
!-----------------------------------------------------------------------
!     sum mass arrays up 1-coordinate
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(comm3d, ierr)
!
      if(ntiles(1) .gt. 1) call sum_over_one_coord(intm_half,local_mass)
!
!
!-----------------------------------------------------------------------
!     pick out interior mass at cell FACES
!-----------------------------------------------------------------------
!
      index = -1
      do i = is, ie+1
       index = index + 2
       intm(i) = intm_half(index)
      enddo
!
      deallocate(intm_half)
      deallocate(dml)
      deallocate(dmu)
!
      return
      end
!
!=======================================================================
      subroutine sum_over_two_coord(intm_half,local_mass)
!
      use real_prec
      use param
      use grid
      use mpiyes
      use mpipar
!
      implicit NONE
!
      real(rl) :: intm_half(2*in)
      real(rl), dimension(:), allocatable :: intm_buf
      real(rl) :: local_mass, loc_buf
!
      integer source, target, tag, i, l, nhalfcells, count
!
      count = 2*in
      allocate(intm_buf(count))
!
      if(coords(2) .ne. 0) then
!
!-----------------------------------------------------------------------
!     process with (coords(1),coords(2)) = (N,M) sends data to 
!     process with (coords(1),coords(2)) = (N,0) 
!-----------------------------------------------------------------------
!
       target = coords(1)*ntiles(2)
       tag    = 1000+coords(1)
       call MPI_SEND(intm_half(1), count, &
                     MPI_DOUBLE_PRECISION, &
                     target, tag, comm3d, ierr)
      else  ! coords(2)
!
!-----------------------------------------------------------------------
!     process with coords(2) = 0 computes theta-summed values of local
!     interior mass
!-----------------------------------------------------------------------
!
       tag    = 1000+coords(1)
       do l = 1, ntiles(2)-1
       target = myid+l
        call MPI_RECV(intm_buf(1), count, &
                     MPI_DOUBLE_PRECISION, &
                     target, tag, comm3d, stat, ierr)
        do i = 1, 2*in
         intm_half(i) = intm_half(i) + intm_buf(i)
        enddo ! i
       enddo ! l
      endif ! coords(2)
!
      call MPI_BARRIER(comm3d, ierr)
!
      if(coords(2) .eq. 0) then
!
!-----------------------------------------------------------------------
!     process   with (coords(1),coords(2)) = (N,0) sends interior mass to
!     processes with (coords(1),coords(2)) = (N,M)
!-----------------------------------------------------------------------
!
       do l = 1, ntiles(2)-1
        target = coords(1)*ntiles(2) + l
        tag    = 2000 + l
        call MPI_SEND(intm_half(1), count, &
                     MPI_DOUBLE_PRECISION, &
                     target, tag, comm3d, ierr)
       enddo ! l
      else  ! coords(2)
!
!-----------------------------------------------------------------------
!     processes with coords(2) > 0 receive incoming messages and update
!     local values of intm_half
!-----------------------------------------------------------------------
!
       source = coords(1)*ntiles(2)
       tag    = 2000 + coords(2)
       call MPI_RECV(intm_buf(1), count, &
                    MPI_DOUBLE_PRECISION, &
                    source, tag, comm3d, stat, ierr)
       do i = 1, 2*in
        intm_half(i) = intm_buf(i)
       enddo
      endif ! coords(2)
!
!-----------------------------------------------------------------------
!     using theta-summed interior masses, each processor updates its
!     value of local_mass, which is now the total mass contained by all
!     processors having a common value of coords(1)
!-----------------------------------------------------------------------
!
      nhalfcells = 2*(ie-is+1)
      local_mass = intm_half(nhalfcells+1)
!
      deallocate(intm_buf)
!
      return
      end
!=======================================================================
      subroutine sum_over_one_coord(intm_half,local_mass)
!
      use real_prec
      use param
      use grid
      use mpiyes
      use mpipar
!
      implicit NONE
!
      real(rl) :: intm_half(2*in)
      real(rl) :: local_mass, loc_buf
!
      integer source, target, tag, i, l, nhalfcells
!
!-----------------------------------------------------------------------
!     all processes with coords(1) = 0 begin upward cascading sum of
!     interior mass, passing data to all processes with same value of 
!     coords(2)
!-----------------------------------------------------------------------
!
      nhalfcells = 2*(ie-is+1)
!
      if(coords(1) .eq. 0) then
!
       tag = 1000*n1p
       call MPI_SEND(local_mass, 1, &
                     MPI_DOUBLE_PRECISION, &
                     n1p, tag, comm3d, ierr) 
!
      else ! coords(1)
!
       tag = 1000*myid
       call MPI_RECV(loc_buf, 1, &
                     MPI_DOUBLE_PRECISION, &
                     n1m, tag, comm3d, stat, ierr)
!
       do i = 1, nhalfcells+1
        intm_half(i) = intm_half(i) + loc_buf
       enddo
!
       local_mass = intm_half(nhalfcells+1)
!
       tag = 1000*n1p
!
       if(coords(1) .ne. ntiles(1)-1) then
        call MPI_SEND(local_mass, 1, &
                      MPI_DOUBLE_PRECISION, &
                      n1p, tag, comm3d, ierr) 
       endif ! coords(1)
!
      endif ! coords(1)
!
      return
      end
