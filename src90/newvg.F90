!=======================================================================
!
      subroutine newvg
!
!  PURPOSE: Computes grid velocities at current timestep to be used in
!  updating the grid variables.  Velocoties are claculated over the
!  range i=is-2,ie+2 and j=js-2,je+2.  The method used to compute the
!  grid velocities depends on the flag igcon and the sign of x1fac,x2fac
!        igcon = 1 gives "lagrangean" tracking in x1 [x2] lines
!        igcon = 2 for input grid boundary speeds
!                vg1(io) = x1fac * central soundspeed
!                vg2(jo) = x2fac * central soundspeed
!        igcon = 3 for uniform translation
!
!  Modified 1:
!     04-05-04 by JCH; modified lagrangean tracking for multi-D: use
!                      V_g(i) = Sum_(j)[rho(i,j)*V1(i,j)] /
!                               Sum_(j)[rho(i,j)]
!                      (i.e. specific radial momentum / specific mass)
!  Modified 2: enclosed MPI call inside of "ifdef MPI_USED" cpp 
!              constructs (JHayes; 11/21/05)
!-----------------------------------------------------------------------
      use config
      use param
      use grid
      use field
      use bndry
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
!
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, dest, l, &
                  source, ncheck, ic3 !, nhalo, nlines(3)
      real(rl) :: qa, qb, rhosum(in), rhovsum(in), buf1(in), &
                          buf2(in), vexpan(nprocs_w), vexpanmax, &
                          vexpan2, ic(nprocs_w), ic2, amtemp
!      real(rl), dimension(128) :: dtable,Ttable,vrtable,rtable
!      common /inflow/ vexpanmax,amtemp,dtable,Ttable,vrtable,rtable,
!     .                nhalo, nlines
      real(rl), dimension(128) :: dwind, twind, vwind
      common /inflow/ vexpanmax,dwind,twind,vwind
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////
!=======================================================================
!
!  igcon =  1: "Lagrangean" tracking of grid lines
!
      if (igcon .eq. 1 .or. igcon .eq. 11) then
       if(niis(1) .eq. 0) then
        ibeg = is-2
       else
        ibeg = is+1
       endif
       if(nois(1) .eq. 0) then
        iend = ie+2
       else
        iend = ie
       endif
!
      if (nx1z .gt. 1 .and. x1fac .lt. 0.0) then
!
! -- zero arrays to be summed
!
       do i = ibeg, iend
        rhovsum(i) = 0.0
        rhosum (i) = 0.0
       enddo
!
! -- compute total specific (radial) momentum and total
! -- specific mass at each radial shell
!
       do k = ks, ke
        do j = js, je
         do i = ibeg, iend
          rhovsum(i) = rhovsum(i) + d(i,j,k)*v1(i,j,k)
          rhosum (i) = rhosum (i) + d(i,j,k)
         enddo
        enddo
       enddo
!
!-----------------------------------------------------------------------
! -- take ratio to compute radial grid velocity; do the
! -- easy case (no parallel decomp along theta) first
!
       if(ntiles(2) .eq. 1) then
        if(igcon .eq. 1) then
         do i = ibeg, iend
          vg1(i) = rhovsum(i)/rhosum(i)
         enddo
        endif
!
!-----------------------------------------------------------------------
!     this a test option for the MOVING_GRID_TEST problem generator
!-----------------------------------------------------------------------
!
        if(igcon .eq. 11) then
         do i = ibeg, iend
          vg1(i) =  -0.1*float(i-is)/float(ie-is+1)
         enddo
        endif
       endif
!-----------------------------------------------------------------------
!
! -- harder case: parallel decomposition along theta axis
!
#ifdef MPI_USED
       if(ntiles(2) .gt. 1) then
        if(igcon .eq. 1) then
         call sum_along_process_rows(rhovsum)
         call sum_along_process_rows(rhosum)
         do i = ibeg, iend
          vg1(i) = rhovsum(i)/rhosum(i)
         enddo
        endif
!
        if(igcon .eq. 11) then
         do i = ibeg, iend
          vg1(i) =  -0.1*float(i-is)/float(ie-is+1)
         enddo
        endif
       endif ! ntiles(2) > 1
#endif /* MPI_USED */
!
!-----------------------------------------------------------------------
!
       if(niis(1) .ne. 0) then
        vg1(is  ) = 0.0
        vg1(is-1) = -vg1(is+1)
        vg1(is-2) = -vg1(is+2)
       endif
       if(nois(1) .ne. 0) then
        vg1(ie+1) = 0.0
        vg1(ie+2) = -vg1(ie  )
       endif
      endif ! nx1z
!
       if (nx2z .gt. 1 .and. x2fac .lt. 0.0) then
        if(nijs(1) .eq. 0) then
         jbeg = js-2
        else
         jbeg = js+1
        endif
        if(nojs(1) .eq. 0) then
         jend = je+2
        else
         jend = je
        endif
        do 110 j=jbeg,jend
          vg2(j) = v2(is,j,ks)
110     continue
        if(nijs(1) .ne. 0) then
         vg2(js  ) = 0.0
         vg2(js-1) = -vg2(js+1)
         vg2(js-2) = -vg2(js+2)
        endif
        if(nojs(1) .ne. 0) then
         vg2(je+1) = 0.0
         vg2(je+2) = -vg2(je  )
        endif
       endif ! nx2z
!
       if (nx3z .gt. 1 .and. x3fac .lt. 0.0) then
        if(niks(1) .eq. 0) then
         kbeg = ks-2
        else
         kbeg = ks+1
        endif
        if(noks(1) .eq. 0) then
         kend = ke+2
        else
         kend = ke
        endif
        do 120 k = kbeg, kend
          vg3(k) = v3(is,js,k)
120     continue
        if(niks(1) .ne. 0) then
         vg3(ks  ) = 0.0
         vg3(ks-1) = -vg3(ks+1)
         vg3(ks-2) = -vg3(ks+2)
        endif
        if(noks(1) .ne. 0) then
         vg3(ke+1) = 0.0
         vg3(ke+2) = -vg3(ke  )
        endif
       endif ! nx3z
       return
      endif ! igcon
!
!  igcon=2:  qa=central sound speed; vg is computed as a
!  linear function of x from x=0 to x(outer boundary).
!
      if (igcon .eq. 2) then
       qa = sqrt(gamma*(gamma-1.0)*e(is,js,ks)/d(is,js,ks))
       if (x1fac .ne. 0.0) then
        qb = x1fac*qa/x1a(ie)
        do 200 i=is+1,ie
          vg1(i)=-qb*x1a(i)
200     continue
        vg1(is  ) = 0.0
        vg1(is-1) = -vg1(is+1)
        vg1(is-2) = -vg1(is+2)
        vg1(ie+1) = 0.0
        vg1(ie+2) = -vg1(ie  )
       endif ! x1fac
!
       if (x2fac .ne. 0.0) then
        qb = x2fac*qa/x2a(je)
        do 210 j=js+1,je
          vg2(j)=-qb*x2a(j)
210     continue
        vg2(js  ) = 0.0
        vg2(js-1) = -vg2(js+1)
        vg2(js-2) = -vg2(js+2)
        vg2(je+1) = 0.0
        vg2(je+2) = -vg2(je  )
       endif ! x2fac
!
       if (x3fac .ne. 0.0) then
        qb = x3fac*qa/x3a(ke)
        do 220 k=ks+1,ke
          vg3(k)=-qb*x3a(k)
220     continue
        vg3(ks  ) = 0.0
        vg3(ks-1) = -vg3(ks+1)
        vg3(ks-2) = -vg3(ks+2)
        vg3(ke+1) = 0.0
        vg3(ke+2) = -vg3(ke  )
       endif ! x3fac
      return
      endif
!
!     igcon = 3.  Uniform translation of the grid at a velocity 
!     x1[2][3]fac.
!
      if (igcon .eq. 3) then
       do 300 i = is-2, ie+2
        vg1(i) = x1fac
300    continue
       do 310 j = js-2, je+2
        vg2(j) = x2fac
310    continue
       do 320 k = ks-2, ke+2
        vg3(k) = x3fac
320    continue
       return
      endif
!
!----------------------------------------------------------------
!------------------------- Begin igcon = 5 control --------------
!----------------------------------------------------------------
      if (igcon .eq. 5) then
      nreq = 0
      do i=1,nprocs_w
      vexpan(i) = 0.
      ic(i)     = 4.
      enddo
      vexpanmax = 0.
      vexpan2   = 0.
      ic3       = 4
      ncheck = int(0.1*(ie-is))
       do k=ks,ke
       do j=js,je
       do i=ie,ie-ncheck,-1
         if ( v1(i,j,k).gt.vexpan(myid_w+1)) then
           vexpan(myid_w+1) = 3.0 * v1(i,j,k)
           ic(myid_w+1)     = i
         endif
!         if (nhy.eq.499) print*,myid_w+1,v1(i,j,k)
       enddo
       enddo
       enddo
#ifdef MPI_USED
!       vexpan2 = vexpan(myid_w+1)
!       ic2     = ic(myid_w+1)
!       nreq = nreq + 1
!       call MPI_Allgather(vexpan2,1,MPI_FLOAT
!     &                   ,vexpan ,1,MPI_FLOAT
!     &                   ,MPI_COMM_WORLD,ierr)
!       call MPI_Allgather(ic2,1,MPI_FLOAT
!     &                   ,ic ,1,MPI_FLOAT
!     &                   ,MPI_COMM_WORLD,ierr)
!       if(nreq .ne. 0) then
!       call mpi_waitall(nreq, req, stat, ierr)
!       nreq = 0
!       endif
#endif /* MPI_USED */
       do i=1,nprocs_w
       if ( vexpan(i).gt.vexpanmax) then
         vexpanmax = vexpan(i)
         ic3       = int(ic(i))
       endif
       enddo
!      print*,"vexpanmax is: ",vexpanmax
       do i=is-2,ie+2
        vg1(i)= vexpanmax * x1a(i) / x1a(ic3) 
       enddo
      endif
!----------------------------------------------------------------
!------------------------- End igcon = 5 control ----------------
!----------------------------------------------------------------
!
      return
      end
!
      subroutine sum_along_process_rows(summed_vec)
#ifdef MPI_USED
!
      use param
      use grid
      use field
      use bndry
      use root
      use mpiyes
      use mpipar
!
      implicit NONE
!
      integer  :: l, nlev, src, dest, tag, pid, npe, ilev, stride, i, &
                  disp, rat
      real(rl) :: summed_vec(in), buf1(in)
!
      if(coords(2) .gt. 0) then
       dest = coords(1)*ntiles(2)
       tag  = myid*4096
       call MPI_SEND(summed_vec(1), in, MPI_FLOAT, &
                     dest, tag, comm3d, ierr)
      endif ! coords(2) > 0
!
      if(coords(2) .eq. 0) then
       do l = 1, ntiles(2)-1
        src = myid + l
        tag = src*4096
        call MPI_RECV(buf1(1), in, MPI_FLOAT, &
                      src, tag, comm3d, stat, ierr)
        do i = is-2, ie+2
         summed_vec(i) = summed_vec(i) + buf1(i)
        enddo
       enddo ! l
!
       do l = 1, ntiles(2)-1
        dest = myid + l
        tag  = dest*4096
        call MPI_SEND(summed_vec(1), in, MPI_FLOAT, &
                      dest, tag, comm3d, ierr)
       enddo ! l
      endif ! coords(2) = 0
      if(coords(2) .gt. 0) then
       src = coords(1)*ntiles(2)
       tag = myid*4096
       call MPI_RECV(buf1(1), in, MPI_FLOAT, &
                     src, tag, comm3d, stat, ierr)
       do i = is-2, ie+2
        summed_vec(i) = buf1(i)
       enddo
      endif ! coords(2) .gt. 0
#endif /* MPI_USED */
!
      return
      end
