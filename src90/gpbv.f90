!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine gpbv
!
!  Written by PSLi (12/4/99)
!
!  PURPOSE: Calculate the gravitational potential at the boundary
!           surfaces contributed from monopole and quadrupole moments
!           of a mass distribution.
!           Boundary potentals are calculated if the flags niis(3),
!           nois(3), nijs(3), nojs(3), niks(3), or noks(3) = 3,
!           respectively.
!
!  EXTERNALS: [none]
!
!  LOCALS:
!  dm(i,j,k)   is mass contained in zone i,j,k
!  lqm         is the local quadrupole moment (qm)
!  lcx,lcy,lcz coords of mass center at local processor (cx,cy,cz)
!  ltm         total mass at local processor (tm)
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use bndry
      use root
      use scratch
      use mpiyes
      use mpipar
!
      implicit NONE
!
      real(rl) ::     errmax
      integer  :: n
!
      real(rl) ::  dm(in,jn,kn),qm(6,nprocs_w),tm(nprocs_w), &
                   cx(nprocs_w), &
                   cy(nprocs_w),cz(nprocs_w),r1,dx,dy,dz,dx2,dy2,dz2, &
                   xdm,ydm,zdm,ltm,lcx,lcy,lcz,lqm(6),cthe,cphi,sthe, &
                   sphi
!
      integer  i,j,k,m,ism1,jsm1,ksm1,iep1,jep1,kep1
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////////
!=======================================================================
!  This routine calculate the monopole and quadrupole moments of a mass
!  distribution and the gravitational potential at the boundary surfaces
!  for solving Poisson equation.
!
!  Positive potential for the sign convention of ZeusMP.
!
      do i=1,6
        do m=1,nprocs_w
          qm(i,m) = 0.0
        enddo
        lqm(i) = 0.0
      enddo
      do k=1,kn
        do j=1,jn
          gpiib(j,k,1)=0.0
          gpoib(j,k,1)=0.0
        enddo
      enddo
      do k=1,kn
        do i=1,in
          gpijb(i,k,1)=0.0
          gpojb(i,k,1)=0.0
        enddo
      enddo
      do j=1,jn
        do i=1,in
          gpikb(i,j,1)=0.0
          gpokb(i,j,1)=0.0
        enddo
      enddo
      xdm=0.0
      ydm=0.0
      zdm=0.0
      ism1=is-1
      jsm1=js-1
      ksm1=ks-1
      iep1=ie+1
      jep1=je+1
      kep1=ke+1
!
!  Determine the center of mass.
!
      do i=1,nprocs_w
        tm(i)=tiny
        cx(i)=0.0
        cy(i)=0.0
        cz(i)=0.0
      enddo
      ltm=0.0
      lcx=0.0
      lcy=0.0
      lcz=0.0
!
      do k=ks,ke
       if(lgeom .eq. 2 .or. lgeom .eq. 3) then
        cphi=cos(x3b(k))
        sphi=sin(x3b(k))
       endif ! CYL or SPHERE
	do j=js,je
         if(lgeom .eq. 3) then
          cthe=cos(x2b(j))
          sthe=sin(x2b(j))
         endif ! SPHERE
          do i=is,ie
!            dm(i,j,k) = g2b(i)*g31b(i)*g32b(j)*dx1a(i)*
!     &                  dx2a(j)*dx3a(k)*d(i,j,k)
            dm(i,j,k) = d(i,j,k)*dvl1a(i)*dvl2a(j)*dvl3a(k)
            ltm=ltm+dm(i,j,k)
           if(lgeom .eq. 1) then
            xdm=xdm+x1b(i)*dm(i,j,k)
            ydm=ydm+x2b(j)*dm(i,j,k)
            zdm=zdm+x3b(k)*dm(i,j,k)
           endif ! CART
           if(lgeom .eq. 2) then
            xdm=xdm+x2b(j)*cphi*dm(i,j,k)
            ydm=ydm+x2b(j)*sphi*dm(i,j,k)
            zdm=zdm+x1b(i)*dm(i,j,k)
           endif ! CYL
           if(lgeom .eq. 3) then
            xdm=xdm+x1b(i)*sthe*cphi*dm(i,j,k)
            ydm=ydm+x1b(i)*sthe*sphi*dm(i,j,k)
            zdm=zdm+x1b(i)*cthe*dm(i,j,k)
           endif ! SPHERE
          enddo
        enddo
      enddo
!
      lcx=xdm/ltm
      lcy=ydm/ltm
      lcz=zdm/ltm
!
!  Compute the 2nd order terms in the expansion.
!
      do k=ks,ke
       if(lgeom .eq. 2 .or. lgeom .eq. 3) then
        sphi=sin(x3b(k))
        cphi=cos(x3b(k))
       endif ! CYL or SPHERE
	do j=js,je
         if(lgeom .eq. 3) then
          cthe=cos(x2b(j))
          sthe=sin(x2b(j))
         endif ! SPHERE
          do i=is,ie
           if(lgeom .eq. 1) then
            dx=x1b(i)-lcx
            dy=x2b(j)-lcy
            dz=x3b(k)-lcz
           endif ! CART
           if(lgeom .eq. 2) then
            dx=x2b(j)*cphi-lcx
            dy=x2b(j)*sphi-lcy
            dz=x1b(i)-lcz
           endif ! CYL
           if(lgeom .eq. 3) then
            dx=x1b(i)*sthe*cphi-lcx
            dy=x1b(i)*sthe*sphi-lcy
            dz=x1b(i)*cthe-lcz
           endif ! SPHERE
            dx2=dx*dx
            dy2=dy*dy
            dz2=dz*dz
            lqm(1)=lqm(1)+(2.*dx2-dy2-dz2)*dm(i,j,k)
            lqm(2)=lqm(2)+(2.*dy2-dx2-dz2)*dm(i,j,k)
            lqm(3)=lqm(3)+(2.*dz2-dx2-dy2)*dm(i,j,k)
            lqm(4)=lqm(4)+3.*dx*dy*dm(i,j,k)
            lqm(5)=lqm(5)+3.*dx*dz*dm(i,j,k)
            lqm(6)=lqm(6)+3.*dy*dz*dm(i,j,k)
          enddo
        enddo
      enddo
!
!  Sum up qm from all processors.
!
      call MPI_Gather(ltm,1,MPI_DOUBLE_PRECISION, &
             tm,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcx,1,MPI_DOUBLE_PRECISION, &
             cx,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcy,1,MPI_DOUBLE_PRECISION, &
             cy,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcz,1,MPI_DOUBLE_PRECISION, &
             cz,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lqm,6,MPI_DOUBLE_PRECISION, &
             qm,6,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(tm,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cx,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cy,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cz,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(qm,nprocs_w*6,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
!
!  Now compute the potential along each boundary surface, starting
!  with the inner i surface.
!
      if (niis(3) .eq. 3) then
        do m=1,nprocs_w
          do k=ks,ke
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif ! CYL OR SPHERE
            do j=js,je
             if(lgeom .eq. 3) then
              cthe=cos(x2b(j))
              sthe=sin(x2b(j))
             endif ! SPHERE
             if(lgeom .eq. 1) then
              dx=x1b(ism1)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(k)-cz(m)
             endif ! CART
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(ism1)-cz(m)
             endif ! CYL
             if(lgeom .eq. 3) then
              dx=x1b(ism1)*sthe*cphi-cx(m)
              dy=x1b(ism1)*sthe*sphi-cy(m)
              dz=x1b(ism1)*cthe-cz(m)
             endif ! SPHERE
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpiib(j,k,1) = gpiib(j,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
!  Initial solution along outer i surface
!
      if (nois(3) .eq. 3) then
        do m=1,nprocs_w
          do k=ks,ke
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do j=js,je
             if(lgeom .eq. 3) then
              cthe=cos(x2b(j))
              sthe=sin(x2b(j))
             endif
             if(lgeom .eq. 1) then
              dx=x1b(iep1)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(iep1)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(iep1)*sthe*cphi-cx(m)
              dy=x1b(iep1)*sthe*sphi-cy(m)
              dz=x1b(iep1)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpoib(j,k,1) = gpoib(j,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
!  Initial solution along inner j surface
!
      if (nijs(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          cthe=cos(x2b(jsm1))
          sthe=sin(x2b(jsm1))
         endif
          do k=ks,ke           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(jsm1)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(jsm1)*cphi-cx(m)
              dy=x2b(jsm1)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpijb(i,k,1) = gpijb(i,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
!  Initial solution along outer j surface
!
      if (nojs(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          cthe=cos(x2b(jep1))
          sthe=sin(x2b(jep1))
         endif
          do k=ks,ke           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(jep1)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(jep1)*cphi-cx(m)
              dy=x2b(jep1)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpojb(i,k,1) = gpojb(i,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
!  Initial solution along inner k surface
!
      if (niks(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          sphi=sin(x3b(ksm1))
          cphi=cos(x3b(ksm1))
         endif
          do j=js,je           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            cthe=cos(x2b(j))
            sthe=sin(x2b(j))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(ksm1)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpikb(i,j,1) = gpikb(i,j,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
!  Initial solution along outer k surface
!
      if (noks(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          sphi=sin(x3b(kep1))
          cphi=cos(x3b(kep1))
         endif
          do j=js,je           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            cthe=cos(x2b(j))
            sthe=sin(x2b(j))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(kep1)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpokb(i,j,1) = gpokb(i,j,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+ &
                             qm(2,m)*dy2+qm(3,m)*dz2)+ &
                             qm(4,m)*dx*dy+qm(5,m)*dx*dz+ &
                             qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
!
      return
      end
