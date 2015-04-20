c
C     CSR (26th June 2003)
C
C     SUBROUTINE TO INCLUDE A KINEMATIC VISCOSITY INTO THE MOMENTUM
C     AND ENERGY EQUATIONS.   MODIFIED VERSION OF CODE OBTAINED FROM
C     JIM STONE (VIA EMAIL ON THE 26TH JUNE 2003).
C
C     ISSUE THE THE USUAL ZEUS-MP INCLUDE COMMANDS
C
      subroutine kin_visc
#ifdef PHYS_VISC
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      implicit none
C
C     SET UP SUBROUTINE VARIABLES AND ARRAYS
C
      integer  :: i,j,k
      real(rl) :: divv,del,divvs,sqr,exy,eyz,ezx,mu
      real(rl) :: v1st(in,jn,kn),v2st(in,jn,kn),v3st(in,jn,kn)
      real(rl) :: total_diss

      common/dissipation/total_diss
c
c     set value for mu
c
      mu=0.01
      namelist /visc/ mu
C
C     PERFORM EXPLICIT UPDATE FOR VELOCITIES
C     SOURCE STEP FOR V1
C
      do k=ks,ke
         do j=js,je
            do i=is,ie
               v1st(i,j,k) =mu*dt*(v1(i+1,j,k)-2*v1(i,j,k)+
     &              v1(i-1,j,k))/
     &              (dx1a(i))**2
     &              + mu*dt*(v1(i,j+1,k)-2*v1(i,j,k)+v1(i,j-1,k))
     &              /(dx2b(j))**2
     &              + mu*dt*(v1(i,j,k+1)-2*v1(i,j,k)+v1(i,j,k-1))
     &              /(dx3b(k))**2
     &              + (mu/3.0)*dt*(v1(i+1,j,k)-2*v1(i,j,k)+v1(i-1,j,k))
     &              /(dx1a(i))**2
     &              + (mu/3.0)*dt*(v2(i,j+1,k)-v2(i,j,k)-v2(i-1,j+1,k)
     &              +v2(i-1,j,k))/
     &              (dx2a(j)*dx1b(i)) + (mu/3.0)*dt*(v3(i,j,k+1)
     &              -v3(i,j,k)-v3(i-1,j,k+1)
     &              +v3(i-1,j,k))/(dx3a(k)*dx1b(i))
            enddo
         enddo
      enddo
C
C     SOURCE STEP FOR V2
C
      do k=ks,ke
         do j=js,je
            do i=is,ie
               v2st(i,j,k) = mu*dt*(v2(i+1,j,k)-2*v2(i,j,k)
     &              +v2(i-1,j,k))/
     &              (dx1b(i))**2
     &              + mu*dt*(v2(i,j+1,k)-2*v2(i,j,k)+v2(i,j-1,k))/
     &              (dx2a(j))**2
     &              + mu*dt*(v2(i,j,k+1)-2*v2(i,j,k)+v2(i,j,k-1))
     &              /(dx3b(k))**2
     &              + (mu/3.0)*dt*(v2(i,j+1,k)-2*v2(i,j,k)+v2(i,j-1,k))
     &              /(dx2a(j))**2
     &              + (mu/3.0)*dt*(v1(i+1,j,k)-v1(i,j,k)-v1(i+1,j-1,k)
     &              +v1(i,j-1,k))/
     &              (dx2b(j)*dx1a(i)) + (mu/3.0)*dt*(v3(i,j,k+1)
     &              -v3(i,j,k)-v3(i,j-1,k+1)
     &              +v3(i,j-1,k))/(dx3a(k)*dx2b(j))
            enddo
         enddo
      enddo
C
C     SOURCE STEP FOR V3
C
      do k=ks,ke
         do j=js,je
            do i=is,ie
              v3st(i,j,k) = mu*dt*(v3(i+1,j,k)-2*v3(i,j,k)+v3(i-1,j,k))/
     &              (dx1b(i))**2
     &              + mu*dt*(v3(i,j+1,k)-2*v3(i,j,k)+v3(i,j-1,k))
     &              /(dx2b(j))**2
     &              + mu*dt*(v3(i,j,k+1)-2*v3(i,j,k)+v3(i,j,k-1))
     &              /(dx3a(k))**2
     &              + (mu/3.0)*dt*(v3(i,j,k+1)-2*v3(i,j,k)+v3(i,j,k-1))
     &              /(dx3a(k))**2
     &              + (mu/3.0)*dt*(v1(i+1,j,k)-v1(i,j,k)-v1(i+1,j,k-1)
     &              +v1(i,j,k-1))/
     &              (dx3b(k)*dx1a(i)) + ( mu/3.0)*dt*(v2(i,j+1,k)
     &              -v2(i,j,k)-v2(i,j+1,k-1)
     &              +v2(i,j,k-1))/(dx3b(k)*dx2a(j))
           enddo
         enddo
      enddo
C
C     UPDATE VELOCITIES
C
      do k=ks,ke
         do j=js,je
            do i=is,ie
               v1(i,j,k) = v1(i,j,k) + v1st(i,j,k)
               v2(i,j,k) = v2(i,j,k) + v2st(i,j,k)
               v3(i,j,k) = v3(i,j,k) + v3st(i,j,k)
           enddo
         enddo
      enddo
C
C     IMPLEMENT BOUNDARY CONDITIONS FOR DENSITY
C
ccc      call bvalv1(3,3,3,3,3,3,v1)
ccc      call bvalv2(3,3,3,3,3,3,v2)
ccc      call bvalv3(3,3,3,3,3,3,v3)
C
C     DISSIPATION TERM
C
        do k=ks,ke
         do j=js,je
            do i=is,ie
               divvs=((v1(i+1,j,k)-v1(i,j,k))/dx1a(i)
     &              +(v2(i,j+1,k)-v2(i,j,k))/dx2a(j)+(v3(i,j,k+1)
     &              -v3(i,j,k))/dx3a(k))**2
               sqr=((v1(i+1,j,k)-v1(i,j,k))/dx1a(i))**2 +
     &              ((v2(i,j+1,k)-v2(i,j,k))/dx2a(j))**2 +
     &              ((v3(i,j,k+1)-v3(i,j,k))/dx3a(k))**2
               exy=(v1(i,j+1,k)+v1(i+1,j+1,k)-(v1(i,j-1,k)
     &              +v1(i+1,j-1,k)))/
     &              4.0*dx2b(j) +
     &              (v2(i+1,j,k)+v2(i+1,j+1,k)-(v2(i-1,j,k)
     &              +v2(i-1,j+1,k)))/4.0*dx1b(i)
               eyz=(v2(i,j,k+1)+v2(i,j+1,k+1)-(v2(i,j,k-1)
     &              +v2(i,j+1,k-1)))/
     &              4.0*dx3b(k) +
     &              (v3(i,j+1,k)+v3(i,j+1,k+1)-(v3(i,j-1,k)
     &              +v3(i,j-1,k+1)))/4.0*dx2b(j)
               ezx=(v1(i,j,k+1)+v1(i+1,j,k+1)-(v1(i,j,k-1)
     &              +v1(i+1,j,k-1)))/
     &              4.0*dx3b(k) +
     &              (v3(i+1,j,k)+v3(i+1,j,k+1)-(v3(i-1,j,k)
     &              +v3(i-1,j,k+1)))/4.0*dx1b(i)
               e(i,j,k)=e(i,j,k)+dt*(mu*(exy**2+eyz**2+ezx**2)
     &              +2.0*mu*sqr
     &              -0.6667*mu*divvs)
               total_diss=total_diss+(mu*(exy**2+eyz**2+ezx**2)
     &              +2.0*mu*sqr-0.6667*mu*divvs)*
     &              dx1b(i)*dx2b(j)*dx3b(k)*dt
           enddo
         enddo
      enddo
C
C     IMPLEMENT BOUNDARY CONDITION FOR INTERNAL ENERGY
C     Chris took this out of his copy to fix things for multiprocessor runs
C      call bvale(3,3,3,3,3,3,e)
C
#endif /* PHYS_VISC */
      return
      end

