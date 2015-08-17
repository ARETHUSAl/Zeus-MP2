!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               S R C S T E P               \\\\\\\\\!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
!=======================================================================
!
subroutine srcstep
!
!    jms:zeus2d.srcstep <------------------------ source step controller
!                                                          october, 1987
!
!    written by: Jim Stone
!    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
!    modified 2: Spring, 1989 by Jim Stone; rewritten
!    modified 3: February, 1990 by David Clarke; incorporated into
!                ZEUS3D
!    modified 4: July, 1990 by David Clarke; because current densities
!                are not needed to compute source terms in the new CT
!                algorithm (MOCCT), workers can be used to store the
!                momenta, thereby saving redundant computations in STV1
!                and STV2.
!    modified 5: June 1992, by David Clarke; added the total energy
!                option originally designed by Byung-IL Jun.
!    modified 6: RAF, 3/27/96, completely rewritten for ZEUS-MP.
!    modified 7: RAF, 1/2/97, added radiation force, radiation
!                diffusion, and radiation-matter coupling terms.
!    modified 8: RAF, 1/22/97, moved modules into driver routines.
!    modified 9: RAF, 2/18/97, added PF and N-R timestep controllers.
!    modified 10: PSLi, 12/30/99, added subcycle of artificial viscosity.
!
!
!  PURPOSE: Controls the update of velocities (v1, v2, v3) and internal
!  energy (e) from source terms in the equation of motion and energy
!  equations respectively.
!
!  LOCAL VARIABLES:      
!    w3da     scratch 1-momentum denisty
!    w3db     scratch 2-momentum denisty
!    w3dc     scratch 3-momentum denisty
!    w3dd     scratch 1-velocity ; scratch density (pdv)
!    w3de     scratch 2-velocity ; scratch e/d     (pdv)
!    w3df     scratch 3-velocity
!    j1       1-current density
!    j2       2-current density
!    j3       3-current density
!
!  EXTERNALS:
!    AVISC_D  ,  FORCES_D  , PDV_D
!
!-----------------------------------------------------------------------
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
      use cons
!
      implicit NONE
!
      real(rl) :: subdt, etot, etot_glb
!
      integer  :: i, j, k, index, n
!
!-----------------------------------------------------------------------
!     EQUATION OF STATE
!
!     In this release, gamma-law ideal gas and isothermal equations of
!     state are included and supported by the EOS parameter LEOS=1.
!     User-added EOS modules will require LEOS parameters of 2 and 
!     higher.
!-----------------------------------------------------------------------
!
      call eos_d
!
!-----------------------------------------------------------------------
! FORCES
!
! Thermal and MHD Pressure, Gravity, and Rotational Pseudo-Forces
! (including gravitational point mass)
!
! Routine "forces" updates the three velocity components.  The arrays
! v1, v2, and v3 save the old values, while w3dd, w3de, w3df get the 
! new ones.
!
      if(lrad .gt. 0) call opac_d
!
      call forces_d (v1,v2,v3,w3dd,w3de,w3df)
!
!       subroutine forces_d
!     1 (v1old, v2old, v3old, v1new, v2new, v3new)
!
!.......................................................................
!
! ARTIFICIAL VISCOSITY
!
!  Update velocity and e from the artificial viscosity source terms.
!  We need 1 layer of updated boundary data for the velocity.
!  We use just the "m" layer of d and, unless ISO is defined, e -- they
!  should be up to date already, since we used them in "forces" but did
!  not update them.  We don't need e BVs with no "linear" viscosity.
!
!  Arrays w3dd, w3de, w3df save the old velcoity values, while v1, v2, 
!  v3 get the updated ones.
!
!  The artificial viscosity routine must also compute momentum densities
!  w3da, w3db, w3dc from v1, v2, v3 for use in the transport step.
!
!PS
!  Subcycle of artificial viscosity calculation
!
      if(xsubav) then
       index=0
       subdt=dt
       avisc_dt=courno / (SQRT(dtqqi2)+tiny)
       if(dtqqi2.eq.0.0) avisc_dt=dt
          buf_in(1) = avisc_dt
          call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                            , MPI_2DOUBLE_PRECISION &
                            , MPI_MINLOC, comm3d, ierr)
          avisc_dt  =   buf_out(1)
       do while (subdt .gt. 0.0)
          if(subdt.gt.avisc_dt) then
             subdt=subdt-avisc_dt
          else
             avisc_dt=subdt
             subdt=-1.0
          endif
          index=index+1
          if(mod(index,2).gt.0) then
             call avisc_d (w3dd,w3de,w3df,v1,v2,v3,w3da,w3db,w3dc)
          else
             call avisc_d (v1,v2,v3,w3dd,w3de,w3df,w3da,w3db,w3dc)
          endif
          avisc_dt=courno / (SQRT(dtqqi2)+tiny)
          buf_in(1) = avisc_dt
          call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1 &
                            , MPI_2DOUBLE_PRECISION &
                            , MPI_MINLOC, comm3d, ierr)
          avisc_dt  =   buf_out(1)
       enddo
!
!  update velocity arrays if necessary.
!
       if(mod(index,2).eq.0) then
          do k=ks,ke
             do j=js,je
                do i=is,ie
                   v1(i,j,k)=w3dd(i,j,k)
                   v2(i,j,k)=w3de(i,j,k)
                   v3(i,j,k)=w3df(i,j,k)
                enddo
             enddo
          enddo
       endif
      else ! xsubav
       call avisc_d (w3dd,w3de,w3df,v1,v2,v3,w3da,w3db,w3dc)
      endif ! xsubav
!
!
!      subroutine avisc_d
!     1            (v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
!
!......................................................................
!     IMPLICIT GREY FLUX-LIMITED DIFFUSION SOLVER
!......................................................................
!
      if(lrad .ne. 0) call rad_solve
!
!......................................................................
! No radiation, not isothermal.
!
!  COMPRESSIONAL WORK TERM (PdV).
!
!  Finally, update the energy with the PdV source term
!  only if the internal energy equation is being solved.
!
!  We need just 1 "p" layer of updated boundary data for v1, v2, and v3,
!  but none for d and e.  
!
!  Routine pdv also saves the density and e/d in (w3dd,w3de) for use
!  in the transport step.  This is why it is being called even when
!  TOTAL_ENERGY is defined.
!
!  NOTE: PDV is not called if radiation transport (lrad > 0) is used
!        because it is assumed that the pdV term is incorporated in
!        the radiation step (as it is in the diffusion solver
!        included here).  A user-supplied radiation option which
!        neglects the pdV step will require an adjustment to the
!        IF condition which follows:
!
      if((xiso .eqv. .false.) .and. (lrad .eq. 0)) then
!       if(xtotnrg) call eos_d
       call pdv_d (w3dd ,w3de )
      endif ! xiso
!
      return
end subroutine srcstep
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               S R C S T E P               \\\\\\\\\!
!=======================================================================
!
