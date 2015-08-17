!=======================================================================
!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!     Purpose: File writer for HDF5 viz data files
!
!     Written by: John Hayes, February 2006
!
      subroutine hdfall(filename)
!
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use scratch
      use cons
      use chem
      use mpiyes
      use mpipar
!
      use hdf5
!
      implicit NONE
!
!-----------------------------------------------------------------------
!     ZEUS-MP - specific file and data descriptors, arrays, etc.
!-----------------------------------------------------------------------
!
      character*16 :: filename
      character*56 :: filename2, incrf
      character*16 :: coordsys
      character*32 :: string,myspecies(12),species_names(12)
!
      integer      :: i,j,k,indx,kp1,l
!
      real(rl4)    :: icoord(ie-is+1)
      real(rl4)    :: jcoord(je-js+1)
      real(rl4)    :: kcoord(ke-ks+1)
!
      real(rl4)    :: tval
      real(rl4)    :: data(ie-is+1,je-js+1,ke-ks+1), tmp
!
!-----------------------------------------------------------------------
!     hdf5-specific parameters, identifiers, etc.
!-----------------------------------------------------------------------
!
      integer        :: rank, error
      integer(hid_t) :: file_id
      integer(hsize_t), dimension(7) :: dims
!
!-----------------------------------------------------------------------
!     Initialize FORTRAN interface.
!-----------------------------------------------------------------------
!
!JMS  Write restart dumps out to a directory to match Enzo for YT.
      incrf = filename(13:)
      CALL system('mkdir -p '//'DD'//trim(incrf))
      filename2 = 'DD'//trim(incrf)//'/'//trim(filename)
      CALL h5open_f (error)
!
!-----------------------------------------------------------------------
!     Create a new file using default properties.
!-----------------------------------------------------------------------
!
      CALL h5fcreate_f(trim(filename2), H5F_ACC_TRUNC_F, file_id, error)
!
!-----------------------------------------------------------------------
!     Create/Write datasets
!-----------------------------------------------------------------------
!
! --- evolution time
!
      RANK = 1
!
      tval      = real(time)
!
      dims(1  ) = 1
      dims(2:7) = 0
!
      call write_viz(file_id,rank,dims,"time",tval)
!
! --- Coordinate arrays
!
      do i = is, ie
       icoord(i-is+1) = x1b(i)
      enddo
!
      dims(1  ) = ie-is+1
      dims(2:7) = 0
!
      call write_viz(file_id,rank,dims,"i_coord",icoord)
!
      do j = js, je
       jcoord(j-js+1) = x2b(j)
      enddo
!
      dims(1  ) = je-js+1
      dims(2:7) = 0
!
      call write_viz(file_id,rank,dims,"j_coord",jcoord)
!
      do k = ks, ke
       kcoord(k-ks+1) = x3b(k)
      enddo
!
      dims(1  ) = ke-ks+1
      dims(2:7) = 0
!
      call write_viz(file_id,rank,dims,"k_coord",kcoord)
!
! --- Species names
!
!      species_names(1) = 'H'
!      species_names(2) = 'H^+'
!      species_names(3) = 'e^-'
!      species_names(4) = 'He'
!      species_names(5) = 'He^+'
!      species_names(6) = 'He^{2+}'
!      species_names(7) = 'H^-'
!      species_names(8) = 'H_2'
!      species_names(9) = 'H_2^+'
!      species_names(10) = 'HD'
!      species_names(11) = 'HD2'
!      species_names(12) = 'HD3'
!      dims(1  ) = 12
!      dims(2:7) = 0
!      call write_viz(file_id,rank,dims,"species_names",species_names)
!
! --- Field arrays
!
      RANK    = 3
      dims(1  ) = ie-is+1
      dims(2  ) = je-js+1
      dims(3  ) = ke-ks+1
      dims(4:7) = 0
!
! --- i_coord
!
!      !write(*,*)'ks,ke=', ks,ke, ' js,je=',js,je,' is,ie=',is,ie
!      do k=ks,ke
!       do j=js,je
!        do i=is,ie
!         data(i-is+1,j-js+1,k-ks+1) = real(icoord(i-is+1))
!        enddo
!       enddo
!      enddo
!!
!      call write_viz(file_id,rank,dims,"i_coord",data)
!!
!! --- j_coord
!!
!      do k=ks,ke
!       do j=js,je
!        do i=is,ie
!         data(i-is+1,j-js+1,k-ks+1) = real(jcoord(j-js+1))
!        enddo
!       enddo
!      enddo
!!
!      call write_viz(file_id,rank,dims,"j_coord",data)
!!
!! --- k_coord
!!
!      do k=ks,ke
!       do j=js,je
!        do i=is,ie
!         data(i-is+1,j-js+1,k-ks+1) = real(kcoord(k-ks+1))
!        enddo
!       enddo
!      enddo
!
!      call write_viz(file_id,rank,dims,"k_coord",data)
!
! --- 1-Velocity
!
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v1(i,j,k)+v1(i+1,j,k)))
        enddo
       enddo
      enddo
!
      call write_viz(file_id,rank,dims,"i_velocity",data)
!
! --- 2-Velocity
!
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v2(i,j,k)+v2(i,j+1,k)))
        enddo
       enddo
      enddo
!
      call write_viz(file_id,rank,dims,"j_velocity",data)
!
! --- 3-Velocity
!
      do k=ks,ke
       if(ldimen .eq. 3) then
        kp1 = k+1
       else
        kp1 = k
       endif
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v3(i,j,k)+v3(i,j,kp1)))
        enddo
       enddo
      enddo
!
      call write_viz(file_id,rank,dims,"k_velocity",data)
!
      if(XMHD) then
!
! --- 1-B Field
!
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b1(i,j,k)+b1(i+1,j,k)))
         enddo
        enddo
       enddo
!
       call write_viz(file_id,rank,dims,"i_mag_field",data)
!
! --- 2-B field
!
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b2(i,j,k)+b2(i,j+1,k)))
         enddo
        enddo
       enddo
!
       call write_viz(file_id,rank,dims,"j_mag_field",data)
!
! --- 3-B Field
!
       do k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b3(i,j,k)+b3(i,j,kp1)))
         enddo
        enddo
       enddo
!
       call write_viz(file_id,rank,dims,"k_mag_field",data)
      endif ! XMHD
!
! --- Gas density
!
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(d(i,j,k))
        enddo
       enddo
      enddo
!
      call write_viz(file_id,rank,dims,"Density",data)
!
! --- Gas energy
!
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(e(i,j,k))
        enddo
       enddo
      enddo
!
      call write_viz(file_id,rank,dims,"GasEnergy",data)
!
! --- radiation energy
!
      if(lrad .gt. 0) then
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(er(i,j,k))
         enddo
        enddo
       enddo
!
       call write_viz(file_id,rank,dims,"rad_energy",data)
      endif ! lrad
!
! --- Gas Temperature
! 
      if (lrad.ne.0) then 
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(tgas(i,j,k))
         enddo
        enddo
       enddo
       call write_viz(file_id,rank,dims,"Temperature",data)
!      write(*,*) tgas(1,1,1)
       endif !lrad
!
! --- Number Density
! 
      if (nspec .ge. 0) then
       do k=ks,ke
        do j=js,je
         do i=is,ie
          tmp = d(i,j,k)/mh
          data(i-is+1,j-js+1,k-ks+1) = real(tmp)
         enddo
        enddo
       enddo
       call write_viz(file_id,rank,dims,"Number_Density",data)
      endif
!
!      if (nspec .ge. 2 .and. nspec .lt. 6) then
!       do k=ks,ke
!        do j=js,je
!         do i=is,ie
!          tmp = (abun(i,j,k,1) + abun(i,j,k,2)) * d(i,j,k)/mh
!          data(i-is+1,j-js+1,k-ks+1) = real(tmp)
!         enddo
!        enddo
!       enddo
!       call write_viz(file_id,rank,dims,"number_density",data)
!      endif
!c
!      if (nspec .ge. 6 .and. nspec .lt. 9) then
!       do k=ks,ke
!        do j=js,je
!         do i=is,ie
!          tmp = (abun(i,j,k,1) + abun(i,j,k,2) + (abun(i,j,k,4)  + 
!     .  abun(i,j,k,5) + abun(i,j,k,6))* qrt) * d(i,j,k)/mh
!          data(i-is+1,j-js+1,k-ks+1) = real(tmp)
!         enddo
!        enddo
!       enddo
!       call write_viz(file_id,rank,dims,"number_density",data)
!      endif
!c
!      if (nspec .ge. 9) then
!       do k=ks,ke
!        do j=js,je
!         do i=is,ie
!          tmp = (abun(i,j,k,1) + abun(i,j,k,2) +  abun(i,j,k,3) +
!     .  abun(i,j,k,4)*qrt + abun(i,j,k,5)*qrt + abun(i,j,k,6)*qrt + 
!     .  abun(i,j,k,7)     + abun(i,j,k,8)*haf + abun(i,j,k,9)*haf)*
!     .  d(i,j,k)/mh
!          data(i-is+1,j-js+1,k-ks+1) = real(tmp)
!         enddo
!        enddo
!       enddo
!       call write_viz(file_id,rank,dims,"number_density",data)
!      endif
!c
!
! --- Species Abundances
! 
      if (nspec .gt. 1) then
      myspecies(1) = 'HI_Density'
      myspecies(2) = 'HII_Density'
      myspecies(3) = 'Electron_Density'
      myspecies(4) = 'HeI_Density'
      myspecies(5) = 'HeII_Density'
      myspecies(6) = 'HeIII_Density'
      myspecies(7) = 'HM_Density'
      myspecies(8) = 'H2I_Density'
      myspecies(9) = 'H2II_Density'
      myspecies(10) = 'DI_Density'
      myspecies(11) = 'DII_Density'
      myspecies(12) = 'HDI_Density'
      do l=1,nspec
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(abun(i,j,k,l))
         enddo
        enddo
       enddo
       call write_viz(file_id,rank,dims,trim(myspecies(l)),data)
      enddo
      endif ! nspec, abun
!
!
!-----------------------------------------------------------------------
!     Terminate access to the file.
!-----------------------------------------------------------------------
!
      CALL h5fclose_f(file_id, error)
!
!-----------------------------------------------------------------------
!     Close FORTRAN interface.
!-----------------------------------------------------------------------
!
      CALL h5close_f(error)
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine write_viz(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_real,dspace_id, &
                       dset_id,error)
!
      call h5dwrite_f(dset_id, h5t_native_real, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
