!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine msave(filename)
!
!  PURPOSE:  Writes [reads] all common block variables to the logical
!  unit 4.  Currently, the common blocks written [read] are:
!            common /gridvarr/ = real    grid     variables
!            common /gridvari/ = integer grid     variables
!            common /fieldr  / = real    field    variables
!            common /fieldi  / = integer field    variables
!            common /bndryr  / = real    boundary variables
!            common /bndryi  / = integer boundary variables
!            common /rootr   / = real    root     variables
!            common /rooti   / = integer root     variables
!            common /gravcomr/ = real    gravity  variables
!            common /gravcomi/ = integer gravity  variables
!
!  The following blocks are NOT written [read]:
!            common /frmlcomr/ = real    formal soln variables
!            common /frmlcomi/ = integer formal soln variables
!            common /mmntcomr/ = real    moment soln variables
!            common /mmntcomi/ = integer moment soln variables
! 
!  EXTERNALS: [none]
!
!  LOCALS:
!
!  MODIFIED: 26 Aug. 1996 by RAF for ZEUS-MP.
!  MODIFIED:  6 Jan. 1999 by JCH
!  MODIFIED:  7 Jan. 1999 by JCH (increased nrootr for root.h;
!                                 makes room for rad t-step var.)
!  MODIFIED: 18 May. 1999 by efh (increased nrootr due to tslice)
!  MODIFIED: December, 2005 by JCH; increased nfieldr for abun array
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use field
      use bndry
      use root
      use gravmod
      use restart_arrays
!
#ifdef USE_HDF5
      use hdf5
#endif
!
      implicit NONE
!
      character*16 :: filename
      character*56 :: filename2
      character*6 :: incrf
#ifdef USE_HDF5
!
!-----------------------------------------------------------------------
!     hdf5-specific parameters, identifiers, etc.
!-----------------------------------------------------------------------
!
      integer        :: rank, error
      integer(hid_t) :: file_id
      integer(hsize_t), dimension(7) :: dims
#endif /* USE_HDF5 */
!
!=======================================================================
!
!JMS  Write restart dumps out to a directory to match Enzo for YT.    
      incrf = filename(13:)
      CALL system('mkdir -p '//'RD'//trim(incrf))
      filename2 = 'RD'//trim(incrf)//'/'//trim(filename)
#ifndef USE_HDF5
      open(unit=4,file=trim(filename2),status='unknown', &
           form='unformatted')
#endif
!=======================================================================
!
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
!
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
!
!      nbdryr = (2*6*nspec) * (jn*kn + in*jn + in*kn) + 6*nbvar
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar &
                        + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
!
      call mapout
!
#ifndef USE_HDF5
      write(4)  rlgrdvr , ntgrdvr &
               ,rlfldvr &
               ,rlbdryvr , ntbdryvr &
               ,rlrtvr   , ntrtvr , chrtvr &
               ,rlgrvvr  , ntgrvvr
      if (myid .eq. 0) &
        write(2,"(/1x,'restart dump written at time=',1pe12.5,' cycle=' &
       ,i6)") time,nhy
      close(unit=4)
#else /* USE_HDF5 */
!
!-----------------------------------------------------------------------
!     Initialize FORTRAN interface.
!-----------------------------------------------------------------------
!
      CALL h5open_f (error)
!
!-----------------------------------------------------------------------
!     Create a new file using default properties.
!-----------------------------------------------------------------------
!
      CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id, error)
!
!-----------------------------------------------------------------------
!     Create/Write datasets
!-----------------------------------------------------------------------
!
      RANK = 1
!
      dims(1  ) = ngridr+mgridr
      dims(2:7) = 0
!
      call write_real_res(file_id,rank,dims,"grid_real",rlgrdvr)
!
      dims(1  ) = ngridi
      dims(2:7) = 0
!
      call write_int_res(file_id,rank,dims,"grid_int",ntgrdvr)
!
      dims(1  ) = nfieldr
      dims(2:7) = 0
!
      call write_real_res(file_id,rank,dims,"field_real",rlfldvr)
!
      dims(1  ) = nbdryr
      dims(2:7) = 0
!
      call write_real_res(file_id,rank,dims,"bndry_real",rlbdryvr)
!
      dims(1  ) = nbdryi
      dims(2:7) = 0
!
      call write_int_res(file_id,rank,dims,"bndry_int",ntbdryvr)
!
      dims(1  ) = nrootr
      dims(2:7) = 0
!
      call write_real_res(file_id,rank,dims,"root_real",rlrtvr)
!
      dims(1  ) = nrooti
      dims(2:7) = 0
!
      call write_int_res(file_id,rank,dims,"root_int",ntrtvr)
!
      dims(1  ) = nrootch
      dims(2:7) = 0
!
      call write_chr_res(file_id,rank,dims,"root_chr",chrtvr)
!
      dims(1  ) = ngravr
      dims(2:7) = 0
!
      call write_real_res(file_id,rank,dims,"grav_real",rlgrvvr)
!
      dims(1  ) = ngravi
      dims(2:7) = 0
!
      call write_int_res(file_id,rank,dims,"grav_int",ntgrvvr)
#endif /* USE_HDF5 */
!
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
!
      return
!
!-----------------------------  MGET  ----------------------------------
!
      entry mget(filename)
      incrf = filename(13:)
      filename2 = 'RD'//trim(incrf)//'/'//trim(filename)
#ifndef USE_HDF5
      open(unit=4,file=trim(filename2),status='old',form='unformatted')
#endif
!
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
!
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
!
!      nbdryr = (2*6*nspec) * (jn*kn + in*jn + in*kn) + 6*nbvar
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar &
                        + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
!
#ifndef USE_HDF5
      read(4)   rlgrdvr , ntgrdvr &
               ,rlfldvr &
               ,rlbdryvr , ntbdryvr &
               ,rlrtvr , ntrtvr , chrtvr &
               ,rlgrvvr , ntgrvvr
      close(unit=4)
#else /* USE_HDF5 */
!
!-----------------------------------------------------------------------
!     Initialize FORTRAN interface.
!-----------------------------------------------------------------------
!
      call h5open_f (error)  ! initialize f90 interface
!
!-----------------------------------------------------------------------
!     Open restart file
!-----------------------------------------------------------------------
!
      call h5fopen_f(trim(filename2),h5f_acc_rdonly_f,file_id,error)
!
!-----------------------------------------------------------------------
!     Read datasets
!-----------------------------------------------------------------------
!
      RANK = 1
!
      dims(1  ) = ngridr+mgridr
      dims(2:7) = 0
!
      call read_real_res(file_id,rank,dims,"grid_real",rlgrdvr)
!
      dims(1  ) = ngridi
      dims(2:7) = 0
!
      call read_int_res(file_id,rank,dims,"grid_int",ntgrdvr)
!
      dims(1  ) = nfieldr
      dims(2:7) = 0
!
      call read_real_res(file_id,rank,dims,"field_real",rlfldvr)
!
      dims(1  ) = nbdryr
      dims(2:7) = 0
!
      call read_real_res(file_id,rank,dims,"bndry_real",rlbdryvr)
!
      dims(1  ) = nbdryi
      dims(2:7) = 0
!
      call read_int_res(file_id,rank,dims,"bndry_int",ntbdryvr)
      dims(1  ) = nrootr
      dims(2:7) = 0
!
      call read_real_res(file_id,rank,dims,"root_real",rlrtvr)
!
      dims(1  ) = nrooti
      dims(2:7) = 0
!
      call read_int_res(file_id,rank,dims,"root_int",ntrtvr)
!
      dims(1  ) = nrootch
      dims(2:7) = 0
!     
      call read_chr_res(file_id,rank,dims,"root_chr",chrtvr)
!
      dims(1  ) = ngravr
      dims(2:7) = 0
!
      call read_real_res(file_id,rank,dims,"grav_real",rlgrvvr)
!
      dims(1  ) = ngravi
      dims(2:7) = 0
!
      call read_int_res(file_id,rank,dims,"grav_int",ntgrvvr)
#endif /* USE_HDF5 */
!
      call mapin
!
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
!
      return
      end
#ifdef USE_HDF5
!=======================================================================
!=======================================================================
      subroutine write_real_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real(kind=8) :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_double,dspace_id, &
                       dset_id,error)
!
      call h5dwrite_f(dset_id, h5t_native_double, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine write_int_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      integer :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_integer,dspace_id, &
                       dset_id,error)
!
      call h5dwrite_f(dset_id, h5t_native_integer, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine read_real_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real(kind=8) :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
!
      call h5dread_f(dset_id, h5t_native_double, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine read_int_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      integer :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
!
      call h5dread_f(dset_id, h5t_native_integer, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine read_chr_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      character :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
!
      call h5dread_f(dset_id, h5t_native_character, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
!=======================================================================
!=======================================================================
      subroutine write_chr_res(file_id,rank,dims,dsetname,dset)
!
      use hdf5
!
      implicit none
!
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      character :: dset
!
!
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
!
      call h5screate_simple_f(rank, dims, dspace_id, error)
!
!                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_character,dspace_id, &
                       dset_id,error)
!
      call h5dwrite_f(dset_id, h5t_native_character, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
!
      return
      end
#endif /* USE_HDF5 */
