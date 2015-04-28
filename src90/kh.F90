!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                K H                        \\\\\\\\\\
!
!=======================================================================
!
subroutine kh
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
      integer  :: i,j,k
      real(rl) :: v0, rho_in, rho_out
      namelist /pgen/ v0, rho_in, rho_out

      v0      = 0.5
      rho_in = 10.0
      rho_out = 1.0
      if(myid.eq.0)then
         write(*,*)"K-H problem setup."
      endif
!
!     read in namelist and write it to log file
!
      if (myid .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
        buf_in( 1)  = v0
        buf_in( 2)  = rho_in
        buf_in( 3)  = rho_out
      endif
      call MPI_BCAST(ibuf_in, 3, MPI_FLOAT , 0, comm3d, ierr )
      if (myid .ne. 0) then
        v0      = buf_in(1)
        rho_in  = buf_in(2)
        rho_out = buf_in(3)
#endif /* MPI_USED */
      endif
!
!
      do k=ks-2,ke+2
        do j=js-2,je+2
          do i=is-2,ie+2
            e(i,j,k)=2.5
            v2(i,j,k)=0.0
            v3(i,j,k)=0.0
            if(x2a(j).le.-0.25)then
              d(i,j,k)=rho_out
              v1(i,j,k)=-v0
            endif
            if(x2a(j).ge.0.25)then
              d(i,j,k)=rho_out
              v1(i,j,k)=-v0 
            endif
            if((x2a(j).gt.-0.25).and.(x2a(j).lt.0.25))then
              d(i,j,k)=rho_in
              v1(i,j,k)=v0 
            endif
          enddo !i
        enddo !j
      enddo !k
!
      return
end subroutine kh
