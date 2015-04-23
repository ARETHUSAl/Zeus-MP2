!=======================================================================
!
!                            Developed by
!                Laboratory for Computational Astrophysics
!                  University of California at San Diego
!
      subroutine grav3D_MG
!
!     3-D gravitational potential solver using James Bordner's
!     MGMPI package for Poisson's equation
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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
#ifdef USE_MGMPI
#define SCALAR_DOUBLE
#include "mgmpif.h"
!
!      real(rl) :: mass
      real(rl), dimension(:,:,:), allocatable, save :: &
                ad, a0, a1, a2, x, b
!
      real(rl), dimension(:,:), allocatable :: &
                gpiibc, gpoibc, gpijbc, gpojbc, gpikbc, gpokbc
!
      integer :: i, j, k
!
!     Object handles
!
      integer, save  :: iParameters  ! Handle to parameters object
      integer, save  :: iConcurrent  ! Handle to concurrent object
      integer, save  :: iDistribute  ! Handle to the data distribute object
      integer, save  :: iGrid        ! Handle to the grid object
      integer, save  :: iBc          ! Handle to the boundary conditions 
!                                    ! object
      integer, save  :: iA           ! Handle to the coefficient matrix
      integer, save  :: iSolver      ! Handle to the linear solver
      integer, save  :: iStopping    ! Handle to the stopping criteria
      integer, save  :: iX           ! Handle to the solution vector
      integer, save  :: iB           ! Handle to the right-hand side vector
!
!     Local variables
!
      integer, save  :: i0,i1,i2       ! Gridpoint loop indices
      integer, save  :: ig0,ig1,ig2    ! Global indices
      integer, save  :: gml0,gml1,gml2 ! Global - local index offsets
      integer, save  :: iface          ! Loop index for boundary conditions
      integer, save  :: ng0,ng1,ng2    ! Global problem size
      integer, save  :: nl0,nl1,nl2    ! Local problem size
      integer, save  :: n0,n1,n2       ! Problem size, either local or global
!
      real(rl), parameter :: EXTENT = 1.2D0
!
      real(rl), save :: x0,x1,x2       ! Global position
      real(rl), save :: gxm,gym,gzm    ! Grid lower extents
      real(rl), save :: gxp,gyp,gzp    ! Grid upper extents
      real(rl), save :: h,h0,h1,h2,h012! Mesh widths
      real(rl), save :: rad            ! Distance of position from (0,0,0)
!
      integer , save :: p0,p1,p2       ! Processor grid size
      integer , save :: per0,per1,per2 ! Whether axes are periodic
      integer , save :: ip             ! Processor rank
!
      LOGICAL , save :: lDumpB         ! Whether to dump B to matlab file B.m
      LOGICAL , save :: lDumpX         ! Whether to dump X to matlab file X.m
      LOGICAL , save :: isSizeLocal    ! Whether problem size is 
!                                      ! processor-local
      LOGICAL , save :: lAlias         ! Whether to use VectorAlias or VectorCopy
!
      CHARACTER, save :: bcType(MGMPI_FACE_FIRST:MGMPI_FACE_LAST)*10
!                                  Boundary condition types
!JB      INTEGER comm              ! Handle to MPI communicator
!
!
!-----------------------------------------------------------------------
!     First timestep only: initialize everything required for MGMPI
!-----------------------------------------------------------------------
!
      if(nhy .eq. 0) then
!     ------------------
!     Initialize MGMPI
!     ------------------
       CALL MGMPI_Init
!     ---------------------------------------
!     Open parametersuration file and read values
!     ---------------------------------------
       CALL MGMPI_ParametersCreate (iParameters)
       CALL MGMPI_ParametersReadFile (iParameters,'mgmpi_inp')
!     ---------------------------------------
!     Create Concurrenct object
!     ---------------------------------------
       CALL MGMPI_ConcurrentCreate (iConcurrent, iParameters, &
                                    comm3d)
!
!JB      comm = MPI_COMM_NULL
!JB      CALL MGMPI_ConcurrentCreate (iConcurrent, iParameters, comm)
!JB      ip = MGMPI_ConcurrentRank (iConcurrent)
!
       ip = myid
!
!     Get processor grid size
!
       p0 = ntiles(1)
       p1 = ntiles(2)
       p2 = ntiles(3)
!
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:P0', p0, 1)
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:P1', p1, 1)
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:P2', p2, 1)
!
!     Determine if processor grid should be periodic
!
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcXM', bcType(MGMPI_FACE_LOWER_X), 'Undefined')
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcYM', bcType(MGMPI_FACE_LOWER_Y), 'Undefined')
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcZM', bcType(MGMPI_FACE_LOWER_Z), 'Undefined')
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcXP', bcType(MGMPI_FACE_UPPER_X), 'Undefined')
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcYP', bcType(MGMPI_FACE_UPPER_Y), 'Undefined')
       CALL MGMPI_ParametersGetValue (iParameters,
     $     'Test:BcZP', bcType(MGMPI_FACE_UPPER_Z), 'Undefined')
       per0 = 0
       per1 = 0
       per2 = 0
       IF ( bcType(MGMPI_FACE_LOWER_X).EQ.'Periodic' .OR.
     $     bcType(MGMPI_FACE_UPPER_X).EQ.'Periodic') per0 = 1
       IF ( bcType(MGMPI_FACE_LOWER_Y).EQ.'Periodic' .OR.
     $     bcType(MGMPI_FACE_UPPER_Y).EQ.'Periodic') per1 = 1
       IF ( bcType(MGMPI_FACE_LOWER_Z).EQ.'Periodic' .OR.
     $     bcType(MGMPI_FACE_UPPER_Z).EQ.'Periodic') per2 = 1
       CALL MGMPI_ConcurrentSetAxisSizes (iConcurrent, p0,p1,p2,
     $     per0, per1, per2)
!     -------------------------------------------------------
!     Create and initialize a 3D Cartesian processor topology
!     -------------------------------------------------------
       CALL MGMPI_DistributeCreate (iDistribute,iParameters,iConcurrent)
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:N0',n0, 0)
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:N1',n1, 0)
!JB      CALL MGMPI_ParametersGetInteger (iParameters,'Test:N2',n2, 0)
       n0 = ie-is+1
       n1 = je-js+1
       n2 = ke-ks+1
       CALL MGMPI_ParametersGetLogical (iParameters,'Test:SizeIsLocal',
     $     isSizeLocal,.FALSE.)
       IF (isSizeLocal) THEN
          nl0 = n0
          nl1 = n1
          nl2 = n2
          CALL MGMPI_DistributeSetLocalSize (iDistribute,nl0,nl1,nl2)
          CALL MGMPI_DistributeGetGlobalSize (iDistribute,ng0,ng1,ng2)
       ELSE
          ng0 = n0
          ng1 = n1
          ng2 = n2
          CALL MGMPI_DistributeSetGlobalSize (iDistribute,ng0,ng1,ng2)
          CALL MGMPI_DistributeGetLocalSize (iDistribute,nl0,nl1,nl2)
       END IF
!     -------------------------------------------------------
!     Define boundary conditions
!     -------------------------------------------------------
       CALL MGMPI_BcCreate (iBc, iParameters, iDistribute)
!     (bcType(:) already read above)
!     Set face types and values
!     (though values not needed for periodic axes)
       DO iface=MGMPI_FACE_FIRST, MGMPI_FACE_LAST
          IF (bcType(iface).EQ.'Dirichlet') THEN
             CALL MGMPI_BcSetFaceType (iBc,iface,MGMPI_BC_DIRICHLET)
          ELSE IF (bcType(iface).EQ.'Neumann') THEN
             CALL MGMPI_BcSetFaceType (iBc,iface,MGMPI_BC_NEUMANN)
          ELSE IF (bcType(iface).EQ.'Periodic') THEN
             CALL MGMPI_BcSetFaceType (iBc,iface,MGMPI_BC_PERIODIC)
          ENDIF
!JH          CALL MGMPI_BcSetFaceValue (iBc, iface, 0.0D0)
       END DO
!     -------------------------------------------------------
!     Create a distributed grid using the given distribution
!     -------------------------------------------------------
       CALL MGMPI_GridCreate (iGrid, iParameters, iDistribute)
!
!     (NOTE: GRID EXTENTS NOT ACCESSED SINCE Pde/Discretize NOT USED)
!
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridXM",gxm, 0.0D0)
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridYM",gym, 0.0D0)
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridZM",gzm, 0.0D0)
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridXP",gxp, 0.0D0)
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridYP",gyp, 0.0D0)
      CALL MGMPI_ParametersGetReal(iParameters,"Test:GridZP",gzp, 0.0D0)
!
!      gxm = -EXTENT ! x1min
!      gym = -EXTENT ! x2min
!      gzm = -EXTENT ! x3min
!      gxp =  EXTENT ! x1max
!      gyp =  EXTENT ! x2max
!      gzp =  EXTENT ! x3max
      CALL MGMPI_GridSetExtents (iGrid, gxm, gym, gzm, gxp, gyp, gzp)
!     -------------------------------------------------------
!     Create matrix, right-hand side, and solution vectors on the grid
!     -------------------------------------------------------
       CALL MGMPI_MatrixCreate (iA, iParameters, iGrid, iGrid)
       CALL MGMPI_VectorCreate (iX, iParameters, iGrid)
       CALL MGMPI_VectorCreate (iB, iParameters, iGrid)
!
       CALL MGMPI_MatrixSetBc (iA, iBc, iBc)
!
!-----------------------------------------------------------------------
!      allocate matrix, solution, and RHS array storage
!-----------------------------------------------------------------------
!
       allocate(ad(IN,JN,KN))
       allocate(a0(IN,JN,KN))
       allocate(a1(IN,JN,KN))
       allocate(a2(IN,JN,KN))
       allocate(x (0:IN-1,0:JN-1,0:KN-1))
       allocate(b (0:IN-1,0:JN-1,0:KN-1))
      endif ! nhy = 0
!
!-----------------------------------------------------------------------
!     compute BC arrays
!-----------------------------------------------------------------------
!
      call gpbv
!
      allocate(gpiibc(nl1,nl2))
      allocate(gpoibc(nl1,nl2))
      allocate(gpijbc(nl2,nl0))
      allocate(gpojbc(nl2,nl0))
      allocate(gpikbc(nl0,nl1))
      allocate(gpokbc(nl0,nl1))
!
      do k = ks, ke
       do j = js, je
        gpiibc(j-js+1,k-ks+1) = -guniv*gpiib(j,k,1)
        gpoibc(j-js+1,k-ks+1) = -guniv*gpoib(j,k,1)
       enddo
      enddo
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_LOWER_X, gpiibc,nl1, &
                                                                 nl2)
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_UPPER_X, gpoibc,nl1, &
                                                                 nl2)
!
!     note the difference between mgmpi and zeusmp conventions on
!     index order for gp[i/o]jbc!
!
      do i = is, ie
       do k = ks, ke
        if(lgeom .eq. 2) then
         gpijbc(k-ks+1,i-is+1) = 0.0D0
        else
         gpijbc(k-ks+1,i-is+1) = -guniv*gpijb(i,k,1)
        endif
        gpojbc(k-ks+1,i-is+1) = -guniv*gpojb(i,k,1)
       enddo
      enddo
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_LOWER_Y, gpijbc,nl2, &
                                                                 nl0)
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_UPPER_Y, gpojbc,nl2, &
                                                                 nl0)
      do j = js, je
       do i = is, ie
        gpikbc(i-is+1,j-js+1) = -guniv*gpikb(i,j,1)
        gpokbc(i-is+1,j-js+1) = -guniv*gpokb(i,j,1)
       enddo
      enddo
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_LOWER_Z, gpikbc,nl0, &
                                                                 nl1)
      CALL MGMPI_BcSetFaceValues(iBc, MGMPI_FACE_UPPER_Z, gpokbc,nl0, &
                                                                 nl1)
!
!-----------------------------------------------------------------------
!     Load arrays for matrix
!-----------------------------------------------------------------------
!
      do i2 = 1,nl2
       do i1 = 1,nl1
        do i0 = 1,nl0
         i = i0+2
         j = i1+2
         k = i2+2
         a0(i0,i1,i2) = dvl2a(j)*dvl3a(k)*g2a(i+1)*g31a(i+1) &
                      * dx1bi(i+1)
         a1(i0,i1,i2) = dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(j+1) &
                      * dx2bi(j+1)
         a2(i0,i1,i2) = dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2 &
                      * dx3bi(k+1)
         ad(i0,i1,i2) = &
                  -dvl2a(j)*dvl3a(k)*(g2a(i+1)*g31a(i+1)*dx1bi(i+1)+ &
                                      g2a(i  )*g31a(i  )*dx1bi(i  )) &
                  -dvl1a(i)*dvl3a(k)*g2bi(i)**2* &
                             (g32a(j)*dx2bi(j)+g32a(j+1)*dx2bi(j+1)) &
                  -dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2* &
                                             (dx3bi(k) + dx3bi(k+1))
         x(i0,i1,i2)  = 0.0D0
        enddo
       enddo
      enddo
!
      if(nhy .eq. 0) then
       do i2 = 1, kn
        do i1 = 1, jn
         do i0 = 1, in
           gp(i0,i1,i2) = 0.0D0
         enddo
        enddo
       enddo
      else
       do i2 = 1, kn
        do i1 = 1, jn
         do i0 = 1, in
           gp(i0,i1,i2) = -gp(i0,i1,i2)
         enddo
        enddo
       enddo
      endif ! nhy
!     Initialize right-hand side array b[]
      DO i2 = 1,nl2
         DO i1 = 1,nl1
            DO i0 = 1,nl0
               b(i0,i1,i2) = 0.0D0
            END DO
         END DO
      END DO
      if(nhy .eq. 0) CALL MGMPI_DistributeGetOffsets &
                                          (iDistribute,gml0,gml1,gml2)
!     Mesh widths (correct for Dirichlet; else replace (ngn+1) with ngn)
      h = (h0 + h1 + h2) / 3.0D0
      h012 = h0*h1*h2
      DO i2=1,nl2
         ig2 = i2 + gml2 - 1
!JH         x2 = -EXTENT + h2*ig2
         x2 = gzm + h2*ig2
         DO i1=1,nl1
            ig1 = i1 + gml1 - 1
!JH            x1 = -EXTENT + h1*ig1
            x1 = gym + h1*ig1
            DO i0=1,nl0
               ig0 = i0 + gml0 - 1
!JH               x0 = -EXTENT + h0*ig0
               x0 = gxm + h0*ig0
               rad = sqrt(x0*x0 + x1*x1 + x2*x2)
               b(i0,i1,i2) = 4.0D0*d(i0+is-1,i1+js-1,i2+ks-1)* &
                             PI*guniv*dvl1a(i0+is-1)*dvl2a(i1+js-1)* &
                                      dvl3a(i2+ks-1)
            END DO
         END DO
      END DO
!     Copy matrix arrays to Matrix class iA
      CALL MGMPI_MatrixStoreDiagonal(iA, 0, 0, 0,ad,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA, 1, 0, 0,a0,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA, 0, 1, 0,a1,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA, 0, 0, 1,a2,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA,-1, 0, 0,a0,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA, 0,-1, 0,a1,IN,JN,KN)
      CALL MGMPI_MatrixStoreDiagonal(iA, 0, 0,-1,a2,IN,JN,KN)
      if(nhy .eq. 0) then
       CALL MGMPI_ParametersGetLogical (iParameters, 'Test:Alias',
     $      lAlias, .FALSE.)
      endif
      IF (lAlias) THEN
         CALL MGMPI_VectorAliasVertices (iX,gp(2,2,2),in,jn,kn)
         CALL MGMPI_VectorAliasVertices (iB,b(0,0,0),in,jn,kn)
      ELSE
         CALL MGMPI_VectorStoreUnknowns (iX,x(1,1,1),in,jn,kn)
         CALL MGMPI_VectorStoreUnknowns (iB,b(1,1,1),in,jn,kn)
      ENDIF
!     Dump B to Matlab file B.m if requested
!      CALL MGMPI_ParametersGetLogical (iParameters,
!     $     'Test:DumpRightHandSide', lDumpB, .FALSE.)
!      IF (lDumpB)  CALL MGMPI_VectorWrite(iB,'B')
!     -------------------------------------------------------
!     Solve the linear system
!     -------------------------------------------------------
      if(nhy .eq. 0) then
       CALL MGMPI_StoppingCreate (iStopping, iParameters)
       CALL MGMPI_SolverCreate (iSolver, iParameters, iStopping)
      endif
      CALL MGMPI_SolverApply (iSolver, iA, iX, iB)
!
!-----------------------------------------------------------------------
!     in theory, the solution is in hand at this point
!-----------------------------------------------------------------------
!
!     -------------------------------------------------------
!     Finish up
!     -------------------------------------------------------
!
!     positive potential, consistent with ZEUS convention
!
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         gp(i,j,k) = -gp(i,j,k)
        enddo
       enddo
      enddo
!
!     scale BC arrays to CGS
!
      do k = ks, ke
       do j = js, je
        gpiib(j,k,1) = guniv*gpiib(j,k,1)
        gpiib(j,k,2) =       gpiib(j,k,1)
        gpoib(j,k,1) = guniv*gpoib(j,k,1)
        gpoib(j,k,2) =       gpoib(j,k,1)
       enddo
      enddo
      do k = ks, ke
       do i = is, ie
        gpijb(i,k,1) = guniv*gpijb(i,k,1)
        gpijb(i,k,2) =       gpijb(i,k,1)
        gpojb(i,k,1) = guniv*gpojb(i,k,1)
        gpojb(i,k,2) =       gpojb(i,k,1)
       enddo
      enddo
      do j = js, je
       do i = is, ie
        gpikb(i,j,1) = guniv*gpikb(i,j,1)
        gpikb(i,j,2) =       gpikb(i,j,1)
        gpokb(i,j,1) = guniv*gpokb(i,j,1)
        gpokb(i,j,2) =       gpokb(i,j,1)
       enddo
      enddo
!
      deallocate(gpiibc)
      deallocate(gpoibc)
      deallocate(gpijbc)
      deallocate(gpojbc)
      deallocate(gpikbc)
      deallocate(gpokbc)
!     Dump X to Matlab file X.m if requested
      CALL MGMPI_ParametersGetLogical (iParameters,'Test:DumpSolution',
     $     lDumpX, .FALSE.)
      IF (lDumpX) CALL MGMPI_VectorWrite(iX,'X')
      if(nhy .eq. nlim-1) then
!JB       CALL MGMPI_ParametersWriteFile (iParameters,
!JB     $      'out.ftest.parameters')
       CALL MGMPI_ParametersDelete (iParameters)
!       CALL MGMPI_TimerDelete (iTimer)
!     @@@ BUG F2! : ~Hierarchy3D() dies when deallocating smoother @@@
!      CALL MGMPI_ConcurrentDelete (iConcurrent)
       CALL MGMPI_DistributeDelete (iDistribute)
       CALL MGMPI_GridDelete (iGrid)
       CALL MGMPI_BcDelete (iBc)
       CALL MGMPI_VectorDelete (iB)
       CALL MGMPI_VectorDelete (iX)
       CALL MGMPI_MatrixDelete (iA)
       CALL MGMPI_StoppingDelete (iStopping)
       CALL MGMPI_SolverDelete (iSolver)
      endif ! nhy
!
#endif /* USE_MGMPI */
      return
      end
