This is the development branch of ZEUS-MP/2.  
--------------------------------------------

It is based off the last stable release by John Hayes and [downloaded
here](http://lca.ucsd.edu/portal/software/zeus-mp2). (The latest stable can
also be downloaded from the Downloads tab.)

ZEUS-MP/2 is a computational fluid dynamics code for the simulation of
astrophysical phenomena, based on ZEUS-3D and parallelized using the MPI
message-passing library.

Project Description:
--------------------

ZEUS-MP is the latest addition to the ZEUS line of community application
codes developed by the Laboratory for Computational Astrophysics. The
"MP" suffix denotes the "multi-physics," "massively parallel," and
"message passing" aspects of the code. The physics suite in this release
of ZEUS-MP includes gas hydrodynamics, ideal MHD, flux-limited radiation
diffusion, self gravity, and multispecies advection. Hydrodynamic,
radiation-hydrodynamic (RHD), and magnetohydrodynamic (MHD) simulations
can be performed on 1, 2, or 3-dimensional grids. Self-gravity is
supported in all dimensions and can be calculated with a variety of
methods, depending on geometry, dimensionality, and boundary conditions.
Methods include: 

1.  GM/r on 1D or 2D spherical grids; 

2. a conjugate-gradient Poisson solver on 2D or 3D spherical and cylindrical
   grids; 

3. a multigrid-based Poisson solver on 3D cartesian grids with
   non-periodic boundary conditions, and 

4. a Fast Fourier Transform
   solver for 3D cartesian grids with triply-periodic boundaries. External
   point-mass potentials are also supported.

Contributing:
-------------

All bug reports and patches are welcome. Please make sure all your changes
successfully merge into the latest tip before you make a pull request.
# Zeus-MP2
