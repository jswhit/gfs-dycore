simplfied version of the NOAA/NCEP global forecast system ([GFS](http://www.emc.ncep.noaa.gov/index.php?branch=GFS)) weather model..

Identical to the operational GFS except that...

1) No MPI.  Runs on a single node with openmp, and uses the blazing fast [shtns](http://users.isterre.fr/nschaeff/SHTns) spherical harmonic library.  This is partly to make the code simpler to read, but also with the idea that when you're running large ensembles you can use MPI to parallelize along the ensemble dimension while running each ensemble member on a single node with shared memory parallelism.  In the GPU era, this may turn out to be the most efficient way to run ensembles.

2) Uses a forward in time integration scheme instead of leapfrog. This makes restarts much easier (only a single snapshot in time needed, and no time filtering since there is no computational mode).  I used a 2nd-order semi-implicit Runge-Kutta scheme (http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-13-00132.1).  You can actually take a time step that is about 1.5 times as long as leapfrog/semi-implicit, but you have to evaluate the tendencies three times per time step.  The physics tendencies are evaluated only once per dynamics time step.

3) not all of the available physics options are implemented (although the current operational configuration is).

4) uses a full Gaussian grid, not a reduced grid (using a reduced grid will speed things up by about 20-25% in the physics, but will not speed up the dynamics significantly).

Since the code is easier to read and modify than the operational version, it can be used a testbed for new parameterization schemes (I'm particularly interested in stochastic parameterization).  It also can serve as a vehicle for exploring the use of GPUs and Intel MIC architectures for accelerating ensemble forecasts.

If you're curious, go ahead and [browse](https://code.google.com/p/gfs-dycore/source/browse/#svn%2Ftrunk%2Fsrc) the source code.