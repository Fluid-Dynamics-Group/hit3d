# HIT3D Conservative Dissipation Control

Control of 3D homogeneous isotropic turbulence with conservative dissipation.

Author: Muralikrishnan Gopalakrishnan Meena, Oak Ridge National Laboratory

The HIT3D code is a modified version of the original code forked from Sergei Chumakov's repository (http://schumakov.info/codes-hit3d.php).

# Hit3DP
HIT3DP is a pseudospectral DNS code, that is, it performs direct numerical
simulation of incompressible isotripic homogeneous turbulence with or without
forcing.  The code has capability of carrying passive scalars, Lagrangian
particles and Large Eddy Simulation


## Extra Packages
The code is written in Fortran 90 and uses the two open libraries:
- Open-MPI  (www.open-mpi.org)
- FFTW3	    (www.fftw.org)

## LICENSING
The code is distributed under the terms of GNU GPLv3 license.  You can read
the full text of the license at http://www.gnu.org/licenses/gpl.html

Copyright (C) 2006-2010 Sergei Chumakov, Natalia Vladimirova, Misha Stepanov


## COMPILING THE CODE
First, edit the Makefile:
- add a section that corresponds to the name of your machine.  Ideally it should
  be a wrapper from your MPI implementation.
- define the name of the F90 compiler
- define FCFLAGS and LDFLAGS.  They should include the include directories, the
  flags that link FFTW3 and MPI implementation.

Run "gmake".

## RUNNING THE CODE
The directory "scripts" provides some examples of the batch job submission files.

The directory "scripts" contains the following files:

00_example.in a sample input file

snapshot.gp	a Gnuplot instruction file that creates two plots that
 		can get attached to the notification emails

coyote.sub Running script for the Coyote cluster at LANL
wcr.sub Example script for WCR cluster at Center for Turbuience Research
		at Stanford University

## MGM's NOTES: COMPLING \& RUNNING THE CODE

Edited: 05/22/2021

1. Edit `MAKEFILE` and add details for your machine. See above and sample machine binaries given in the Makefile.
1. Check input file `sample_inp.in`. File name (also run name) should have 10 characters.
2. Compile:
```
make clean
make
```
3. Run code:
```
mpirun -np <nproc> ./hit3d <run-name> "nosplit"
```
Add this to submission file to run in HPC.

Notes:
1. File extension edited! Check `io_write_4.f90`, `m_io.f90`, `restart_io.f90`
1. Scalar statistics: `stat1.gp` \& `stat2.gp` (see `m_stat.f90`)
1. E-spectra in time: `es.gp` (see `m_stat.f90`)

## THE INPUT FILE
`NX`,`NY`,`NZ`  Number of grid points in one dimension.  The grid will be NX x NY x NZ.
	  The physical dimensions will be 2*pi x 2*pi x 2*pi

`ITMIN`	The timestep number of the restart file.  The restart files have names
	such as "test__this.64.123456".  Here, "test__this" is the run name,
	"64" signifies that the file is written with double precision and
	"123456" is the timestep number.  If the ITMIN is set to 0, the
	subroutine that defines	the initial conditionis for the flow is called.

`ITMAX`	The maximum number of timesteps in the simulation.

`IPRNT1` How often to generate the statistics.

`IPRNT2` How often to write restart files

`IWRITE4` How often to write the real*4 files that are used for post-processing.

`TMAX` The runtime of the simulation (not the wallclocok time)

`TRESCALE` The time at which to rescale the velocity.  This is used in decaying
	  simulations when we want to establish some correlations first and
	  then rescale the velocity field so it has higher kinetic energy.

`TSCALAR` When to start moving the passive scalars.

`flow_type` Parameter that switches the flow type
	  0 - decaying turbulence
	  1 - forced turbulence

`RE` The local Reynolds number (1/nu, where nu is viscosity)

`DT` The timestep.
	  If DT is negative, then the timestep is fixed to be (-DT)
	  If DT is positive, the timestep is found from the stability
	  criteria for the time-stepping scheme that is used.

`ISPCV1` Initial spectrum type (see init_velocity.f90)
`mv1` initial infrared exponent in the spectrum
`wm0v1` initial peak wavenumber in the spectrum


`force_type` The type of the forcing that is applied for the case of forced turbulence.
	* 1 - forcing from Michaels PRL paper (PRL #79(18) p.3411)
	* So far no other forcing has been implemented

`KFMAX`	The upper bound for the forcing band in the Fourier space.

`FAMP` The magnitude of the forcing (usually set to 0.5)

`det_rand` The parameter that switches the random generation for the random seeds for the code.

DEFUNCTIONAL.  In the current version of the code, the seeds for the
random number generator are fixed and are taken from the input file.
The fixed seeds have the effect of producing the initial data that
looks similar for different resolutions (the large features of
initial flow in 32^3 simulation will look similar to the large features
of a 1024^3 simulation if the seeds are the same).

`RN1`, `RN2`, `RN3` - random number seeds


`DEALIAS` The parameter that switches between the dealiasing algorithms.

* 0 - the standard 3/2-rule (or 2/3 rule).  Faster computations, but fewer modes.
* 1 - the phase shift combined with truncation.  This retains much more
	modes than the 2/3-rule, while increasing the computations 60% or so.
	The most economical mode for DNS in terms of flops per the number of
	Fourier modes in the resulting data.

`np` The number of Lagrangian particles

`particle tracking mechanism`:

* 0 - trilinear interpolations
* 1 - 4-point cubic interpolation

`time_p` time in the simulation when to release the particles in the flow

`particle_filter_size`

The particles can be advected by fully resolved field or by locally averaged
field.  The filter size determines the size of the filter that is applied
to the velocity field before computing the particles' velocities.

`les_model` The LES model.  See m_les.f90 for list of the current models.

`NUMS` The number of passive scalars to carry around

The last section contains the parameters of the passive scalars.  Each scalar
must have the type, Schmidt number, infrared exponent, peak wavenumber and
reaction rate.

`TYPE`:
*0	The scalar that is forced by the mean gradient.

* 1-9 The initial conditions for the scalar are generated using Fourier space.
	* 1: Exponential spectrum
	* 2: von Karman spectrum
	* 3: double-delta PDF

* 10 The initial conditions for the scalar are generated in the real space.
	* 11: single slab of the scalar.
	* 12: two slabs of the scalar
	* 13: 1-D sinusoidal wave of the scalar
	* 14: The Aditya gradient

The reaction rate parameter is defunctional in this version of the code.

## Example
The provided example in scripts/00_example.in is somewhat difficult to run out of the box. The
solver requires that the file name (not including the .in extension) is 10 characters in length.
Therefore, I have moved the example file to the project root (./input_file.in).

The solver may be run on the sample input file with

```
./hit3d.x input_file
```

If there is no output some diagnostic information can be found in a newly created ./d00000.txt.

Also note that the solver really doesnt like input file names that start with "NUM_...". I have
no idea why, probably related to how fortran parses command line arguments.

After running the above command a `stat1.gp` and `stat2.gp` file will be produced. These
files can be plotted with `./scripts/snapshot.gp`:

```
./scripts/snapshot.gp
```

which will produce `1.png` and `2.png` for viewing.

## Questions

email Sergei Chumakov at sergei@chumakov.info
