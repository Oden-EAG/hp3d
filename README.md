# hp3d
A Scalable MPI/OpenMP hp-Adaptive Finite Element Software Library
for Complex Multiphysics Applications

## Downloading the library
1. Clone the repository
- via HTTPS: `git clone https://github.com/Oden-EAG/hp3d.git`
- via SSH: `git clone git@github.com:Oden-EAG/hp3d.git`
2. Access the main directory: `cd hp3d/trunk`

## Compiling the library
1. Create `m_options` file in `hp3d/trunk/`:
Use the default `m_options` file in `hp3d/trunk/` or copy one of the existing `m_options` files from `hp3d/trunk/m_options_files/` into `hp3d/trunk/`.
For example: `cp m_options_files/m_options_TACC_intel18 m_options`
2. Modify `m_options` file to set the correct path to the main directory:
Set the `HP3D_BASE_PATH` to the path of the `hp3d/trunk/`
3. To compile the library, type `make` in `hp3d/trunk/`. **Before compiling**, you **must** link to the external libraries and set compiler options by modifying the `m_options` file as described below.

## Linking to external libraries
The `m_options` file must link to the correct paths for external libraries. The following external libraries are used:
- Intel MKL [optional]
- X11
- PETSc (all following packages can be installed with PETSc)
- HDF5/pHDF5
- MUMPS
- Metis/ParMetis
- Scotch/PT-Scotch
- PORD
- Zoltan

## Compiler options
Compilation is governed by preprocessing flags `COMPLEX` and `DEBUG`.
- `COMPLEX = 0` , stiffness matrix, load vector(s) and solution DOFs are real-valued
- `COMPLEX = 1` , stiffness matrix, load vector(s) and solution DOFs are complex-valued
- `DEBUG   = 0` , compiler uses optimization flags and the library performs only minimal checks during the computation
- `DEBUG   = 1` , compiler uses debug flags, and the library performs additional checks during the computation

Library will be created under either `hp3d/complex/` or `hp3d/real/`.

Additional preprocessing flags for enabling/disabling third-party libraries:
- `HP3D_USE_INTEL_MKL = 0` , dependency on Intel MKL package is disabled
- `HP3D_USE_INTEL_MKL = 1` , additional solver options are available to the user (e.g., Intel MKL PARDISO)

## Compiling a problem
Projects are implemented in `hp3d/trunk/problems/`. A few projects have been implemented and can serve as an example. For example, `/problems/MPI_POISSON/GALERKIN/` is a Galerkin implementation for the classical variational Poisson problem. To compile and run the problem, type `make`  in the project folder, i.e., `cd problems/MPI_POISSON/GALERKIN; make; ./run.sh`.

## Citing hp3d
Please add the following citation to any paper, technical report, or article that incorporated the `hp3d` library:
```bibtex
@book{hpbook2,
      Author = {Demkowicz, L. and Kurtz, J. and Pardo, D. and Paszy\'{n}ski, M. and Rachowicz, W. and Zdunek, A.},
      Publisher = {Chapman \& Hall/CRC},
      Title = {Computing with $hp$ Finite Elements. II. Frontiers: Three-Dimensional Elliptic and Maxwell Problems with Applications},
      Year = {2007}}
```
And, optionally,
```bibtex
@book{hpbook,
      Author = {Demkowicz, L.},
      Publisher = {Chapman \& Hall/CRC Press, Taylor and Francis},
      Title = {Computing with $hp$ Finite Elements. I. One- and Two-Dimensional Elliptic and Maxwell Problems},
      Year = {2006}}
```
```bibtex
@article{fuentes2015shape,
    Author = {Fuentes, F. and Keith, B. and Demkowicz, L. and Nagaraj, S.},
    Journal = {Comput. Math. Appl.},
    Number = {4},
    Pages = {353--458},
    Publisher = {Elsevier},
    Title = {Orientation embedded high order shape functions for the exact sequence elements of all shapes},
    Volume = {70},
    Year = {2015}}
```

## User Guide
... in development

## Support
Contact: ``stefan@oden.utexas.edu``, ``leszek@oden.utexas.edu``
