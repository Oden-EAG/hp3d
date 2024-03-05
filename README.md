# hp3D
A Scalable MPI/OpenMP hp-Adaptive Finite Element Software Library
for Complex Multiphysics Applications

[![DOI](https://joss.theoj.org/papers/10.21105/joss.05946/status.svg)](https://doi.org/10.21105/joss.05946)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10763374.svg)](https://doi.org/10.5281/zenodo.10763374)

## Downloading the library
1. Clone the repository
- via HTTPS: `git clone https://github.com/Oden-EAG/hp3d.git`
- via SSH: `git clone git@github.com:Oden-EAG/hp3d.git`
2. Access the main directory: `cd hp3d/trunk`

## Compiling the library
1. Create `m_options` file in `hp3d/trunk/`:
Copy one of the existing `m_options` files from `hp3d/trunk/m_options_files/` into `hp3d/trunk/`.
For example: `cp m_options_files/m_options_linux ./m_options`
2. Modify `m_options` file to set the correct path to the main directory:
Set the `HP3D_BASE_PATH` to the path of the `hp3d/trunk/`
3. To compile the library, type `make` in `hp3d/trunk/`. **Before compiling**, you **must** link to the external libraries and set compiler options by modifying the `m_options` file as described below.

- Note: We recommend configuring PETSc with all of hp3D's dependencies (see below). Then, the default `m_options` files `m_options_files/m_options_linux` and `m_options_files/m_options_macos` only require setting the corresponding values of `PETSC_DIR` and `PETSC_ARCH` to link to the required external libraries.

## Linking to external libraries
The `m_options` file must link to the correct paths for external libraries. The following external libraries are used:
- Intel MKL [optional]
- X11 [optional]
- PETSc (all following packages can be installed with PETSc)
- HDF5/pHDF5
- MUMPS
- Metis/ParMetis
- Scotch/PT-Scotch
- PORD
- Zoltan

For example, assuming MPI libraries are installed,
PETSc configure may look like this:
```
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-mpiexec=mpirun \
            --download-fblaslapack=yes \
            --download-scalapack=yes \
            --download-mumps=yes \
            --download-metis=yes \
            --download-parmetis=yes \
            --download-ptscotch=yes \
            --download-zoltan=yes \
            --download-hdf5=yes \
            --with-hdf5-fortran-bindings=1 \
            --with-shared-libraries=0 \
            --with-debugging=0 \
            --with-scalar-type=real \
            --PETSC_ARCH=arch-debian-real
```
Note: PETSc can also install MPI libraries if needed,
e.g. `--download-openmpi=yes` or `--download-mpich=yes`.

## Compiler options
Compilation is governed by preprocessing flags `HP3D_COMPLEX` and `HP3D_DEBUG`.
- `HP3D_COMPLEX = 0` , stiffness matrix, load vector(s) and solution DOFs are real-valued
- `HP3D_COMPLEX = 1` , stiffness matrix, load vector(s) and solution DOFs are complex-valued
- `HP3D_DEBUG   = 0` , compiler uses optimization flags, and the library performs only minimal checks during the computation
- `HP3D_DEBUG   = 1` , compiler uses debug flags, and the library performs additional checks during the computation

Library will be created under either `hp3d/trunk/complex/` or `hp3d/trunk/real/`.

Additional preprocessing flags for enabling/disabling dependencies on third-party libraries:
- `HP3D_USE_INTEL_MKL = 0/1` , disable/enable dependency on Intel MKL package
- `HP3D_USE_MPI_F08   = 0/1` , disable/enable MPI Fortran 2008 binding (module `mpi_f08`)
- `HP3D_USE_OPENMP    = 0/1` , disable/enable OpenMP threading
- `HP3D_USE_X11       = 0/1` , disable/enable dependency on X11

## Verifying build
In addition to the default `make` that builds and installs the hp3D library, the `makefile` provides various targets which can be viewed via `make help`. For example, use `make check` to run a quick check after building the library, or run more extensive tests using `make test`.

## Compiling a problem
Projects are implemented in `hp3d/trunk/problems/`. A few projects have been implemented and can serve as an example. For example, `/problems/POISSON/GALERKIN/` is a Galerkin implementation for the classical variational Poisson problem. To compile and run the problem, type `make`  in the project folder, i.e., `cd problems/POISSON/GALERKIN; make; ./run.sh`.

## Citing hp3D
Please add the following citation to any paper, technical report, or article that incorporated the `hp3D` library:
```bibtex
@article{hp2024joss,
         Author = {Henneking, Stefan and Petrides, Socratis and Fuentes, Federico and Badger, Jacob and Demkowicz, Leszek},
         Title = {{$hp$3D: A Scalable MPI/OpenMP $hp$-Adaptive Finite Element Software Library for Complex Multiphysics Applications}},
         Publisher = {The Open Journal},
         Journal = {Journal of Open Source Software},
         Year = {2024},
         Volume = {9}, 
         Number = {95}, 
         Pages = {5946},
         doi = {10.21105/joss.05946}}
```
And, optionally,
```bibtex
@article{hpUserManual,
         Author = {Henneking, Stefan and Demkowicz, Leszek},
         Title = {{$hp$3D User Manual}},
         Journal={arXiv preprint arXiv:2207.12211},
         Year = {2022}}
```
```bibtex
@book{hpBook3,
      Author = {Henneking, Stefan and Demkowicz, Leszek and Petrides, Socratis and Fuentes, Federico and Keith, Brendan and Gatto, Paolo},
      Title = {{Computing with $hp$ Finite Elements. III. Parallel $hp$3D Code}},
      Publisher = {In preparation},
      Year = {2024}}
```
```bibtex
@book{hpBook2,
      Author = {Demkowicz, Leszek and Kurtz, Jason and Pardo, David and Paszy\'{n}ski, Maciej and Rachowicz, Waldemar and Zdunek, Adam},
      Title = {{Computing with $hp$ Finite Elements. II. Frontiers: Three-Dimensional Elliptic and Maxwell Problems with Applications}},
      Publisher = {Chapman \& Hall/CRC},
      Year = {2007}}
```
```bibtex
@book{hpBook1,
      Author = {Demkowicz, Leszek},
      Title = {{Computing with $hp$ Finite Elements. I. One- and Two-Dimensional Elliptic and Maxwell Problems}},
      Publisher = {Chapman \& Hall/CRC Press, Taylor and Francis},
      Year = {2006}}
```
```bibtex
@article{shape2015,
         Author = {Fuentes, Federico and Keith, Brendan and Demkowicz, Leszek and Nagaraj, Sriram},
         Title = {Orientation embedded high order shape functions for the exact sequence elements of all shapes},
         Publisher = {Elsevier},
         Journal = {Computers \& Mathematics with Applications},
         Year = {2015},
         Volume = {70},
         Number = {4},
         Pages = {353--458}}
```

## User Guide
The user manual is continuously updated and maintained here:
https://github.com/Oden-EAG/hp3d_user_guide (LaTeX source)

A PDF version of the user manual is available on arXiv: https://arxiv.org/abs/2207.12211

## How to Contribute
hp3D is distributed under the terms of the BSD-3 license. All new contributions must be made under this license.

Contributions of all kinds are welcome, including bug fixes, code optimizations, new capabilities, improved documentation, and model problems or applications. The new feature or contribution should be developed on a properly named feature branch based off of `hp3d:master` in a forked repository. If you would like to propose a contribution, please use a pull request (PR) toward the `hp3d:master` branch from your forked hp3D repository. Before starting a PR or working on a new feature, we encourage opening an [issue](https://github.com/Oden-EAG/hp3d/issues) first and discussing the new feature or contribution with the developer team.

## Support
The [issue tracker](https://github.com/Oden-EAG/hp3d/issues) serves as the primary tool for resolving questions related to code features, etc.

For other inquiries, please contact:
``stefan@oden.utexas.edu``, ``leszek@oden.utexas.edu``
