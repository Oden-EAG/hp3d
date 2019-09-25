# par_hp3d
A Scalable MPI/OpenMP hp-Adaptive Finite Element Software Library
for Complex Multiphysics Applications

## Downloading the library
1. Clone the repository
via HTTPS: `git clone https://github.com/ICES-EAG/par_hp3d.git`
via SSH: `git clone git@github.com:ICES-EAG/par_hp3d.git`
2. Access the main directory
`cd par_hp3d/trunk`

## Compiling the library
1. Create `m_options` file in `par_hp3d/trunk/`
Use the default `m_options` file in `par_hp3d/trunk/` or copy one of the existing `m_options` files from `par_hp3d/trunk/m_options_files/` into `par_hp3d/trunk/`. For example: `cp m_options_files/m_options_TACC_intel18 m_options`
2. Modify `m_options` file to set the correct path to the main directory
Set the `HP3D_BASE_PATH` to the path of the `par_hp3d/trunk/`
3. To compile the library, type `make` in `par_hp3d/trunk/`. **Before compiling**, you **must** link to the external libraries and set compiler options by modifying the `m_options` file as described below.

## Linking to external libraries
The `m_options` file must link to the correct paths for external libraries. The following external libraries are requried:
- Intel MKL
- MUMPS
- Metis/ParMetis
- Scotch/PT-Scotch
- PORD
- Zoltan

## Compiler options
Compilation is governed by preprocessing flags `COMPLEX`, `SHAPE`, and `DEBUG`.
- `COMPLEX  = 0` , stiffness matrix, load vector(s) and solution dofs are real-valued
- `COMPLEX  = 1` , stiffness matrix, load vector(s) and solution dofs are complex-valued
- `SHAPE = 0` , shape functions of Gatto, Demkowicz, Kim, 2008 (H1,H(curl) for selected element shapes only)
- `SHAPE = 1` , shape functions of Fuentes et al., 2014 (all spaces, all shapes)
- `DEBUG = YES` , compiler uses debug flags, and the library performs additional checks during the computation
- `DEBUG = NO` , compiler uses optimization flags and the library performs only minimal checks during the computation

Library will be created under either `par_hp3d/complex/` or `par_hp3d/real/`.

## Compiling a problem
Projects are implemented in `par_hp3d/trunk/problems/`. A few projects have been implemented and can serve as an example. For example, `/problems/MPI_POISSON/GALERKIN/` is a Galerkin implementation for the classical variational Poisson problem. To compile and run the problem, type `make`  in the project folder, i.e., `cd problems/MPI_POISSON/GALERKIN; make; ./run.sh`.

## Citing par_hp3D
Please add the following citation to any paper, technical report, presentation, or article that incorporated the `par_hp3D` library:
```bibtex
@book{hpbook2,
      Author = {Demkowicz, L. and Kurtz, J. and Pardo, D. and Paszy\'{n}ski, M. and Rachowicz, W. and Zdunek, A.},
      Publisher = {Chapman \& Hall/CRC},
      Title = {Computing with $hp$ Finite Elements. II. Frontiers: Three-Dimensional Elliptic and Maxwell Problems with Applications},
      Year = {2007}
      }
```
And, optionally,
```bibtex
@book{hpbook,
      Author = {Demkowicz, L.},
      Publisher = {Chapman \& Hall/CRC Press, Taylor and Francis},
      Title = {Computing with $hp$ Finite Elements. I. One- and Two-Dimensional Elliptic and Maxwell Problems},
      Year = {2006}
   }
```
```bibtex
@article{fuentes2015shape,
    Author = {Fuentes, F. and Keith, B. and Demkowicz, L. and Nagaraj, S.},
    Journal = CAMWA,
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
Contact: ``stefan@ices.utexas.edu``, ``leszek@ices.utexas.edu``
