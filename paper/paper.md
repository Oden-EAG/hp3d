---
title: '$hp\mathrm{3D}$: A Scalable MPI/OpenMP $hp$-Adaptive Finite Element Software Library for Complex Multiphysics Applications'
tags:
  - finite element
  - hp-adaptivity
  - DPG method
authors:
  - name: Stefan Henneking
    orcid: 0000-0003-2177-8519
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: Jacob Badger
    orcid: 0000-0001-6482-105X
    equal-contrib: true
    affiliation: 1
  - name: Socratis Petrides
    orcid: 0000-0002-1284-5495
    equal-contrib: true
    affiliation: 2
  - name: Federico Fuentes
    orcid: 0000-0002-4039-082X
    equal-contrib: true
    affiliation: 3
  - name: Leszek Demkowicz
    orcid: 0000-0001-7839-8037
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Oden Institute, The University of Texas at Austin, USA
   index: 1
 - name: Lawrence Livermore National Laboratory, USA
   index: 2
 - name: Institute for Mathematical and Computational Engineering, Pontificia Universidad Católica de Chile, Chile
   index: 3
date: 2 October 2023
bibliography: paper.bib

---

# Summary

The $hp\mathrm{3D}$ finite element (FE) library is a tool for computational modeling of engineering applications. The library provides a framework for discretization of three-dimensional multiphysics problems described by systems of partial differential equations. $hp\mathrm{3D}$ can be compiled in real or complex mode to accommodate the need for real- or complex-valued physics variables, respectively. The library is written entirely in Fortran; user applications interfacing with the library are also written in Fortran. The $hp\mathrm{3D}$ software can be installed and runs efficiently on various compute architectures, from laptops and single workstations to state-of-the-art supercomputers.

# Statement of need

$hp\mathrm{3D}$ combines a list of unique features and algorithms setting it apart from other publicly available finite element libraries. The software supports hybrid meshes combining elements of "all shapes:" hexahedra, tetrahedra, prisms, and pyramids. The internal Geometric Modeling Package (GMP) provides support for \emph{exact geometry} elements and \emph{isoparametric} elements. Exact geometry elements directly use the parametrizations provided by GMP resulting in computations with no geometry error; isoparametric elements approximate the geometry maps with polynomials spanning the element space of shape functions. $hp\mathrm{3D}$'s shape functions package provides compatible discretization of energy spaces forming the $H^1-H(\text{curl})-H(\text{div})-L^2$ exact sequence [@fuentes2015shape; @demkowicz2023fem]. Additionally, the $hp\mathrm{3D}$ FE code sets itself apart from other advanced FE libraries (e.g., MFEM [@anderson2021mfem] or deal.II [@bangerth2007deal]) by focusing on \emph{$hp$-adaptive} solutions; in particular, the code supports \emph{anisotropic refinements} in both element size $h$ and polynomial order $p$. Such $hp$-adaptive methods are the most efficient way to converge to difficult solutions, e.g., resolving solutions with boundary layers, adapting toward geometric singularities, etc. [@chakraborty2023hp]. $hp\mathrm{3D}$ features a number of unique algorithms, including \emph{constrained approximation} routines for assembling elements with hanging nodes [@demkowicz1989-hp1; @demkowicz1989-hp2; @demkowicz1989-hp3], and \emph{projection-based interpolation} for computation of nodal constraints [@demkowicz2008interp]. Besides discretization with the standard Bubnov-Galerkin FE method, the $hp\mathrm{3D}$ library also supports discretization with the discontinuous Petrov-Galerkin (DPG) method [@demkowicz2017dpg]. The $hp\mathrm{3D}$ software leverages hybrid MPI/OpenMP parallelism to run efficiently on large-scale computing facilities [@henneking2021phd; @badger2023scalable].

# Dependencies

$hp\mathrm{3D}$ interfaces with several well-established third-party libraries: for mesh partitioning (ParMETIS [@parmetis], PT-Scotch [@pt-scotch]), for dynamic load balancing (Zoltan [@ZoltanOverviewArticle2002]), for linear solvers (MUMPS [@amestoy2001mumps], PETSc [@petsc-user-ref]), and for I/O (pHDF5 [@hdf5]). We note that all of these dependencies can be directly installed via PETSc. Additionally, the geometry mesh and the solution can be exported to VTK and visualized with ParaView [@ahrens2005paraview].

# Examples of applications

Some examples of 3D applications that have been implemented in $hp\mathrm{3D}$ include modeling of acoustic wave propagation in the human head [@gatto2012phd], modeling of electromagnetic waves with thermal effects in the human head [@kim2013phd], electromagnetic and acoustic scattering and high-frequency beams [@petrides2019phd; @petrides2021adaptive], modeling of insulators in high-energy density electric motors (thermo-viscoelasticity) [@fuentes2017viscoelasticity], and modeling of optical amplifiers (nonlinear Maxwell equations coupled with the heat equation) [@nagaraj2018raman; @henneking2021fiber; @henneking2022parallel].

# Further reading

Instructions on installing and using the code are available in the $hp\mathrm{3D}$ user manual [@hpUserManual]. $hp\mathrm{3D}$'s underlying algorithms are described in various published articles and books. Many parts of the current version of the $hp\mathrm{3D}$ software are based on the algorithms described in the two-volume $hp$ book series [@hpbook; @hpbook2] on the former 2D and 3D versions of the code (which were not published as open-source libraries). A third $hp$ book volume detailing the additions and modifications of the newest version, including its MPI/OpenMP parallel algorithms, will be published soon [@hpbook3]. The orientation-embedded shape functions package is described in [@fuentes2015shape]. Details on the FE methodology and conforming discretization of exact-sequence elements are given in [@demkowicz2023fem].

# Acknowledgements

We acknowledge contributions from Ankit Chakraborty, Paolo Gatto, Brendan Keith, Kyungjoo Kim, Jaime D. Mora, and Sriram Nagaraj for this project. The development and open-sourcing of the $hp\mathrm{3D}$ FE code are supported by NSF award 2103524. This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344, LLNL-JRNL-855288.

# References

