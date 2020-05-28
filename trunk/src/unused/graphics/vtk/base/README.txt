This is the new code for using VTK to visualize graphics output files
produced by hp3d. Everything is in C++. The code is structured in the
form of two classes:

The implementation of class geometry is separated into the source files
geometry.cxx, geometry_hexa.cxx, geometry_tetra.cxx, geometry_prism.cxx
and geometry_pyramid.cxx, with the prototype in the header file
geometry.h. This class is used to render mesh geometry.

Class solution (solution.cxx and solution.h) inherits much of the
functionality of class geometry, and adds the ability to read and
display components of the solution.

There are some simple utility routines defined in util.h and implemented
in util.cxx. The main program (main.cxx) allows the user to select
between geometry and solution graphics, opens the appropriate file, and
invokes methods of the appropriate class.
