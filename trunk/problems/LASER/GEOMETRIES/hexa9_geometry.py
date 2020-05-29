# ----------------------------------------------------------------------------
#
# geometry.py
#
# ----------------------------------------------------------------------------
# latest revision:  - Jun 2018
#
# purpose:          - generates geometry file for full fiber (core + cladding)
#                     with 9 hexa elements (5 in core, 4 in cladding).
#
# ----------------------------------------------------------------------------

import math

# Length of fiber (z-direction)
LZ = 1.2
# Core and cladding size parameters
core = 0.9
clad = 9.0
r_core = core*math.sqrt(2)
r_clad = clad*math.sqrt(2)
print("Fiber length: %.1f" % LZ)
print("Core radius : %.5f" % r_core)
print("Clad radius : %.5f" % r_clad)
inner = core/4

# Open geometry file
f = open("fiber_hexa9","w+")
# Dimension
f.write("3 3 NDIM,MANDIM\n")
f.write("\n")
# Surfaces
f.write("2 NRSURFS\n")
f.write("\n")
f.write("Cylinder             surface 1: infinite cylinder\n")
f.write("0.0d0 0.0d0 0.0d0    point on the cylinder axis\n")
f.write("0.0d0 0.0d0 1.0d0    cylinder axis vector\n")
f.write("%.15fD0  cylinder radius\n" % r_core)
f.write("\n")
f.write("Cylinder             surface 2: infinite cylinder\n")
f.write("0.0d0 0.0d0 0.0d0    point on the cylinder axis\n")
f.write("0.0d0 0.0d0 1.0d0    cylinder axis vector\n")
f.write("%.15fD0  cylinder radius\n" % r_clad)
f.write("\n")
# Elements
f.write("3 NRDOMAIN\n")
f.write("\n")
# Points
f.write("24 NRPOINT\n")
f.write("\n")
f.write("Regular             point 1\n")
f.write("-%fd0 -%fd0 0.0d0\n" % (inner, inner))
f.write("\n")
f.write("Regular             point 2\n")
f.write("%fd0 -%fd0 0.0d0\n" % (inner, inner))
f.write("\n")
f.write("Regular             point 3\n")
f.write("%fd0 %fd0 0.0d0\n" % (inner, inner))
f.write("\n")
f.write("Regular             point 4\n")
f.write("-%fd0 %fd0 0.0d0\n" % (inner, inner))
f.write("\n")
f.write("Regular             point 5\n")
f.write("-%fd0 -%fd0 0.0d0\n" % (core, core))
f.write("\n")
f.write("Regular             point 6\n")
f.write("%fd0 -%fd0 0.0d0\n" % (core, core))
f.write("\n")
f.write("Regular             point 7\n")
f.write("%fd0 %fd0 0.0d0\n" % (core, core))
f.write("\n")
f.write("Regular             point 8\n")
f.write("-%fd0 %fd0 0.0d0\n" % (core, core))
f.write("\n")
f.write("Regular             point 9\n")
f.write("-%fd0 -%fd0 0.0d0\n" % (clad, clad))
f.write("\n")
f.write("Regular             point 10\n")
f.write("%fd0 -%fd0 0.0d0\n" % (clad, clad))
f.write("\n")
f.write("Regular             point 11\n")
f.write("%fd0 %fd0 0.0d0\n" % (clad, clad))
f.write("\n")
f.write("Regular             point 12\n")
f.write("-%fd0 %fd0 0.0d0\n" % (clad, clad))
f.write("\n")
f.write("Regular             point 13\n")
f.write("-%fd0 -%fd0 %.1fd0\n" % (inner, inner, LZ))
f.write("\n")
f.write("Regular             point 14\n")
f.write("%fd0 -%fd0 %.1fd0\n" % (inner, inner, LZ))
f.write("\n")
f.write("Regular             point 15\n")
f.write("%fd0 %fd0 %.1fd0\n" % (inner, inner, LZ))
f.write("\n")
f.write("Regular             point 16\n")
f.write("-%fd0 %fd0 %.1fd0\n" % (inner, inner, LZ))
f.write("\n")
f.write("Regular             point 17\n")
f.write("-%fd0 -%fd0 %.1fd0\n" % (core, core, LZ))
f.write("\n")
f.write("Regular             point 18\n")
f.write("%fd0 -%fd0 %.1fd0\n" % (core, core, LZ))
f.write("\n")
f.write("Regular             point 19\n")
f.write("%fd0 %fd0 %.1fd0\n" % (core, core, LZ))
f.write("\n")
f.write("Regular             point 20\n")
f.write("-%fd0 %fd0 %.1fd0\n" % (core, core, LZ))
f.write("\n")
f.write("Regular             point 21\n")
f.write("-%fd0 -%fd0 %.1fd0\n" % (clad, clad, LZ))
f.write("\n")
f.write("Regular             point 22\n")
f.write("%fd0 -%fd0 %.1fd0\n" % (clad, clad, LZ))
f.write("\n")
f.write("Regular             point 23\n")
f.write("%fd0 %fd0 %.1fd0\n" % (clad, clad, LZ))
f.write("\n")
f.write("Regular             point 24\n")
f.write("-%fd0 %fd0 %.1fd0\n" % (clad, clad, LZ))
f.write("\n")
# Curves
f.write("52 NRCURVE\n")
f.write("\n")
f.write("Seglin              curve 1\n")
f.write("1 2\n")
f.write("\n")
f.write("Seglin              curve 2\n")
f.write("2 3\n")
f.write("\n")
f.write("Seglin              curve 3\n")
f.write("4 3\n")
f.write("\n")
f.write("Seglin              curve 4\n")
f.write("1 4\n")
f.write("\n")
f.write("Seglin              curve 5\n")
f.write("5 1\n")
f.write("\n")
f.write("Seglin              curve 6\n")
f.write("6 2\n")
f.write("\n")
f.write("Seglin              curve 7\n")
f.write("7 3\n")
f.write("\n")
f.write("Seglin              curve 8\n")
f.write("8 4\n")
f.write("\n")
f.write("QuaCir              curve 9\n")
f.write("5 6\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 10\n")
f.write("6 7\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 11\n")
f.write("8 7\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 12\n")
f.write("5 8\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("Seglin              curve 13\n")
f.write("9 5\n")
f.write("\n")
f.write("Seglin              curve 14\n")
f.write("10 6\n")
f.write("\n")
f.write("Seglin              curve 15\n")
f.write("11 7\n")
f.write("\n")
f.write("Seglin              curve 16\n")
f.write("12 8\n")
f.write("\n")
f.write("QuaCir              curve 17\n")
f.write("9 10\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 18\n")
f.write("10 11\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 19\n")
f.write("12 11\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("QuaCir              curve 20\n")
f.write("9 12\n")
f.write(".0d0, .0d0, .0d0\n")
f.write("\n")
f.write("Seglin              curve 21\n")
f.write("13 14\n")
f.write("\n")
f.write("Seglin              curve 22\n")
f.write("14 15\n")
f.write("\n")
f.write("Seglin              curve 23\n")
f.write("16 15\n")
f.write("\n")
f.write("Seglin              curve 24\n")
f.write("13 16\n")
f.write("\n")
f.write("Seglin              curve 25\n")
f.write("17 13\n")
f.write("\n")
f.write("Seglin              curve 26\n")
f.write("18 14\n")
f.write("\n")
f.write("Seglin              curve 27\n")
f.write("19 15\n")
f.write("\n")
f.write("Seglin              curve 28\n")
f.write("20 16\n")
f.write("\n")
f.write("QuaCir              curve 29\n")
f.write("17 18\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 30\n")
f.write("18 19\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 31\n")
f.write("20 19\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 32\n")
f.write("17 20\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("Seglin              curve 33\n")
f.write("21 17\n")
f.write("\n")
f.write("Seglin              curve 34\n")
f.write("22 18\n")
f.write("\n")
f.write("Seglin              curve 35\n")
f.write("23 19\n")
f.write("\n")
f.write("Seglin              curve 36\n")
f.write("24 20\n")
f.write("\n")
f.write("QuaCir              curve 37\n")
f.write("21 22\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 38\n")
f.write("22 23\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 39\n")
f.write("24 23\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("\n")
f.write("QuaCir              curve 40\n")
f.write("21 24\n")
f.write(".0d0, .0d0, %.1fd0\n" % (LZ))
f.write("Seglin              curve 41\n")
f.write("1 13\n")
f.write("\n")
f.write("Seglin              curve 42\n")
f.write("2 14\n")
f.write("\n")
f.write("Seglin              curve 43\n")
f.write("3 15\n")
f.write("\n")
f.write("Seglin              curve 44\n")
f.write("4 16\n")
f.write("\n")
f.write("Seglin              curve 45\n")
f.write("5 17\n")
f.write("\n")
f.write("Seglin              curve 46\n")
f.write("6 18\n")
f.write("\n")
f.write("Seglin              curve 47\n")
f.write("7 19\n")
f.write("\n")
f.write("Seglin              curve 48\n")
f.write("8 20\n")
f.write("\n")
f.write("Seglin              curve 49\n")
f.write("9 21\n")
f.write("\n")
f.write("Seglin              curve 50\n")
f.write("10 22\n")
f.write("\n")
f.write("Seglin              curve 51\n")
f.write("11 23\n")
f.write("\n")
f.write("Seglin              curve 52\n")
f.write("12 24\n")
f.write("\n")
# Triangles
f.write("0  NRTRIAN\n")
f.write("\n")
# Rectangles
f.write("38 NRRECTA\n")
f.write("\n")
f.write("TraQua              rectangle 1\n")
f.write("1 2 3 4\n")
f.write("\n")
f.write("TraQua              rectangle 2\n")
f.write("5 6 2 1\n")
f.write("\n")
f.write("TraQua              rectangle 3\n")
f.write("6 7 3 2\n")
f.write("\n")
f.write("TraQua              rectangle 4\n")
f.write("8 7 3 4\n")
f.write("\n")
f.write("TraQua              rectangle 5\n")
f.write("5 8 4 1\n")
f.write("\n")
f.write("TraQua              rectangle 6\n")
f.write("9 10 6 5\n")
f.write("\n")
f.write("TraQua              rectangle 7\n")
f.write("10 11 7 6\n")
f.write("\n")
f.write("TraQua              rectangle 8\n")
f.write("12 11 7 8\n")
f.write("\n")
f.write("TraQua              rectangle 9\n")
f.write("9 12 8 5\n")
f.write("\n")
f.write("TraQua              rectangle 10\n")
f.write("13 14 15 16\n")
f.write("\n")
f.write("TraQua              rectangle 11\n")
f.write("17 18 14 13\n")
f.write("\n")
f.write("TraQua              rectangle 12\n")
f.write("18 19 15 14\n")
f.write("\n")
f.write("TraQua              rectangle 13\n")
f.write("20 19 15 16\n")
f.write("\n")
f.write("TraQua              rectangle 14\n")
f.write("17 20 16 13\n")
f.write("\n")
f.write("TraQua              rectangle 15\n")
f.write("21 22 18 17\n")
f.write("\n")
f.write("TraQua              rectangle 16\n")
f.write("22 23 19 18\n")
f.write("\n")
f.write("TraQua              rectangle 17\n")
f.write("24 23 19 20\n")
f.write("\n")
f.write("TraQua              rectangle 18\n")
f.write("21 24 20 17\n")
f.write("\n")
f.write("TraQua              rectangle 19\n")
f.write("1 2 14 13\n")
f.write("\n")
f.write("TraQua              rectangle 20\n")
f.write("2 3 15 14\n")
f.write("\n")
f.write("TraQua              rectangle 21\n")
f.write("4 3 15 16\n")
f.write("\n")
f.write("TraQua              rectangle 22\n")
f.write("1 4 16 13\n")
f.write("\n")
f.write("TraQua              rectangle 23\n")
f.write("5 1 13 17\n")
f.write("\n")
f.write("TraQua              rectangle 24\n")
f.write("6 2 14 18\n")
f.write("\n")
f.write("TraQua              rectangle 25\n")
f.write("7 3 15 19\n")
f.write("\n")
f.write("TraQua              rectangle 26\n")
f.write("8 4 16 20\n")
f.write("\n")
f.write("TraQua              rectangle 27\n")
f.write("5 6 18 17\n")
f.write("\n")
f.write("TraQua              rectangle 28\n")
f.write("6 7 19 18\n")
f.write("\n")
f.write("TraQua              rectangle 29\n")
f.write("8 7 19 20\n")
f.write("\n")
f.write("TraQua              rectangle 30\n")
f.write("5 8 20 17\n")
f.write("\n")
f.write("TraQua              rectangle 31\n")
f.write("9 5 17 21\n")
f.write("\n")
f.write("TraQua              rectangle 32\n")
f.write("10 6 18 22\n")
f.write("\n")
f.write("TraQua              rectangle 33\n")
f.write("11 7 19 23\n")
f.write("\n")
f.write("TraQua              rectangle 34\n")
f.write("12 8 20 24\n")
f.write("\n")
f.write("TraQua              rectangle 35\n")
f.write("9 10 22 21\n")
f.write("\n")
f.write("TraQua              rectangle 36\n")
f.write("10 11 23 22\n")
f.write("\n")
f.write("TraQua              rectangle 37\n")
f.write("12 11 23 24\n")
f.write("\n")
f.write("TraQua              rectangle 38\n")
f.write("9 12 24 21\n")
f.write("\n")
# Prisms
f.write("0 NRPRISM\n")
f.write("\n")
# Hexahedra
f.write("9 NRHEXAS\n")
f.write("\n")
f.write("TraHex              hexa 1\n")
f.write("1    1 2 3 4 13 14 15 16\n")
f.write("\n")
f.write("TraHex              hexa 2\n")
f.write("2    5 6 2 1 17 18 14 13\n")
f.write("\n")
f.write("TraHex              hexa 3\n")
f.write("2    6 7 3 2 18 19 15 14\n")
f.write("\n")
f.write("TraHex              hexa 4\n")
f.write("2    8 7 3 4 20 19 15 16\n")
f.write("\n")
f.write("TraHex              hexa 5\n")
f.write("2    5 8 4 1 17 20 16 13\n")
f.write("\n")
f.write("TraHex              hexa 6\n")
f.write("3    9 10 6 5 21 22 18 17\n")
f.write("\n")
f.write("TraHex              hexa 7\n")
f.write("3    10 11 7 6 22 23 19 18\n")
f.write("\n")
f.write("TraHex              hexa 8\n")
f.write("3    12 11 7 8 24 23 19 20\n")
f.write("\n")
f.write("TraHex              hexa 9\n")
f.write("3    9 12 8 5 21 24 20 17\n")
f.write("\n")
# Tetrahedra
f.write("0 NRTETRA\n")
f.write("\n")
# Pyramids
f.write("0 NRPYRAM\n")
# Close geometry file
f.close()
