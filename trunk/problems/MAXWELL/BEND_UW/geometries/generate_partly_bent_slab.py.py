# ----------------------------------------------------------------------------
#
# generate_partly_bent_slab.py
#
# ----------------------------------------------------------------------------
# latest revision:  - Oct 2024
#
# purpose:          - generates geometry file for bent slab waveguide, 
#                     with bending radius R at a square of side length 2a
#
# ----------------------------------------------------------------------------

import math

# Set side half length, bending radius of cross section's center 
# and spanning angle in degrees
a = 2.0;
R = 210.0*a;
TH = 30.0; 
sl = R*TH*math.pi/180.0;     # straight portion length
# compute cosine and sine of angle
cth = math.cos(TH*math.pi/180.0); sth = math.sin(TH*math.pi/180.0);
# and open geometry file
f = open("partly_bent_slab_d2_cR420_DEG30","w+")

# Dimension
f.write("3 3 NDIM,MANDIM\n")
f.write("\n")
# Surfaces
f.write("0 NRSURFS\n")
f.write("\n")
# Materials
f.write("1 NRDOMAIN\n")
f.write("\n")
# Boundary domains flag
f.write("1  ISURF_FLAG\n")
f.write("\n")
# Points
f.write("12 NRPOINT\n")
f.write("\n")
f.write("Regular             point 1\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R-a, 0.0))
f.write("\n")
f.write("Regular             point 2\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, R-a, 0.0))
f.write("\n")
f.write("Regular             point 3\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, R+a, 0.0))
f.write("\n")
f.write("Regular             point 4\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R+a, 0.0))
f.write("\n")
f.write("Regular             point 5\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R-a)*cth, (R-a)*sth))
f.write("\n")
f.write("Regular             point 6\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, (R-a)*cth, (R-a)*sth))
f.write("\n")
f.write("Regular             point 7\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, (R+a)*cth, (R+a)*sth))
f.write("\n")
f.write("Regular             point 8\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R+a)*cth, (R+a)*sth))
f.write("\n")
f.write("Regular             point 9\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R-a, -sl))
f.write("\n")
f.write("Regular             point 10\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, R-a, -sl))
f.write("\n")
f.write("Regular             point 11\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (a, R+a, -sl))
f.write("\n")
f.write("Regular             point 12\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R+a, -sl))
f.write("\n")

# Curves
f.write("20 NRCURVE\n")
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
f.write("5 6\n")
f.write("\n")
f.write("Seglin              curve 6\n")
f.write("6 7\n")
f.write("\n")
f.write("Seglin              curve 7\n")
f.write("8 7\n")
f.write("\n")
f.write("Seglin              curve 8\n")
f.write("5 8\n")
f.write("\n")
f.write("CylCoord            curve 9\n")
f.write("1 5\n")
f.write("\n")
f.write("CylCoord            curve 10\n")
f.write("2 6\n")
f.write("\n")
f.write("CylCoord            curve 11\n")
f.write("3 7\n")
f.write("\n")
f.write("CylCoord            curve 12\n")
f.write("4 8\n")
f.write("\n")
f.write("Seglin              curve 13\n")
f.write("9 10\n")
f.write("\n")
f.write("Seglin              curve 14\n")
f.write("10 11\n")
f.write("\n")
f.write("Seglin              curve 15\n")
f.write("12 11\n")
f.write("\n")
f.write("Seglin              curve 16\n")
f.write("9 12\n")
f.write("\n")
f.write("Seglin              curve 17\n")
f.write("9 1\n")
f.write("\n")
f.write("Seglin              curve 18\n")
f.write("10 2\n")
f.write("\n")
f.write("Seglin              curve 19\n")
f.write("11 3\n")
f.write("\n")
f.write("Seglin              curve 20\n")
f.write("12 4\n")
f.write("\n")
# Triangles
f.write("0  NRTRIAN\n")
f.write("\n")
# Rectangles
f.write("11 NRRECTA\n")
f.write("\n")
f.write("BilQua              rectangle 1\n")     # bent-straight interface
f.write("0   1 2 3 4\n")
f.write("\n")
f.write("BilQua              rectangle 2\n")     # top face
f.write("2   5 6 7 8\n")
f.write("\n")
f.write("CylRec              rectangle 3\n")
f.write("0   1 2 6 5\n")
f.write("\n")
f.write("CylRec              rectangle 4\n")     # normal to x
f.write("1   2 3 7 6\n")
f.write("\n")
f.write("CylRec              rectangle 5\n")
f.write("0   4 3 7 8\n")
f.write("\n")
f.write("CylRec              rectangle 6\n")     # normal to x
f.write("1   1 4 8 5\n")
f.write("\n")
f.write("BilQua              rectangle 7\n")     # bottom face
f.write("0   9 10 11 12\n")
f.write("\n")
f.write("BilQua              rectangle 8\n")
f.write("0   9 10 2 1\n")
f.write("\n")
f.write("BilQua              rectangle 9\n")     # normal to x
f.write("1   10 11 3 2\n")
f.write("\n")
f.write("BilQua              rectangle 10\n")
f.write("0   12 11 3 4\n")
f.write("\n")
f.write("BilQua              rectangle 11\n")     # normal to x
f.write("1   9 12 4 1\n")
f.write("\n")
# Prisms
f.write("0 NRPRISM\n")
f.write("\n")
# Hexahedra
f.write("2 NRHEXAS\n")
f.write("\n")
f.write("CylHex              hexa 1\n")
f.write("1    1 2 3 4 5 6 7 8\n")
f.write("\n")
f.write("Linear              hexa 2\n")
f.write("1    9 10 11 12 1 2 3 4\n")
f.write("\n")
# Tetrahedra
f.write("0 NRTETRA\n")
f.write("\n")
# Pyramids
f.write("0 NRPYRAM\n")
# Close geometry file
f.close()
