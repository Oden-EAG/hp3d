# ----------------------------------------------------------------------------
#
# generate_bent_slab_Rcenter.py
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
a = 0.5;
b = 5.0;

# R = 1300.0;
# keff= 2.17655321571770E+02-2.17455000000000E+02;
# f = open("bentstepslab_R1300_4wl_m00","w+")

R = 2600.0;
keff= 2.17543433386383E+02-2.1743000000000E+02;
f = open("verif_bentstepslab_R2600","w+")

# R = 2600.0;
# keff= 2.17655321571770E+02-2.17455000000000E+02;
# f = open("bentstepslab_R2600_4wl_m00","w+")

# R = 2600.0;
# keff= 2.17543433386383E+02-2.17480E+02;
# f = open("bentstepslab_R2600_4wl_m02","w+")

# R = 2600.0;
# keff= 2.17611E+02-2.17480E+02;
# f = open("verif_bentstepslab_R2600_m01","w+")

# R = 1300.0;
# keff= 2.17611E+02-2.17480E+02;
# f = open("verif_bentstepslab_R1300_m01","w+")


wleff=2.0*math.pi/keff;
nwl_bent = 4;
TH = nwl_bent*wleff/R;
print("wleff=",wleff)
print("theta_end=",TH)

# TH = 3.0*math.pi/180.0; 

# compute cosine and sine of angle
cth = math.cos(TH); sth = math.sin(TH);
# and open geometry file

# Dimension
f.write("3 3 NDIM,MANDIM\n")
f.write("\n")
# Surfaces
f.write("0 NRSURFS\n")
f.write("\n")
# Materials
f.write("3 NRDOMAIN\n")
f.write("\n")
# Boundary domains flag
f.write("1  ISURF_FLAG\n")
f.write("\n")
# Points
f.write("16 NRPOINT\n")
f.write("\n")
f.write("Regular             point 1\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R-b, 0.0))
f.write("\n")
f.write("Regular             point 2\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, R-b, 0.0))
f.write("\n")
f.write("Regular             point 3\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R-a, 0.0))
f.write("\n")
f.write("Regular             point 4\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, R-a, 0.0))
f.write("\n")
f.write("Regular             point 5\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R+a, 0.0))
f.write("\n")
f.write("Regular             point 6\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, R+a, 0.0))
f.write("\n")
f.write("Regular             point 7\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, R+b, 0.0))
f.write("\n")
f.write("Regular             point 8\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, R+b, 0.0))
f.write("\n")
f.write("Regular             point 9\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R-b)*cth, (R-b)*sth))
f.write("\n")
f.write("Regular             point 10\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, (R-b)*cth, (R-b)*sth))
f.write("\n")
f.write("Regular             point 11\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R-a)*cth, (R-a)*sth))
f.write("\n")
f.write("Regular             point 12\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, (R-a)*cth, (R-a)*sth))
f.write("\n")
f.write("Regular             point 13\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R+a)*cth, (R+a)*sth))
f.write("\n")
f.write("Regular             point 14\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, (R+a)*cth, (R+a)*sth))
f.write("\n")
f.write("Regular             point 15\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (-a, (R+b)*cth, (R+b)*sth))
f.write("\n")
f.write("Regular             point 16\n")
f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % ( a, (R+b)*cth, (R+b)*sth))
f.write("\n")
# Curves
f.write("28 NRCURVE\n")
f.write("\n")
f.write("Seglin              curve 1\n")
f.write("1 2\n")
f.write("\n")
f.write("Seglin              curve 2\n")
f.write("3 4\n")
f.write("\n")
f.write("Seglin              curve 3\n")
f.write("5 6\n")
f.write("\n")
f.write("Seglin              curve 4\n")
f.write("7 8\n")
f.write("\n")
f.write("Seglin              curve 5\n")
f.write("1 3\n")
f.write("\n")
f.write("Seglin              curve 6\n")
f.write("2 4\n")
f.write("\n")
f.write("Seglin              curve 7\n")
f.write("3 5\n")
f.write("\n")
f.write("Seglin              curve 8\n")
f.write("4 6\n")
f.write("\n")
f.write("Seglin              curve 9\n")
f.write("5 7\n")
f.write("\n")
f.write("Seglin              curve 10\n")
f.write("6 8\n")
f.write("\n")
f.write("Seglin              curve 11\n")
f.write("9 10\n")
f.write("\n")
f.write("Seglin              curve 12\n")
f.write("11 12\n")
f.write("\n")
f.write("Seglin              curve 13\n")
f.write("13 14\n")
f.write("\n")
f.write("Seglin              curve 14\n")
f.write("15 16\n")
f.write("\n")
f.write("Seglin              curve 15\n")
f.write("9 11\n")
f.write("\n")
f.write("Seglin              curve 16\n")
f.write("10 12\n")
f.write("\n")
f.write("Seglin              curve 17\n")
f.write("11 13\n")
f.write("\n")
f.write("Seglin              curve 18\n")
f.write("12 14\n")
f.write("\n")
f.write("Seglin              curve 19\n")
f.write("13 15\n")
f.write("\n")
f.write("Seglin              curve 20\n")
f.write("14 16\n")
f.write("\n")
f.write("CylCoord            curve 21\n")
f.write("1 9\n")
f.write("\n")
f.write("CylCoord            curve 22\n")
f.write("2 10\n")
f.write("\n")
f.write("CylCoord            curve 23\n")
f.write("3 11\n")
f.write("\n")
f.write("CylCoord            curve 24\n")
f.write("4 12\n")
f.write("\n")
f.write("CylCoord            curve 25\n")
f.write("5 13\n")
f.write("\n")
f.write("CylCoord            curve 26\n")
f.write("6 14\n")
f.write("\n")
f.write("CylCoord            curve 27\n")
f.write("7 15\n")
f.write("\n")
f.write("CylCoord            curve 28\n")
f.write("8 16\n")
f.write("\n")
# Triangles
f.write("0  NRTRIAN\n")
f.write("\n")
# Rectangles
f.write("16 NRRECTA\n")
f.write("\n")
f.write("BilQua              rectangle 1\n")     # bottom face
f.write("0   1 2 4 3\n")
f.write("\n")
f.write("BilQua              rectangle 2\n")     # top face
f.write("0   9 10 12 11\n")
f.write("\n")
f.write("BilQua              rectangle 3\n")     # bottom face
f.write("0   3 4 6 5\n")
f.write("\n")
f.write("BilQua              rectangle 4\n")     # top face
f.write("0   11 12 14 13\n")
f.write("\n")
f.write("BilQua              rectangle 5\n")     # bottom face
f.write("0   5 6 8 7\n")
f.write("\n")
f.write("BilQua              rectangle 6\n")     # top face
f.write("0   13 14 16 15\n")
f.write("\n")
f.write("CylRec              rectangle 7\n")     # inner radius boundary
f.write("1   1 2 10 9\n")
f.write("\n")
f.write("CylRec              rectangle 8\n")     # left cladding-core interface
f.write("0   3 4 12 11\n")
f.write("\n")
f.write("CylRec              rectangle 9\n")     # right cladding-core interface
f.write("0   5 6 14 13\n")
f.write("\n")
f.write("CylRec              rectangle 10\n")     # outer radius boundary
f.write("0   7 8 16 15\n")
f.write("\n")
f.write("CylRec              rectangle 11\n")     # normal to x
f.write("0   1 3 11 9\n")
f.write("\n")
f.write("CylRec              rectangle 12\n")     # normal to x
f.write("0   2 4 12 10\n")
f.write("\n")
f.write("CylRec              rectangle 13\n")     # normal to x
f.write("0   3 5 13 11\n")
f.write("\n")
f.write("CylRec              rectangle 14\n")     # normal to x
f.write("0   4 6 14 12\n")
f.write("\n")
f.write("CylRec              rectangle 15\n")     # normal to x
f.write("0   5 7 15 13\n")
f.write("\n")
f.write("CylRec              rectangle 16\n")     # normal to x
f.write("0   6 8 16 14\n")
f.write("\n")
# Prisms
f.write("0 NRPRISM\n")
f.write("\n")
# Hexahedra
f.write("3 NRHEXAS\n")
f.write("\n")
f.write("CylHex              hexa 1\n")
f.write("2    1 2 4 3 9 10 12 11\n")
f.write("\n")
f.write("CylHex              hexa 2\n")
f.write("1    3 4 6 5 11 12 14 13\n")
f.write("\n")
f.write("CylHex              hexa 3\n")
f.write("3    5 6 8 7 13 14 16 15\n")
f.write("\n")
# Tetrahedra
f.write("0 NRTETRA\n")
f.write("\n")
# Pyramids
f.write("0 NRPYRAM\n")
# Close geometry file
f.close()
