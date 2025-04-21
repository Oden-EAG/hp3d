# ----------------------------------------------------------------------------
#
# generate_bent_fiber.py
#
# ----------------------------------------------------------------------------
# latest revision:  - Nov 2024
#
# purpose:          - generates geometry file for bent optical fiber
#
# ----------------------------------------------------------------------------

import math
import numpy as np

def toroidal2cartesian(rmaj,rmin,phi,theta):
	x = rmin*math.sin(phi)
	y = (rmaj+rmin*math.cos(phi))*math.cos(theta)
	z = (rmaj+rmin*math.cos(phi))*math.sin(theta)
	return x, y, z

def print_point(x,y,z,np):
	f.write("Regular            point %3i \n" % np)
	f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (x, y, z))
	f.write("\n")

def print_curve(type,p1,p2,nc,cc=[0.0, 0.0, 0.0]):
	f.write("%-8s              curve %3i\n" % (type,nc))
	f.write("%3i %3i\n" % (p1, p2))
	if (type=='QuaCir'):
		f.write("%18.14fd0 %18.14fd0 %18.14fd0\n" % (cc[0], cc[1], cc[2]) )
	f.write("\n")

def print_triangle(type,p1,p2,p3,nt):
	f.write("%-8s              triangle %3i\n" % (type,nt))
	f.write("%3i %3i %3i\n" % (p1, p2, p3))
	f.write("\n")

def print_rectangle(type,p1,p2,p3,p4,nr,cl=[0.0, 0.0, 0.0],cu=[0.0, 0.0, 0.0]):
	f.write("%-8s              rectangle %3i\n" % (type,nr))
	f.write("%3i %3i %3i %3i\n" % (p1, p2, p3, p4))
	if (type=='TorRec'):
		f.write("%18.14fd0 %18.14fd0 %18.14fd0\t\t center of arc v1v2\n" % (cl[0], cl[1], cl[2]) )
		f.write("%18.14fd0 %18.14fd0 %18.14fd0\t\t center of arc v3v4\n" % (cu[0], cu[1], cu[2]) )
	f.write("\n")

def print_prism(type,ndom,p1,p2,p3,p4,p5,p6,npri):
	f.write("%-8s              prism %3i\n" % (type,npri))
	f.write("%1i \t %3i %3i %3i %3i %3i %3i\n" % (ndom, p1, p2, p3, p4, p5, p6))
	f.write("\n")

def print_hexa(type,ndom,p1,p2,p3,p4,p5,p6,p7,p8,nh,nro=0,nri=0):
	f.write("%-8s              hexahedron %3i\n" % (type,nh))
	f.write("%1i \t %3i %3i %3i %3i %3i %3i %3i %3i\n" % (ndom, p1, p2, p3, p4, p5, p6, p7, p8))
	if (type=='TorHex'):
		f.write("%3i %3i\t\t indices of outer and inner toroidal quad faces\n" % (nro,nri))
	f.write("\n")




# Set fiber radii, bending radius of cross section's center 
# and spanning angle in radians (start angle assumed 0)
r_core = 0.5;
r_prism = 0.5*r_core;
r_inner_clad = 10.0*r_core;
r_outer_clad = 20.0*r_core;
R = 2600.0*r_core;

keff=1.100582; # for mode LP01 and usual data
wleff=2.0*math.pi/keff;
nwl_bent = 6;
nwl_strt = 2;
theta_end = nwl_bent*wleff/R;
print("wleff=",wleff)
print("theta_end=",theta_end)
#
# angle (wrt plane yz) to place first point in cross section. 
# Typical values 0 or pi/4
phi1 = math.pi/4.0; 
hexlayers_core = 1; # don't include the inner prisms layer
hexlayers_inner_clad = 1; # 
hexlayers_outer_clad = 1; # 
theta_subdiv = 3;



dr_core = (r_core-r_prism)/hexlayers_core
dr_inner_clad = (r_inner_clad-r_core)/max(hexlayers_inner_clad,1)
dr_outer_clad = (r_outer_clad-r_inner_clad)/max(hexlayers_outer_clad,1)

# and open geometry file
f = open("partly_bent_fiber_test_4wl","w+")

# Dimension
f.write("3 3 NDIM,MANDIM\n")
f.write("\n")
# Surfaces
f.write("0 NRSURFS\n")
f.write("\n")
# Materials
f.write("4 NRDOMAIN\n")
f.write("\n")


# Points

# adjust subdivisions if there is a straight portion
if nwl_strt>0:
	theta_subdiv += 1
	zp = -nwl_strt*wleff;


# nr of points per cross section, without the one in the axis
npcs = 4*(1+hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad)
# total nr of points
nr_point = (1+npcs)*(theta_subdiv+1)


f.write("%-3i NRPOINT\n" % nr_point)
f.write("\n")

ip = 0
cc = np.zeros((theta_subdiv+1,3))
ci = np.zeros(theta_subdiv+1,dtype=int)

for div in range(theta_subdiv+1):
	if nwl_strt>0:
		if div==0:
			th = 0
		else:
			th = (div-1)*theta_end/(theta_subdiv-1)
	else:
		th = div*theta_end/theta_subdiv
	# prism layer in core
	r = r_prism
	for t in range(4):
		ip += 1
		phi = phi1 + t*math.pi/2.0
		x,y,z = toroidal2cartesian(R,r,phi,th); 
		if nwl_strt>0 and div==0:
			z = zp
		print_point(x,y,z,ip)
	# hex layers in core
	for l in range(hexlayers_core):
		r = r_prism + dr_core*(l+1)
		for t in range(4):
			ip += 1
			phi = phi1 + t*math.pi/2.0
			x,y,z = toroidal2cartesian(R,r,phi,th); 
			if nwl_strt>0 and div==0:
				z = zp
			print_point(x,y,z,ip)
	# hex layers in inner clad
	for l in range(hexlayers_inner_clad):
		r = r_core + dr_inner_clad*(l+1)
		for t in range(4):
			ip += 1
			phi = phi1 + t*math.pi/2.0
			x,y,z = toroidal2cartesian(R,r,phi,th); 
			if nwl_strt>0 and div==0:
				z = zp
			print_point(x,y,z,ip)
	# hex layers in outer clad
	for l in range(hexlayers_outer_clad):
		r = r_inner_clad + dr_outer_clad*(l+1)
		for t in range(4):
			ip += 1
			phi = phi1 + t*math.pi/2.0
			x,y,z = toroidal2cartesian(R,r,phi,th);
			if nwl_strt>0 and div==0:
				z = zp
			print_point(x,y,z,ip)
	# 
# points in cross section's axis  
for div in range(theta_subdiv+1):
	if nwl_strt>0:
		if div==0:
			th = 0
		else:
			th = (div-1)*theta_end/(theta_subdiv-1)
	else:
		th = div*theta_end/theta_subdiv
	r = 0.0; phi = 0.0
	ip += 1
	x,y,z = toroidal2cartesian(R,r,phi,th);
	if nwl_strt>0 and div==0:
		z = zp
	print_point(x,y,z,ip)
	# save in memory coordinates and index of this center point
	cc[div] = [x,y,z]
	ci[div] = ip

	# print("div,ci[div],cc[div]=",div,ci[div],cc[div])
	#


# Curves
# nr of seglin per cross section
nseglin_cs = 4*(1+hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad)
# nr of quacir per cross section
nquacir_cs = nseglin_cs
# nr of curves joining two consecutive cross sections, cylcoord
npcyl = 1+nquacir_cs
# total nr of curves
nr_curve = (nseglin_cs+nquacir_cs)*(theta_subdiv+1) + npcyl*theta_subdiv
f.write("%-3i NRCURVE\n" % nr_curve)
f.write("\n")
ic = 0
for div in range(theta_subdiv+1):
	# prism layer in core
	for t in range(1,5):
		ic += 1
		p1 = npcs*div + t
		p2 = ci[div]
		print_curve('Seglin',p1,p2,ic)
	for t in range(1,5):
		ic += 1
		p1 = npcs*div + t
		p2 = npcs*div + t%4 +1
		print_curve('QuaCir',p1,p2,ic,cc[div])
	# hex layers in core
	for l in range(hexlayers_core):
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+l) + t
			p2 = p1 - 4
			print_curve('Seglin',p1,p2,ic)
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+l) + t
			p2 = npcs*div + 4*(1+l) + t%4 +1
			print_curve('QuaCir',p1,p2,ic,cc[div])
	# hex layers in inner clad
	for l in range(hexlayers_inner_clad):
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+hexlayers_core+l) + t
			p2 = p1 - 4
			print_curve('Seglin',p1,p2,ic)
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+hexlayers_core+l) + t
			p2 = npcs*div + 4*(1+hexlayers_core+l) + t%4 +1
			print_curve('QuaCir',p1,p2,ic,cc[div])
	# hex layers in outer clad
	for l in range(hexlayers_inner_clad):
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+hexlayers_core+hexlayers_inner_clad+l) + t
			p2 = p1 - 4
			print_curve('Seglin',p1,p2,ic)
		for t in range(1,5):
			ic += 1
			p1 = npcs*div + 4*(1+hexlayers_core+hexlayers_inner_clad+l) + t
			p2 = npcs*div + 4*(1+hexlayers_core+hexlayers_inner_clad+l) + t%4 +1
			print_curve('QuaCir',p1,p2,ic,cc[div])
# curves between cross sections
for div in range(theta_subdiv):
	for p in range(1,npcs+1):
		ic += 1
		p1 = p + npcs*div
		p2 = p1+ npcs
		if nwl_strt>0 and div==0:
			print_curve('Seglin',p1,p2,ic)	
		else:
			print_curve('CylCoord',p1,p2,ic)
# curves of central axis
for div in range(theta_subdiv):
	ic += 1
	p1 = ci[div]
	p2 = ci[div+1]
	if nwl_strt>0 and div==0:
		print_curve('Seglin',p1,p2,ic)	
	else:
		print_curve('CylCoord',p1,p2,ic)


# triangles
nr_trian = 4*(theta_subdiv+1)
f.write("%-2i         NRTRIAN\n" % nr_trian)
f.write("\n")
# 
it = 0
for div in range(theta_subdiv+1):
	for t in range(1,5):
		it += 1
		p1 = npcs*div + t
		p2 = ci[div]
		p3 = npcs*div + t%4 + 1
		print_triangle('TransTri',p1,p2,p3,it)

# rectangles
# nr of transfinite rectangles per cross section
nqtracs = 4*(hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad)
# nr of cylindrical mapped rectangles per subdivision
nqcyldiv = nqtracs + 4
# nr of toroidal mapped rectangles per subdivision
nqtordiv = nqcyldiv
# total
nr_quad = nqtracs*(theta_subdiv+1) + (nqcyldiv+nqtordiv)*theta_subdiv
f.write("%-3i         NRRECTA\n" % nr_quad)
f.write("\n")
ir = 0
# first, the toroidal rectangles
for div in range(theta_subdiv):
	for l in range(1+hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad):
		for t in range(1,5):
			ir += 1
			p1 = npcs*div +4*l + t
			p2 = npcs*div +4*l + t%4 +1
			p3 = npcs*(div+1) +4*l + t%4 +1
			p4 = npcs*(div+1) +4*l + t
			if nwl_strt>0 and div==0:
				print_rectangle('TraQua',p1,p2,p3,p4,ir)
			else:
				print_rectangle('TorRec',p1,p2,p3,p4,ir,cc[div],cc[div+1])
# now, the cylindrical rectangles
for div in range(theta_subdiv):
	for t in range(1,5):
		ir += 1
		p1 = ci[div]
		p2 = npcs*div + t
		p3 = npcs*(div+1) + t
		p4 = ci[div+1]
		if nwl_strt>0 and div==0:
			print_rectangle('BilQua',p1,p2,p3,p4,ir)
		else:
			print_rectangle('CylRec',p1,p2,p3,p4,ir)
	for l in range(hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad):
		for t in range(1,5):
			ir += 1
			p1 = npcs*div +4*l + t
			p2 = npcs*div +4*(l+1) + t
			p3 = npcs*(div+1) +4*(l+1) + t
			p4 = npcs*(div+1) +4*l + t
			if nwl_strt>0 and div==0:
				print_rectangle('BilQua',p1,p2,p3,p4,ir)
			else:
				print_rectangle('CylRec',p1,p2,p3,p4,ir)
# at last, the transfinite rectangles
for div in range(theta_subdiv+1):
	for l in range(hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad):
		for t in range(1,5):
			ir += 1
			p1 = npcs*div +4*l + t
			p2 = npcs*div +4*l + t%4 +1
			p3 = npcs*div +4*(l+1) + t%4 +1
			p4 = npcs*div +4*(l+1) + t
			print_rectangle('TraQua',p1,p2,p3,p4,ir)

# prisms
nr_prism = 4*theta_subdiv
f.write("%-2i         NRPRISM\n" % nr_prism)
f.write("\n")
# 
ipri = 0
for div in range(theta_subdiv):
	for t in range(1,5):
		ipri += 1
		p1 = npcs*div + t
		p2 = ci[div]
		p3 = npcs*div + t%4 + 1
		p4 = npcs*(div+1) + t
		p5 = ci[div+1]
		p6 = npcs*(div+1) + t%4 + 1
		ndom = 1
		print_prism('TIprism',ndom,p1,p2,p3,p4,p5,p6,ipri)

# hexahedra
nr_hexa = 4*theta_subdiv*(hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad)
f.write("%-3i         NRHEXAS\n" % nr_hexa)
f.write("\n")
# 
ih = 0
for div in range(theta_subdiv):
	for l in range(hexlayers_core+hexlayers_inner_clad+hexlayers_outer_clad):
		for t in range(1,5):
			ih += 1
			p1 = npcs*div +4*l + t
			p2 = npcs*div +4*l + t%4 +1
			p3 = npcs*div +4*(l+1) + t%4 +1
			p4 = npcs*div +4*(l+1) + t
			p5 = npcs*(div+1) +4*l + t
			p6 = npcs*(div+1) +4*l + t%4 +1
			p7 = npcs*(div+1) +4*(l+1) + t%4 +1
			p8 = npcs*(div+1) +4*(l+1) + t
			nqi = nqtordiv*div + 4*l + t
			nqo = nqi + 4
			if l<hexlayers_core:
				ndom = 2
			elif l<hexlayers_core+hexlayers_inner_clad:
				ndom = 3
			else:
				ndom = 4
			if nwl_strt>0 and div==0:
				print_hexa('TraHex',ndom,p1,p2,p3,p4,p5,p6,p7,p8,ih)
			else:
				print_hexa('TorHex',ndom,p1,p2,p3,p4,p5,p6,p7,p8,ih,nqo,nqi)

# Tetrahedra
f.write("0 NRTETRA\n")
f.write("\n")
# Pyramids
f.write("0 NRPYRAM\n")
# Close geometry file
f.close()