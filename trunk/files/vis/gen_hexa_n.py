import sys
import numpy as np

if len(sys.argv) < 2:
    print "usage: python gen_hexa_n.py [n]"
else:
    i = int(sys.argv[1])
    n = 2**i
    f = open("hexa_%d" % i, "w")
    xi = np.linspace(0, 1, n+1)
    vertices = [(xi[i], xi[j], xi[k]) for k in xrange(n+1) for j in xrange(n+1) for i in xrange(n+1)]
    def ivert(i,j,k):
        return i + (n+1)*j + (n+1)**2*k
    elems = [(ivert(i,  j,  k  ),
              ivert(i+1,j,  k  ),
              ivert(i+1,j+1,k  ),
              ivert(i,  j+1,k  ),
              ivert(i,  j,  k+1),
              ivert(i+1,j,  k+1),
              ivert(i+1,j+1,k+1),
              ivert(i,  j+1,k+1)) for k in xrange(n) for j in xrange(n) for i in xrange(n)]
    string = "%d" % len(vertices)
    for vertex in vertices:
        string = "".join([string, """
%f %f %f""" % vertex])
    string = "".join([string, """
%d""" % len(elems)])
    for elem in elems:
        string = "".join([string, """
12 %d %d %d %d %d %d %d %d""" % elem])
    f.write(string)
    f.close()
