#!/usr/bin/python

import string
import re
import sys
import math
import os

# prints out the center-of-mass for a obj file 
# syntax: com <input>
if (len(sys.argv) < 3):
  print "syntax: transform <input filename> <output filename> <transform>"
  print "    ex: transform 'test.obj' 'out' '1 0 0 0; 0 1 0 0; 0 0 1 10; 0 0 0 1'"
  exit()

# open the OBJ file
infile = open(sys.argv[1],'r')

# open the output file
outfile = open(sys.argv[2], 'w')

# determine the transform
matrix = []
rows = sys.argv[3].split(';')
for i in rows:
  col = i.split()
  matrix.append([ float(col[0]), float(col[1]), float(col[2]), float(col[3]) ])
  
# transform all vertices
inp = infile.readline()
while (inp != ""):
	
	# look for a vertex
	if (re.compile("\s*v\s+.*").search(inp) != None):

		# found a vertex; split it
		vert = inp.split()

		# get the vertex
		v = [ float(vert[1]), float(vert[2]), float(vert[3]) ]

		# transform the vertex
		v1 = matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2]*v[2] + matrix[0][3]
		v2 = matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]*v[2] + matrix[1][3]
		v3 = matrix[2][0]*v[0] + matrix[2][1]*v[1] + matrix[2][2]*v[2] + matrix[2][3]

		# output the new vertex
		outfile.write("v ")
		outfile.write(str(v1))
		outfile.write(" ")
		outfile.write(str(v2))
		outfile.write(" ")
		outfile.write(str(v3))
		outfile.write("\n")
	else:
		outfile.write(inp)
	

	# read another line
	inp = infile.readline()

# close the file
infile.close()
outfile.close()

