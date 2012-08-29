#!/usr/bin/python

# extracts contact points and normals from an output log and writes them to
# stdout as VRML.  Contact points are output as spheres and normals are 
# output as arrows.  Contact visualization sizes and colors are selectable.

import sys
from optparse import OptionParser
import os

# structure for holding contact data
class ContactData:
  px = 0.0
  py = 0.0
  pz = 0.0
  nx = 0.0
  ny = 1.0
  nz = 0.0

# parse command line options
usage = "usage: extract-contacts [options] filename"
parser = OptionParser(usage=usage)
parser.add_option("-c", "--color", dest="color", help="The RGB color to write contact data", default="1 0 0")
parser.add_option("-s", "--scale", dest="scale", help="The scaling factor to apply to the contact data", type="float", default="1.0")
parser.add_option("-i", dest="iteration", help="The iteration from which to capture the contact data", type="int", default="-1")
(options, args) = parser.parse_args()

# get filename arguments
if len(args) != 1:
  parser.print_help()
  sys.exit(-1)

# make sure an iteration was specified
if options.iteration == -1:
  print 'No iteration specified!  You may need to use -i option.'

# get color
color = options.color
color_chk = color.split()
if len(color_chk) != 3:
  print 'Invalid color specified; using red'
  color = "1 0 0"

# get scale
scale = float(options.scale)

# determine contact points and normals
try:
  inp = open(args[0], 'r')
except IOError:
  print 'Unable to open ' + sys.argv[1] + ' for reading!'
  sys.exit(-2)

# read until file is exhausted
contacts = []
if options.iteration == -1:
  capturing = True
else:
  capturing = False
pstring = ""
istring = inp.readline()
while istring != "":
  # look for iteration string
  if options.iteration > -1 and istring.find('iteration:') != -1 and istring.find('simulation time:') != -1:
    # if capturing is already activated, we're done..  break out
    if capturing:
      break

    # get the iteration
    spl = istring.split()
    if len(spl) != 8:
      print 'Unexpected iteration string!'
      sys.exit(-3)
    iter = int(spl[4])

    # otherwise, see whether to activate capturing
    if iter == options.iteration:
      capturing = True

  elif istring.find('contact geometry:  (point) --> [') != -1:
    # see whether we capture
    if capturing == 0:
      pstring = istring
      istring = inp.readline()
      continue

    # replace ',', '[', and ']' in istring and pstring
    istring = istring.replace(',', ' ')
    istring = istring.replace('[', ' ')
    istring = istring.replace(']', ' ')
    pstring = pstring.replace(',', ' ')
    pstring = pstring.replace('[', ' ')
    pstring = pstring.replace(']', ' ')

    # create a new ContactData object
    cd = ContactData() 

    # get the normal from the previous string
    spl1 = pstring.split()
    cd.nx = float(spl1[1])
    cd.ny = float(spl1[2])
    cd.nz = float(spl1[3])
  
    # get the contact point
    spl2 = istring.split()
    cd.px = float(spl2[4])
    cd.py = float(spl2[5])
    cd.pz = float(spl2[6])

    # add this data to the list
    contacts.append(cd)

  # read the next line and save the current line
  pstring = istring
  istring = inp.readline()


# close the file
inp.close()

# print out a message if nothing found
if len(contacts) == 0:
  print 'No contact data found?!  Did you output logging data from Moby using the'
  print ' -l=1 option from driver?'
  sys.exit(-4)

# write VRML header
print '#VRML V2.0 utf8\n'

# write contact points and normals
for i in contacts:
  print 'Transform {'
  print '  translation ' + str(i.px) + ' ' + str(i.py) + ' ' + str(i.pz)
  print '  children'
  print '    Shape {'
  print '      appearance Appearance { material Material { diffuseColor ' + color + ' } }'
  print '      geometry Sphere { radius ' + str(scale) + ' } } }'
  print 'Shape {'
  print '  appearance Appearance { material Material { diffuseColor ' + color + ' } }'
  print '  geometry IndexedLineSet {'
  qx = str(i.px + i.nx*scale*10)
  qy = str(i.py + i.ny*scale*10)
  qz = str(i.pz + i.nz*scale*10)
  print '    coord Coordinate { point [ ' + str(i.px) + ' ' + str(i.py) + ' ' + str(i.pz) + ', ' + qx + ' ' + qy + ' ' + qz + '] }'
  print '    coordIndex [ 0 1 -1 ] } }\n' 

