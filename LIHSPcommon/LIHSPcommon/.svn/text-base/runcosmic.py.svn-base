# The following two lines are only needed as cosmic.py is not in this directory nor in the python path.
# They would not be required if you copy cosmics.py in this directory.
import sys
sys.path.append("../.") # The directory that contains cosmic.py

if (len(sys.argv)<2 or sys.argv[1].strip() == '-h' or sys.argv[1].strip() == 'h'):
    print """usage:
    file input path, file input image, begin img number, end img number.
    names must be in the form 'j0923_4sec_gain50_1_1294312143.fits[dddd].fits'
    and you must omit the [dddd].fits ending
    """
    exit(0)

import cosmics
import os
# Read the FITS :
path = sys.argv[1]
nameroot = sys.argv[2]
print path,nameroot

cleandir = path+'/unclean/'
#str = 'rm -rf %s' %cleandir
#os.system(str)
str = 'mkdir %s' %cleandir
os.system(str)

inname = '%s/%s' %(path,nameroot)
print inname
if os.path.exists(inname) == 0:
    print "no image"
    sys.exit()
array, header = cosmics.fromfits(inname)
# array is a 2D numpy array

# Build the object :
c = cosmics.cosmicsimage(array, gain=200.0,readnoise=30000.0, sigclip = 13.0, sigfrac = 1.0, objlim = 5.0)
# There are other options, check the manual...

# Run the full artillery :
c.run(maxiter = 4)

# Write the cleaned image into a new FITS file, conserving the original header :
strg = 'mv %s %s'%(inname, cleandir)
os.system(strg)	
fname = '%s/%s' %(path,nameroot)
cosmics.tofits(fname, c.cleanarray, header)

# If you want the mask, here it is :
cosmics.tofits("mask.fits", c.mask, header)
# (c.mask is a boolean numpy array, that gets converted here to an integer array)
