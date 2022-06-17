#!/usr/bin/env python

"""
Median combine multiple extensions of multiple fits files

This could be much more memory efficient by using fitsio to read
sub-images and medianize as we go.
"""

import sys
import pyfits
import numpy as N
from time import time

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-o", "--output",  type="string",  help="output file name")
parser.add_option("-n", "--numimages",  type="int",  help="number of images to include")

opts, args = parser.parse_args()

if opts.output is None:
    print "You must specify and output file (-o)"
    sys.exit(1)

if len(args) == 0:
    print "You must specify some input files"
    sys.exit(1)

#- First check how big these files are

nimages = dict()
image_shape = None
for filename in args:
    fx = pyfits.open(filename)
    nimages[filename] = len(fx)
    # Make the filelist in header  KZ
    if filename == args[0]:
        header_all = fx[len(fx)-1].header
        nfile_all = int(header_all["NFILE"])
    else:
        header_this = fx[len(fx)-1].header
        nfile_this = int(header_this["NFILE"])
        for file_loop in range(0,int(nfile_this)):
            temp=str(int(nfile_all) + file_loop)
            header_all.append(card="FILE{0}".format(temp))
            header_all["FILE{0}".format(temp)]=header_this['FILE{0}'.format(file_loop)]   
        nfile_all = nfile_all + nfile_this
        header_all["NFILE"]=nfile_all

    fx.close()
    
    if image_shape is None:
        d = pyfits.getdata(filename, 0)
        image_shape = d.shape

ny,nx = image_shape
if opts.numimages is not None:
    ntot = opts.numimages
else:
    ntot = N.sum(nimages.values())

t0 = time()
print 'Allocating memory'
data = N.zeros( (ntot, ny, nx), dtype='float32')
print '  --> %.1f' % (time() - t0, )

i = 0
for filename in args:
    if i >= ntot: break
    for j in range(nimages[filename]):                
        t0 = time()
        d = pyfits.getdata(filename, j).astype('float32')
        d[d != d] = 1  #- reset NaN
        data[i] = d
        print '%s %d %.1f' % (filename, i, time() - t0)
        i += 1
        if i >= ntot: break
        
print 'medianizing'
t0 = time()

#- medianize in chunks
median_image = N.zeros(image_shape)
step = 500
for xmin in range(0, nx, step):
    print xmin
    for ymin in range(0, ny, step):
        xmax = min(xmin+step, nx)
        ymax = min(ymin+step, ny)
        ix = slice(xmin, xmax)
        iy = slice(ymin, ymax)
        
        median_image[iy, ix] = N.median(data[:, iy, ix], axis=0)
        
print '%.1f' % (time() - t0, )

print 'writing output'
pyfits.writeto(opts.output, median_image, header_all,clobber=True)


