#!/usr/bin/env python

"""
Perform a windowed median flattening of the images in an input fits file.
X and Y are medianed separately.

Stephen Bailey, LBL, 2013
"""

import numpy as N
from scipy.ndimage import median_filter
import pyfits

import optparse
parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--input", type="string",  help="input data")
parser.add_option("-o", "--output", type="string",  help="output data")
parser.add_option("-n", "--numimages", type="int",  help="number of images to process", default=1000)
parser.add_option("-w", "--window", type="int",  help="median window size (%default)", default=51)
opts, args = parser.parse_args()

nw = opts.window  #- Median filtering window size

fx = pyfits.open(opts.input)
if opts.numimages < len(fx):
    print "Processing %d of %d input images" % (opts.numimages, len(fx))

nmax = min(opts.numimages, len(fx))

#- start at 1 since HDU 0 is the mask
for i in range(1, nmax):
    print i
    d = fx[i].data
    ny, nx = d.shape
    
    #- First pass: xy median filtering
    dx = median_filter(d, size=(1,nw))  #- x filter
    ddx = d/dx
    ddx[ddx != ddx] = 0.0  #- reset NaN from 0.0/0.0
    dxy = median_filter(ddx, size=(nw,1))  #- y filter
    pxy = ddx/dxy

    #- Second pass: replace pixels which may have pulled the median, and repeat
    ibad = N.where(pxy < 0.5)
    d[ibad] = dx[ibad]

    dx = median_filter(d, size=(1,nw))
    ddx = d/dx
    ddx[ddx != ddx] = 0.0
    dxy = median_filter(ddx, size=(nw,1))
    pixflat = (ddx/dxy)

    #- Use original values for the masked regions
    pixflat[ibad] = pxy[ibad]
    
    #- Cleanup NaN and extreme values
    pixflat[pixflat != pixflat] = 0.0
    pixflat[pixflat < 0.0] = 0.0
    pixflat[pixflat > 10.0] = 0.0

    #- Special case: r1/r2 central rows
    if opts.input.count('r1') == 1 or opts.input.count('r2') == 1:
        print "Fixing central rows for R channel"
        yy = range(2050,2060) + range(2068,2078)
        ymid = range(2060, 2068)
        for ix in range(0, nx):
            p = N.polyfit(yy, d[yy, ix], 5)
            pix = N.polyval(p, ymid)
            pixflat[ymid, ix] = d[ymid, ix] / pix
     # Make hdu primary and hdulist
    if i == 1:
        pyfits.writeto(opts.output, pixflat,fx[i].header, clobber=True)
    else:
        pyfits.append(opts.output, pixflat,fx[i].header)
fx.close()
    

