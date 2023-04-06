#!/usr/bin/env python3

from astropy.io import fits
import argparse
import os.path as ptt
from sdss_access import SDSSPath
from glob import glob
from os import makedirs, chmod, getcwd, chdir, environ, remove
import subprocess
import pandas as pd

path = SDSSPath(release='sdss5', preserve_envvars=True)

def run(dir_, cmd, log=None, indir = None):
    cwd = getcwd()
    chdir(ptt.join(dir_))
    chmod(ptt.join('.',cmd), 0o755)
    my_env = environ.copy()
    if indir is not None:
        my_env["BOSS_SPECTRO_DATA"] = indir
    if log is not None: 
        with open(ptt.join('.',log), 'w') as f:
            subprocess.call(ptt.join('.',cmd), env=my_env, shell=True, stdout = f, stderr = subprocess.STDOUT)
    else: subprocess.call(ptt.join('.',cmd), env=my_env, shell=True)
    chdir(cwd)


def findgain(mjd, expid1, ver='', obs='apo', indir=None):
    makedirs(ptt.join('.', 'gains', mjd+ver),exist_ok = True)
    if indir is None:
        sdR = path.full('sdR', mjd=mjd, br='b', id='*', frame='*')
        indir = ptt.dirname(sdR).replace('apo',obs)
    else:
        indir = ptt.join(indir, mjd)
    with open(ptt.join('.', 'gains', mjd+ver, 'boss_gain_'+obs+'.cmd'), 'w') as cmdfile:
        t = cmdfile.write("#!/bin/bash"+"\n")
        if obs == 'apo': cams = "['b1','r1']"
        else: cams = "['b2','r2']"
        cmd = "boss_gain, "+expid1+", docams="+cams+", indir='"+indir+"'"
        t = cmdfile.write('idl -e "'+cmd+'"'+' > boss_gain.log'+"\n")
    run(ptt.join('.', 'gains', mjd+ver), 'boss_gain_'+obs+'.cmd')



def LossyPixFlats(mjd, expid1, expid2, outmjd=None, ver='', obs='apo', blue_exptime = '500', red_exptime = '150',indir=None):
    if outmjd is None: outmjd = mjd
    makedirs(ptt.join('.', 'pixflats', mjd+ver),exist_ok = True)

    if indir is  None:
        sdR = path.full('sdR', mjd=mjd, br='b', id='*', frame='*')
        indir = ptt.dirname(ptt.dirname(sdR).replace('apo',obs))
    with open(ptt.join('.', 'pixflats', mjd+ver, 'lossy_fiber_'+obs+'.cmd'), 'w') as cmdfile:
        t = cmdfile.write("#!/bin/bash"+"\n")
        t = cmdfile.write("export BOSS_SPECTRO_DATA='"+indir+"'"+"\n")
        cwd = getcwd()
        outdir = ptt.join(cwd, 'pixflats', mjd+ver)
        t = cmdfile.write("export outdir='"+outdir+"/'"+"\n")

        if obs == 'apo': ccds = ['r1', 'b1']
        else: ccds = ['r2', 'b2']

        for ccd in ccds:
            cmd = "preproc_pixflats, "+mjd+", '"+ccd+"', "+str(expid1)+", "+str(expid2)+", outdir='$outdir'"
            if indir is not None:
                cmd = cmd+", indir='"+indir+"/'"
            t = cmdfile.write('idl -e "'+cmd+'"'+"\n")
        t = cmdfile.write("\n")

        if blue_exptime != red_exptime:
            t = cmdfile.write("rm "+ptt.join(outdir, 'imgflat-r?-'+blue_exptime+'.fits')+"\n")
            t = cmdfile.write("\n")

        t = cmdfile.write('for file in $outdir/imgflat*.fits; do'+"\n")
        t = cmdfile.write('    python $SPECFLAT_DIR/bin/median_flatten.py -i $file -o "${file/imgflat/medflat}"'+"\n")
        t = cmdfile.write('done'+"\n")
        t = cmdfile.write("\n")
    
        for ccd in ccds:
            t = cmdfile.write("if ls $outdir/medflat-"+ccd+"-*.fits 1> /dev/null 2>&1; then"+"\n")
            t = cmdfile.write("    python $SPECFLAT_DIR/bin/median_stack.py -o $outdir/pixflat-"+mjd+"-"+ccd+".fits $outdir/medflat-"+ccd+"-*.fits"+"\n")
            t = cmdfile.write("fi"+"\n")
        t = cmdfile.write("\n")

        if obs == 'apo': cams = "['b1','r1']"
        else: cams = "['b2','r2']"
        cmd = "spflatave, mjd="+mjd+", mjout="+outmjd+", docam="+cams+", outdir ='./', indir = './'"
        t = cmdfile.write('idl -e "'+cmd+'"')

    run(ptt.join('.', 'pixflats', mjd+ver), 'lossy_fiber_'+obs+'.cmd', log = "lossy_fiber_"+obs+".log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup and run Lossy flat tests scripts')
    parser.add_argument('--mjd', '-m', required=True, type=str,  help='mjd of Lossy Flat Data')
    parser.add_argument('--outmjd', '-o', required=False, type=str,  help='mjd of to apply Lossy Flat Data')
    parser.add_argument('--expid1', '-f', required=True, help='First Lossy flat exposure')
    parser.add_argument('--expid2', '-l', required=True, help='Last Lossy flat exposure')
    parser.add_argument('--gain_expid1', '-g', required=False, help='First Lossy flat exposure for BOSS Gain test (use CCD=R targeted fiber flat)')
    parser.add_argument('--APO', '--apo', action='store_true', help='Run for APO')
    parser.add_argument('--LCO', '--lco', action='store_true', help='Run for LCO')
    parser.add_argument('--blue_exptime', type=str, help='Blue Camera Flat Exposure Time', default='500')
    parser.add_argument('--red_exptime', type=str, help='Red Camera Flat Exposure Time', default='150')
    parser.add_argument('--ver', '-v', type=str, help='Run version for nights with multiple sets', default = '')
    parser.add_argument('--outdir', '-d', type=str, help='Output Working Directory', default='./')
    parser.add_argument('--indir', type=str, help='Manual Input Directory', default=None)
    args = parser.parse_args()

    cwd = getcwd()
    makedirs(args.outdir,exist_ok = True)
    chdir(args.outdir)

    if args.gain_expid1 is None: args.gain_expid1 = args.expid1
    if args.outmjd is None: args.outmjd = args.mjd

    print(pd.Series(vars(args)).to_string())
    print('-------------------')

    if args.APO is True: 
        obs='apo'
        LossyPixFlats(args.mjd, args.expid1, args.expid2, ver=args.ver, obs=obs, outmjd=args.outmjd, blue_exptime = args.blue_exptime, red_exptime = args.red_exptime,indir=args.indir)
        findgain(args.mjd, args.gain_expid1, ver=args.ver, obs=obs, indir=args.indir)
    if args.LCO is True: 
        obs='lco'
        LossyPixFlats(args.mjd, args.expid1, args.expid2, ver=args.ver, obs=obs, outmjd=args.outmjd, blue_exptime = args.blue_exptime, red_exptime = args.red_exptime,indir=args.indir)
        findgain(args.mjd, args.gain_expid1, ver=args.ver, obs=obs, indir=args.indir)

    chdir(cwd)
