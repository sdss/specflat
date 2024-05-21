#!/usr/bin/env python3

from astropy.io import fits
import argparse
import os.path as ptt
from sdss_access.path import Path
from glob import glob
from os import makedirs, chmod, getcwd, chdir, environ, remove
import subprocess
import pandas as pd

path = Path(release='sdsswork', preserve_envvars=True)

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



def LossyPixFlats(mjd, expid1, expid2, mjd2=None, expid2_1=None, expid2_2=None, outmjd=None, 
                  ver='', obs='apo', blue_exptime = '500', red_exptime = '150',indir=None):
    if outmjd is None: outmjd = mjd
    makedirs(ptt.join('.', 'pixflats', mjd+ver),exist_ok = True)
    indir2 = indir
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
        
        if mjd2 is not None:
            if indir2 is None:
                sdR2 = path.full('sdR', mjd=mjd2, br='b', id='*', frame='*')
                indir2 = ptt.dirname(ptt.dirname(sdR2).replace('apo',obs))
            for ccd in ccds:
                cmd = "preproc_pixflats, "+mjd2+", '"+ccd+"', "+str(expid2_1)+", "+str(expid2_2)+", outdir='$outdir'"
                if indir2 is not None:
                    cmd = cmd+", indir='"+indir2+"/'"
                t = cmdfile.write('idl -e "'+cmd+'"'+"\n")
            t = cmdfile.write("\n")

        if blue_exptime != red_exptime:
            t = cmdfile.write("rm "+ptt.join(outdir, 'imgflat-r?-'+str(blue_exptime)+'.fits')+"\n")
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
    parser.add_argument('--mjd2', required=False, default=None, type=str, help="Second MJD for Lossy Flat Data (if data taken over 2 days)")
    parser.add_argument('--expid2_1', required=False, default=None, help='First Lossy flat exposure of MJD2')
    parser.add_argument('--expid2_2', required=False, default=None,  help='Last Lossy flat exposure of MJD2')
    parser.add_argument('--gain_mjd', required=False, default=None, type=str, help='MJD of the gain exposure (defaults to mjd)')
    parser.add_argument('--gain_expid1', '-g', required=False, help='First Lossy flat exposure for BOSS Gain test (use CCD=R targeted fiber flat)')
    parser.add_argument('--APO', '--apo', action='store_true', help='Run for APO')
    parser.add_argument('--LCO', '--lco', action='store_true', help='Run for LCO')
    parser.add_argument('--blue_exptime', type=str, help='Blue Camera Flat Exposure Time', default=None)
    parser.add_argument('--red_exptime', type=str, help='Red Camera Flat Exposure Time', default=None)
    parser.add_argument('--ver', '-v', type=str, help='Run version for nights with multiple sets', default = '')
    parser.add_argument('--outdir', '-d', type=str, help='Output Working Directory', default='./')
    parser.add_argument('--indir', type=str, help='Manual Input Directory', default=None)
    args = parser.parse_args()

    cwd = getcwd()
    makedirs(args.outdir,exist_ok = True)
    chdir(args.outdir)

    if args.gain_expid1 is None: args.gain_expid1 = args.expid1
    if args.outmjd is None: args.outmjd = args.mjd
    if args.gain_mjd is None: args.gain_mjd = args.mjd

    print(pd.Series(vars(args)).to_string())
    print('-------------------')
    
    obs = []
    if args.APO: obs.append('apo')
    if args.LCO: obs.append('lco')

    for ob in obs:
        if ob == 'apo':
            blue_exptime = 500 if args.blue_exptime is None else args.blue_exptime
            red_exptime  = 150 if args.red_exptime  is None else args.red_exptime
        else:
            blue_exptime = 300 if args.blue_exptime is None else args.blue_exptime
            red_exptime  = 45  if args.red_exptime  is None else args.red_exptime

        LossyPixFlats(args.mjd, args.expid1, args.expid2, ver=args.ver, obs=ob, 
                      mjd2=args.mjd2, expid2_1 = args.expid2_1, expid2_2 = args.expid2_2,
                      outmjd=args.outmjd, blue_exptime = blue_exptime, 
                      red_exptime = red_exptime,indir=args.indir)

        findgain(args.gain_mjd, args.gain_expid1, ver=args.ver, obs=ob, indir=args.indir)

    chdir(cwd)
