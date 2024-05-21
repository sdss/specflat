#!/usr/bin/env python3

from astropy.io import fits
import argparse
import os.path as ptt
from sdss_access import Path
from glob import glob
from os import makedirs, chmod, getcwd, chdir, environ
import subprocess
import pandas as pd

path = Path(release='sdsswork', preserve_envvars=True)

def makepixBias(mjd, biasrange, obs='apo', ver=''):
    makedirs(ptt.join('.', 'pixBias', mjd+ver),exist_ok = True)
    with open(ptt.join('.', 'pixBias',mjd+ver,'pixBias_'+obs+'.cmd'), 'w') as cmdfile:
        t = cmdfile.write("#!/bin/bash"+"\n")
        
        sdR = path.full('sdR', mjd=mjd, br='b', id='*', frame='*').replace('apo',obs)+'*'
        cmd = "idl -e 'bolton_biasgen, ["
        fileList = glob(sdR)
        for i, sdR_f in enumerate(fileList):
            expid = int(ptt.splitext(ptt.splitext(ptt.basename(sdR_f))[0])[0].split('-')[-1])
            if ((expid < int(biasrange[0])) or (expid > int(biasrange[1]))): continue
            hdr=fits.getheader(sdR_f,0)
            hdr=fits.getheader(sdR_f,0)
            if hdr['FLAVOR'] == 'bias': cmd = cmd+'"'+sdR_f+'",'
        
        cmd = cmd[:-1] +"],/outbias' > pixBias_"+obs+".log"
        t = cmdfile.write(cmd+'\n')
        t = cmdfile.write('\n')

        sdR = path.full('sdR', mjd=mjd, br='r', id='*', frame='*').replace('apo',obs)+'*'
        cmd = "idl -e 'bolton_biasgen, ["
        fileList = glob(sdR)
        for i, sdR_f in enumerate(fileList):
            expid = int(ptt.splitext(ptt.splitext(ptt.basename(sdR_f))[0])[0].split('-')[-1])
            if ((expid < int(biasrange[0])) or (expid > int(biasrange[1]))): continue
            hdr=fits.getheader(sdR_f,0)
            if hdr['FLAVOR'] == 'bias': cmd = cmd+'"'+sdR_f+'",'
        cmd = cmd[:-1] +"],/outbias' >> pixBias_"+obs+".log"
        t = cmdfile.write(cmd+'\n')
        t = cmdfile.write('\n')
    cwd = getcwd()
    chdir(ptt.join('.','pixBias',mjd+ver))
    chmod(ptt.join('.','pixBias_'+obs+'.cmd'), 0o755)
    my_env = environ.copy()
    subprocess.call(ptt.join('.','pixBias_'+obs+'.cmd'), env=my_env, shell=True)
    chdir(cwd)


def makeBPM(mjd, darkrange, biasrange, obs='apo', ver=''):
    makedirs(ptt.join('.', 'badPixelMask', mjd+ver),exist_ok = True)
    sdR = path.full('sdR', mjd=mjd, br='b', id='*', frame='*')
    indir = ptt.dirname(sdR).replace('apo',obs)

    with open(ptt.join('.', 'badPixelMask', mjd+ver, 'bpm_'+obs+'.cmd'), 'w') as cmdfile:
            t = cmdfile.write("#!/bin/bash"+"\n")
            if obs.lower() == 'lco': ccds = ['r2','b2']
            else: ccds = ['r1','b1']

            for ccd in ccds:
                cmd = 'combinebpm,docams="'+ccd+'",nsig=10,ncount=25,darkstart='+str(darkrange[0])+',darkend='+str(darkrange[1])+',biasstart='+str(biasrange[0])+',biasend='+str(biasrange[1])+',/output, indir="'+indir+'"'
                cmd = "idl -e '"+cmd+"' > bpm_"+ccd+"_"+obs+".log"
                t = cmdfile.write(cmd+'\n')
                t = cmdfile.write('\n')
    cwd = getcwd()
    chdir(ptt.join('.','badPixelMask',mjd+ver))
    chmod(ptt.join('.','bpm_'+obs+'.cmd'), 0o755)
    my_env = environ.copy()
    subprocess.call(ptt.join('.','bpm_'+obs+'.cmd'), env=my_env, shell=True)
    chdir(cwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup and run Monthly Check scripts')
    parser.add_argument('--mjd', '-m', required=True, type=str,  help='mjd of monthly checkout')
    parser.add_argument('--biasrange', '-b', required=True, nargs=2, help="Start & End Bias Exposure Numbers")
    parser.add_argument('--darkrange', '-d', required=True, nargs=2, help="Start & End Dark Exposure Numbers")
    parser.add_argument('--APO', '--apo', action='store_true', help='Run for APO')
    parser.add_argument('--LCO', '--lco', action='store_true', help='Run for LCO')
    parser.add_argument('--type', required=False, nargs='*', help="build cal type default=['bpm','bias']", default=['bpm','bias'])
    parser.add_argument('--ver', '-v', type=str, help='Run version for nights with multiple sets', default = '')
    parser.add_argument('--outdir', '-o', type=str, help='Output Working Directory', default='./')
    args = parser.parse_args()

    cwd = getcwd()
    makedirs(args.outdir,exist_ok = True)
    chdir(args.outdir)

    print('-------------------')
    print(pd.Series(vars(args)).to_string())
    print('-------------------')
    args.type = [tp.lower() for tp in args.type]
    if args.APO is True: 
        if 'bpm' in args.type: makeBPM(args.mjd, args.darkrange, args.biasrange, obs='apo', ver=args.ver)
        if 'bias' in args.type: makepixBias(args.mjd, args.biasrange, obs='apo', ver=args.ver)
    if args.LCO is True: 
        if 'bpm' in args.type: makeBPM(args.mjd, args.darkrange, args.biasrange, obs='lco', ver=args.ver)
        if 'bias' in args.type: makepixBias(args.mjd, args.biasrange, obs='lco', ver=args.ver)

    chdir(cwd)
