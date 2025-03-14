#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from tqdm.auto import tqdm
import os.path as ptt
import glob
import pandas as pd
import numpy as np
from os import makedirs, getenv

from shutil import rmtree
import tarfile
from datetime import date
import astropy.time
from IPython.display import display, HTML

try:
    from matplotlib import cm
except:
    pass
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from sdss_access import Path
path = Path(release='sdsswork', preserve_envvars=True)


def get_cmap():
    if True:#try:
        cmap16 = plt.colormaps['seismic']
        cmap16 = cmap16(np.linspace(0,1,13))

        viridis = plt.colormaps['viridis']
        newcolors = viridis(np.linspace(0, 1, 386))

        newcolors[0:109, :] = cmap16[0]  #32
        newcolors[109:147, :] = cmap16[1] #16
        newcolors[147:171, :] = cmap16[2] #8
        newcolors[171:183, :] = cmap16[3] #4
        newcolors[183:189, :] = cmap16[4] #2
        newcolors[189:193, :] = cmap16[5] #1
        newcolors[193:194, :] = cmap16[6] #0
        
        newcolors[194:198, :] = cmap16[7] #1
        newcolors[198:204, :] = cmap16[8] #2
        newcolors[204:216, :] = cmap16[9] #4
        newcolors[216:240, :] = cmap16[10] #8
        newcolors[240:278, :] = cmap16[11] #16
        newcolors[278:, :] = cmap16[12] #32

    else:#except:
        cmap16 = cm.get_cmap('seismic', 13)
        cmap16.colors=cmap16(range(0,13))

        viridis = cm.get_cmap('viridis', 253)
        newcolors = viridis(np.linspace(0, 1, 386))
   
    
        newcolors[0:109, :] = cmap16.colors[0]  #32
        newcolors[109:147, :] = cmap16.colors[1] #16
        newcolors[147:171, :] = cmap16.colors[2] #8
        newcolors[171:183, :] = cmap16.colors[3] #4
        newcolors[183:189, :] = cmap16.colors[4] #2
        newcolors[189:193, :] = cmap16.colors[5] #1
        newcolors[193:194, :] = cmap16.colors[6] #0

        newcolors[194:198, :] = cmap16.colors[7] #1
        newcolors[198:204, :] = cmap16.colors[8] #2
        newcolors[204:216, :] = cmap16.colors[9] #4
        newcolors[216:240, :] = cmap16.colors[10] #8
        newcolors[240:278, :] = cmap16.colors[11] #16
        newcolors[278:, :] = cmap16.colors[12] #32

    newcmp = ListedColormap(newcolors)
    return(newcmp)
get_cmap()
################################################################################################
################################################################################################

def run_analysis(savdir = 'bpm', OBS='APO', clobber=False, caltype='bpm', show=False, specflat_product=None, print_filelist=False, 
                 cmap='seismic', fcolors=['#ff7f0e','#2ca02c'], test=False, mjd=None, bpm_raw=False, term=False,mjd_only=None,
                 bias_range_red =None, bias_range_blue =None, bias_range =None, percentile=None,vrange=None,
                 bias_diff_range=None, bias_diff_range_blue=None, bias_diff_range_red=None, ccds=None, stage_dir = '..'):
    if caltype.lower()!='dark': make_struct(savdir, clobber=clobber)
    if ccds is None:
        if OBS.upper()=='APO': ccds=['r1', 'b1']
        else: ccds=['r2','b2']
 

    bias_range_bkup = bias_range
    bias_diff_range_bkup = bias_diff_range
    
    if term is False: 
        if type_of_script() != 'jupyter': term = True

    if mjd is not None: 
        if np.isscalar(mjd): mjd=[mjd]
        
    if specflat_product is not None: desc= 'SpecFlat '
    else: desc= ''
    for ccd in ccds:
        bias_range = bias_range_bkup
        bias_diff_range = bias_diff_range_bkup
        
        if ccd in ['r1', 'r2']:
            if bias_range_red is not None: 
                bias_range = bias_range_red
            if bias_diff_range_red is not None:
                bias_diff_range = bias_diff_range_red
        elif ccd in ['b1', 'b2']:
            if bias_range_blue is not None:
                bias_range = bias_range_blue
            if bias_diff_range_blue is not None:
                bias_diff_range = bias_diff_range_blue
        
        if   caltype.lower()=='bpm' : fpath1 = ptt.join(stage_dir, 'badPixelMask','?????','badpixels-*'+ccd+'*')
        elif caltype.lower()=='dark': fpath1 = ptt.join(stage_dir, 'badPixelMask','?????','badpixels-*'+ccd+'*')
        elif caltype.lower()=='flat': fpath1 = ptt.join(stage_dir, 'pixflats','?????','pixflatave-*'+ccd+'*')
        elif caltype.lower()=='bias': fpath1 = ptt.join(stage_dir, 'pixBias','?????','boss_pixbias*'+ccd+'*')
        if specflat_product is not None: 
            if   caltype.lower()=='bpm' : fpath2 = ptt.join(specflat_product, 'flats', 'badpixels-*'+ccd+'*.fits.gz')
            elif caltype.lower()=='dark': fpath2 = ptt.join(specflat_product, 'flats', 'badpixels-*'+ccd+'*.fits.gz')
            elif caltype.lower()=='flat': fpath2 = ptt.join(specflat_product, 'flats', 'pixflatave-*'+ccd+'*.fits.gz')
            elif caltype.lower()=='bias': fpath2 = ptt.join(specflat_product, 'biases', 'boss_pixbias-*'+ccd+'*.fits.gz')
        else: fpath2= None
        
        filt, files = find_product_files(fpath1, fpath2=fpath2)
        if filt is None or files is None: continue
        if print_filelist is True:
            for file in files: print(file)
        if caltype.lower()=='dark': disable=True
        else: disable=False
        if term: disable =True
        if mjd_only is not None:
            files1 = []
            for f in files:
                if int(getint2(f)) in mjd_only:
                    files1.append(f)
            files = files1
        for i, f in enumerate(tqdm(files, desc=desc+caltype+' '+ccd, disable=disable)):
            if disable:
                print(f'{i}/{len(files)}: {f}')
            if mjd is not None:
                if int(getint2(f)) not in mjd: continue
            if caltype.lower()=='bpm' and bpm_raw is True: bpm(f.replace(filt,ccd), cmap, fcolors, ccd = ccd, show=show, savdir=savdir)
            if i == 0 and caltype.lower()!='dark': continue
            if   caltype.lower()=='bpm' : bpm_diff( files[i-1].replace(filt, ccd), f.replace(filt,ccd), cmap, fcolors, ccd=ccd, show=show, savdir=savdir)
            elif caltype.lower()=='dark': get_darks(f.replace(filt, ccd), getint2(f), OBS, ccd=ccd)
            elif caltype.lower()=='flat': flat_diff(files[i-1].replace(filt, ccd), f.replace(filt,ccd), cmap, ccd=ccd, show=show, savdir=savdir, 
                                                    percentile=percentile, vrange=vrange)
            elif caltype.lower()=='bias': Bias_diff(files[i-1].replace(filt, ccd), f.replace(filt,ccd), cmap, ccd=ccd, show=show, savdir=savdir, 
                                                    bias_range=bias_range, bias_diff_range=bias_diff_range)
            if test is True and i == 1: break
    if caltype.lower()!='dark':  tar_figs(savdir, term=term)



    
################################################################################################
################################################################################################    
def get_darks(file1, mjd, obs, ccd = 'b1'):
    print('MJD:'+ mjd+ ' ccd:'+ ccd)
    if mjd == '58737': mjd='58736'
    f0 = fits.getheader(file1,0)['DARK STA']
    f1 = fits.getheader(file1,0)['DARK END']
    for i, frame in enumerate(range(int(f0),int(f1)+1)): 
        sdR = path.full('sdR', mjd=mjd, br=ccd[0], id=ccd[1], frame=str(frame).zfill(8))
        if not ptt.exists(sdR): sdR = sdR+'.gz'
        sdR = sdR.replace('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/boss/spectro/', 'https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/data/')
        print('    '+sdR)

    
################################################################################################
################################################################################################
def make_struct(folder, clobber=False):
    if clobber is True:
        if ptt.exists(folder): rmtree(folder)
    makedirs(ptt.join(folder, 'b1'), mode=0o777, exist_ok=True)
    makedirs(ptt.join(folder, 'b2'), mode=0o777, exist_ok=True)
    makedirs(ptt.join(folder, 'r1'), mode=0o777, exist_ok=True)
    makedirs(ptt.join(folder, 'r2'), mode=0o777, exist_ok=True)


    
################################################################################################
################################################################################################    
def type_of_script():
    try:
        ipy_str = str(type(get_ipython()))
        if 'zmqshell' in ipy_str:
            return 'jupyter'
        if 'terminal' in ipy_str:
            return 'ipython'
    except:
        return 'terminal'
################################################################################################
################################################################################################
def getint(name):
    num = ptt.basename(ptt.normpath(name))
    if '_' in num: num, ver = num.split('_')
    return int(num)

################################################################################################
################################################################################################
def find_files(fpath):
    files=(glob.glob(fpath))
    for i, f in enumerate(files): files[i]=ptt.dirname(f)
    files=np.unique(files).tolist()
    files.sort(key=getint)
    return(files)

################################################################################################
################################################################################################
def getint2(name):
    name = ptt.basename(name)
    return name.split('-')[1]

################################################################################################
################################################################################################
def find_product_files(fpath, fpath2=None):
    files=(glob.glob(fpath))
    if fpath2 is not None: 
        files2 = glob.glob(fpath2)
        filenames_list1 = {ptt.basename(file) for file in files}
        files2 = [file for file in files2 if (ptt.basename(file) not in filenames_list1) and (ptt.splitext(ptt.basename(file))[0] not in filenames_list1)]
        files.extend(files2)
    files=np.unique(files).tolist()
    if len(files) == 0: return(None, None)
    files.sort(key=getint2)
    filt= (ptt.basename(files[0]).split('-')[-1].split('.')[0])
    return filt,files

################################################################################################
################################################################################################
def tar_figs(folder, term=False):
    ftype = ptt.normpath(folder)
    mjd=(int(float(astropy.time.Time( str(date.today())).jd)-2400000.5))
    print('creating '+ftype+'_'+str(mjd)+'.tar.gz archive')
    out = tarfile.open(ftype+'_'+str(mjd)+'.tar.gz', mode='w')
    try:
        print('adding '+ftype)
        out.add(ftype, arcname= ptt.basename(ptt.normpath(folder)))
    finally:
        print('closing tar archive')
        out.close()
    ftype =ptt.basename(ptt.normpath(folder))
    if term is False: 
        display(HTML("Download: <a href='https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/test/sean/specflat/Analysis/"+ftype+"_"+str(mjd)+".tar.gz' download> "+ftype+"_"+str(mjd)+".tar.gz </a>"))
    else: 
        tqdm.write("https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/test/sean/specflat/Analysis/"+ftype+"_"+str(mjd)+".tar.gz")
    
################################################################################################
################################################################################################
def bpm_diff(file1, file2, cmap, fcolors, ccd = 'b1', show=True, savdir='bpm'):
    masks = {'hot': 1, 'saturated': 2, 'weak': 4, 'bad_col': 8, 'bad_warm_col': 16, 'hand_added': 32}
    if ptt.basename(file1) == ptt.basename(file2): return
    if not (ptt.exists(file1) and ptt.exists(file2)): return 
    if (ptt.basename(file1+'.gz') == ptt.basename(file2)) or (ptt.basename(file2+'.gz') == ptt.basename(file1)): return
    bpm1 = fits.getdata(file1, ext = 0)
    bpm2 = fits.getdata(file2, ext = 0)

    try: diff = bpm1-bpm2
    except: return 

    def finddiff(bpm, diff):
        m = np.zeros_like(np.ravel(bpm))
        m[diff] = 1
        m = np.reshape(m, bpm.shape)
        return(np.where(m == 1))
    
    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax = plt.subplots(3,2, figsize=(15, 25))
    ax = ax.ravel()

    
    mask1_bit=get_bits(bpm1)
    mask2_bit=get_bits(bpm2)
   
    
    for i, bit in enumerate(masks.keys()):
        
        amax = 32
        img = ax[i].imshow(diff.astype(int), cmap=cmap, vmax=amax, vmin=-1*amax)
        plt.colorbar(img, ax = ax[i])
        mask1 = np.where(mask1_bit[bit] == 1)
        mask2 = np.where(mask2_bit[bit] == 1)
        ax[i].plot(mask1[1],mask1[0], ls ='', marker='x', color=fcolors[0], label = 'mask1')
        ax[i].plot(mask2[1],mask2[0], ls ='', marker='+', color=fcolors[0], label = 'mask2')
        fmask1 = np.where(np.ravel(bpm1) == masks[bit])[0]
        fmask2 = np.where(np.ravel(bpm2) == masks[bit])[0]        
        
        mdiff = finddiff(bpm1, np.intersect1d(fmask1, fmask2))
        ax[i].plot(mdiff[1], mdiff[0],  ls ='', marker = 'D', color=fcolors[1],  mfc='None', label = 'and', alpha = 1)
        ax[i].title.set_text(bit+ '\n' +'MJD: '+str(getint2(file1))+' '+str(getint2(file2))+'\n'+'Nmask: '+str(len(fmask1))+ ' '+str(len(fmask2)) +'\n' +'%mask: '+str(round(len(fmask1)/len(np.ravel(bpm1))*100,4))+'%'+ ' '+str(round(len(fmask2)/len(np.ravel(bpm2))*100,4))+'%')
        ax[i].legend(loc='best')
    plt.suptitle(str(getint2(file1))+'-'+str(getint2(file2))+ ' '+ ccd )
    plt.tight_layout()
    plt.savefig(ptt.join( savdir,ccd, str(getint2(file1))+'-'+str(getint2(file2))+ '_'+ ccd+'.png' ))
    if show is True: plt.show()
    else: plt.close()
    
    bpm1=0
    bpm2=0

################################################################################################
################################################################################################
def get_bits(bpm):
    bpm_r = bpm.ravel()
    masks = {'hot': 1, 'saturated': 2, 'weak': 4, 'bad_col': 8, 'bad_warm_col': 16, 'hand_added': 32}
    invmasks = dict((v, k) for k, v in masks.items())
    mask_bit={}
    for bit_val in (sorted(masks.values(), reverse=True)):
        bit = invmasks[bit_val]
        mask_bit[bit] = np.zeros_like(bpm_r)
        mask_bit[bit][np.where(bpm_r >= masks[bit])] = 1
        bpm_r = bpm_r -  masks[bit] * mask_bit[bit]
        mask_bit[bit]=(mask_bit[bit]).reshape(*bpm.shape)
        
    return(mask_bit)

################################################################################################
################################################################################################
def bpm(file1, cmap, fcolors, ccd = 'b1', show=True, savdir='bpm'):
    masks = {'hot': 1, 'saturated': 2, 'weak': 4, 'bad_col': 8, 'bad_warm_col': 16, 'hand_added': 32}
    if not (ptt.exists(file1)): return 
    bpm1 = fits.getdata(file1, ext = 0)
    
    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax = plt.subplots(3,2, figsize=(15, 25))
    ax = ax.ravel()

    
    mask1_bit=get_bits(bpm1)
    
    for i, bit in enumerate(masks.keys()):
        amax = 32
        img = ax[i].imshow(bpm1, cmap=cmap, vmax=amax, vmin=-1*amax)
        plt.colorbar(img, ax = ax[i])
        mask1 = np.where(mask1_bit[bit] == 1)
        ax[i].plot(mask1[1],mask1[0], ls ='', marker='x', color=fcolors[1], label = bit)
        fmask1 = np.where(np.ravel(bpm1) == masks[bit])[0]
        
        ax[i].title.set_text(bit+ '\n' +str(len(fmask1)) +'\n' +str(round(len(fmask1)/len(np.ravel(bpm1))*100,4))+'%')
        ax[i].legend(loc='best')    
    plt.suptitle(str(getint2(file1))+' '+ ccd )
    plt.tight_layout()
    plt.savefig(ptt.join(  savdir,ccd, str(getint2(file1))+ '_'+ ccd+'.png' ))
    if show is True: plt.show()
    else: plt.close()
   
    
################################################################################################
################################################################################################
def flat_diff(file1, file2,cmap, ccd = 'b1', show=True, savdir='Flats', percentile=None, vrange=None):
    if ptt.basename(file1) == ptt.basename(file2): return
    if not (ptt.exists(file1) and ptt.exists(file2)): return 
    if (ptt.basename(file1+'.gz') == ptt.basename(file2)) or (ptt.basename(file2+'.gz') == ptt.basename(file1)): return
    bpm1 = fits.getdata(file1, ext = 0)
    bpm2 = fits.getdata(file2, ext = 0)

    try: diff = bpm1-bpm2
    except: return 

    if percentile is not None:
        if len(percentile) == 1:
            perc = [50-(percentile/2),50+(percentile/2)]
        else:
            perc = percentile
    
        vrange = [np.nanmin([np.percentile(bpm1, perc[0]),np.percentile(bpm2, perc[0])]),
                  np.nanmax([np.percentile(bpm1, perc[1]),np.percentile(bpm2, perc[1])])]
        amax = np.max([np.abs(np.nanpercentile(diff,perc[0])), np.nanpercentile(diff,perc[1])])
    elif vrange is not None:
        perc = [20,80]
        amax = np.max([np.abs(np.nanpercentile(diff,perc[0])), np.nanpercentile(diff,perc[1])])
    else:
        amax = np.max([np.abs(np.nanmin(diff)), np.nanmax(diff)])
        vrange = [None, None]

    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax = plt.subplots(2,2, figsize=(15, 15))
    ax = ax.ravel()
    i=0
    img = ax[i].imshow(bpm1, cmap=cmap, vmin = vrange[0], vmax = vrange[1],origin='lower' )
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(getint2(file1))

    i=1
    img = ax[i].imshow(bpm2, cmap=cmap, vmin = vrange[0], vmax = vrange[1],origin='lower' )
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(getint2(file2))
    
    i=2
    
    img = ax[i].imshow(diff, cmap=cmap, vmax=amax, vmin=-1*amax,origin='lower')
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(str(getint2(file1))+'-'+str(getint2(file2)))

    i=3
    histogram1, bin_edges1 = np.histogram(bpm1, bins=256)
    histogram2, bin_edges2 = np.histogram(bpm2, bins=256)
    
    ax[i].plot(bin_edges1[0:-1], histogram1,label=getint2(file1))
    ax[i].plot(bin_edges2[0:-1], histogram2,label=getint2(file2))
    ax[i].set_yscale('log')
    ax[i].legend(loc='best')    
    plt.suptitle(str(getint2(file1))+'-'+str(getint2(file2))+ ' '+ ccd )
    plt.tight_layout()
    plt.savefig(ptt.join( savdir,ccd, str(getint2(file1))+'-'+str(getint2(file2))+ '_'+ ccd+'.png' ))
    if show is True: plt.show()
    else: plt.close()
    

    
################################################################################################
################################################################################################
def Bias_diff(file1, file2, cmap, ccd = 'b1', show=True, savdir='bias', bias_range= None, bias_diff_range= None):
    if ptt.basename(file1) == ptt.basename(file2): return
    if not (ptt.exists(file1) and ptt.exists(file2)): return 
    if (ptt.basename(file1+'.gz') == ptt.basename(file2)) or (ptt.basename(file2+'.gz') == ptt.basename(file1)): return
    bpm1 = fits.getdata(file1, ext = 0)
    bpm2 = fits.getdata(file2, ext = 0)

    try: diff = bpm1-bpm2
    except: return 

    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax = plt.subplots(2,2, figsize=(15, 15))
    ax = ax.ravel()
    i=0
    if bias_range is not None:
        img = ax[i].imshow(bpm1, cmap=cmap, vmin=bias_range[0], vmax=bias_range[1])
    else:
        img = ax[i].imshow(bpm1, cmap=cmap)
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(getint2(file1))

    i=1
    if bias_range is not None:
        img = ax[i].imshow(bpm2, cmap=cmap, vmin=bias_range[0], vmax=bias_range[1])
    else:
        img = ax[i].imshow(bpm2, cmap=cmap)
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(getint2(file2))

    i=2
    amax = np.max([np.abs(np.nanmin(diff)), np.nanmax(diff)])
    if bias_diff_range is not None:
        amax = bias_diff_range[1]
    img = ax[i].imshow(diff, cmap=cmap, vmax=amax, vmin=-1*amax)
    plt.colorbar(img, ax = ax[i])
    ax[i].title.set_text(str(getint2(file1))+'-'+str(getint2(file2)))

    i=3
    if bias_range is not None:
        histogram1, bin_edges1 = np.histogram(bpm1, bins=256, range = bias_range)
        histogram2, bin_edges2 = np.histogram(bpm2, bins=256, range = bias_range)
    else:
        histogram1, bin_edges1 = np.histogram(bpm1, bins=256)
        histogram2, bin_edges2 = np.histogram(bpm2, bins=256)
    
    ax[i].plot(bin_edges1[0:-1], histogram1,label=file1.split('/')[-2])
    ax[i].plot(bin_edges2[0:-1], histogram2,label=file2.split('/')[-2])
    ax[i].set_yscale('log')
    ax[i].legend(loc='best')   
    plt.suptitle(str(getint2(file1))+'-'+str(getint2(file2))+ ' '+ ccd )
    plt.tight_layout()
    plt.savefig(ptt.join( savdir,ccd, str(getint2(file1))+'-'+str(getint2(file2))+ '_'+ ccd+'.png' ))
    if show is True: plt.show()
    else: plt.close()
    bpm1 = None
    bpm2 = None
    diff = None

    
################################################################################################
################################################################################################
def Gain_analysis(gainfiles, show=False, term=False, savdir = '.'):
    gains=pd.DataFrame()
    plt.close('all')
    gains = pd.concat([gains, pd.DataFrame.from_dict({'mjd': [58023], 'b2_0':[0.9787], 'b2_1':[1.0040], 'b2_2':[0.9647], 'b2_3':[1.0040], 
                                                  'r2_0':[1.5773], 'r2_1':[1.6527], 'r2_2':[1.5740], 'r2_3':[1.5840]})], ignore_index=True)
    

    for fi in gainfiles:
        with open(ptt.join(fi, 'boss_gain.log')) as f:
            ngain = pd.DataFrame({'mjd':[ptt.basename(fi)]})

            for line in f.readlines():
                cam = ['r2', 'b2']
                for c in cam:
                    if 'Camera = '+c in line: 
                        ngain[c+'_'+line[line.find('Amp = ')+6]] = [line[line.find('Gain =  ')+8:line.find(' from')]]
            if ngain.equals(pd.DataFrame.from_dict({'mjd':[ptt.basename(fi)]})): continue
        gains = pd.concat([gains,ngain], ignore_index=True)

    display(gains)

    
    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax=plt.subplots(2,1,figsize=(15, 15))

    ax[0].grid(False)
    if 'b2_0' in gains.columns:
        ax[0].plot(gains.mjd, gains.b2_0, label = 'Amp 0', marker='.')
        ax[0].plot(gains.mjd, gains.b2_1, label = 'Amp 1', marker='.')
        ax[0].plot(gains.mjd, gains.b2_2, label = 'Amp 2', marker='.')
        ax[0].plot(gains.mjd, gains.b2_3, label = 'Amp 3', marker='.')
        ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                               ncol=4, mode="expand", borderaxespad=0.)
    ax[0].set_xlabel('MJD')
    ax[0].set_ylabel('b2 gain')


    ax[1].grid(False)
    if 'r2_0' in gains.columns:
        ax[1].plot(gains.mjd, gains.r2_0, label = 'Amp 0', marker='.')
        ax[1].plot(gains.mjd, gains.r2_1, label = 'Amp 1', marker='.')
        ax[1].plot(gains.mjd, gains.r2_2, label = 'Amp 2', marker='.')
        ax[1].plot(gains.mjd, gains.r2_3, label = 'Amp 3', marker='.')
        ax[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                              ncol=4, mode="expand", borderaxespad=0.)
    ax[1].set_xlabel('MJD')
    ax[1].set_ylabel('r2 gain')

    plt.suptitle('BOSS-S CCD Gains')
    plt.tight_layout()

    plt.savefig(ptt.join(savedir,'gains_S.png'))
    if show is True: plt.show()
    else: plt.close()
    if term is False: display(HTML("Download: <a href='https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/test/sean/specflat/Analysis/gains_S.png' download> gains_S.png </a>"))
  
    gains_lco = gains

    gains=pd.DataFrame()
    gains = pd.concat([gains,pd.DataFrame.from_dict({'mjd': [58023], 'b1_0':[1.0447], 'b1_1':[1.0460], 'b1_2':[1.0140], 'b1_3':[1.0053], 
                                                  'r1_0':[2.046], 'r1_1': [1.6513], 'r1_2':[1.5913], 'r1_3':[1.5533]})], ignore_index=True)
    

    for fi in gainfiles:
        with open(ptt.join(fi, 'boss_gain.log')) as f:
            ngain = pd.DataFrame({'mjd':[ptt.basename(fi)]})

            for line in f.readlines():
                cam = ['r1', 'b1']
                for c in cam:
                    if 'Camera = '+c in line: 
                        ngain[c+'_'+line[line.find('Amp = ')+6]] = [line[line.find('Gain =  ')+8:line.find(' from')]]
            if ngain.equals(pd.DataFrame.from_dict({'mjd':[ptt.basename(fi)]})): continue
        gains = pd.concat([gains,ngain], ignore_index=True)

    display(gains)

    
    if show is False: plt.ioff()
    else: plt.ion()
    fig, ax=plt.subplots(2,1,figsize=(15, 15))

    ax[0].grid(False)
    if 'b1_0' in gains.columns:
        ax[0].plot(gains.mjd, gains.b1_0, label = 'Amp 0', marker='.')
        ax[0].plot(gains.mjd, gains.b1_1, label = 'Amp 1', marker='.')
        ax[0].plot(gains.mjd, gains.b1_2, label = 'Amp 2', marker='.')
        ax[0].plot(gains.mjd, gains.b1_3, label = 'Amp 3', marker='.')
        ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                               ncol=4, mode="expand", borderaxespad=0.)
    ax[0].set_xlabel('MJD')
    ax[0].set_ylabel('b2 gain')


    ax[1].grid(False)
    if 'r1_0' in gains.columns:
        ax[1].plot(gains.mjd, gains.r1_0, label = 'Amp 0', marker='.')
        ax[1].plot(gains.mjd, gains.r1_1, label = 'Amp 1', marker='.')
        ax[1].plot(gains.mjd, gains.r1_2, label = 'Amp 2', marker='.')
        ax[1].plot(gains.mjd, gains.r1_3, label = 'Amp 3', marker='.')
        ax[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                              ncol=4, mode="expand", borderaxespad=0.)
    ax[1].set_xlabel('MJD')
    ax[1].set_ylabel('r1 gain')

    plt.suptitle('BOSS-N CCD Gains')
    plt.tight_layout()

    plt.savefig(ptt.join(savedir,'gains_N.png'))
    if show is True: 
        plt.show()
        fig.canvas.draw()
    else: plt.close()
    if term is False: display(HTML("Download: <a href='https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/test/sean/specflat/Analysis/gains_N.png' download> gains_N.png </a>"))

    return(gains, gains_lco)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog=ptt.basename(sys.argv[0]),
            description='Analyze a library of specflat library')
    parser.add_argument('--staged_dir', help='Location of the staged specflat product products (default:  $SPECFLAT_WORK_DIR)', default=getenv('SPECFLAT_WORK_DIR'))
    parser.add_argument('--savedir', help='Location to save the Analysis figure (default: $SPECFLAT_WORK_DIR/Analysis)', default = None)
    parser.add_argument('--lco', '-l', action='store_true', required=False, help='Run for LCO data')
    parser.add_argument('--tagged_loc', '-t', help ='Location of the tagged specflat product for comparison')
    parser.add_argument('--tagged_loaded', action='store_true',help='Load the tagged specflat production location from "SPECFLAT_DIR"')
    parser.add_argument('--clobber', action='store_true', required=False, help='Clobber the existing figure')
    parser.add_argument('--bias', action='store_true',help='Analyze the Pixel Bias Frames')
    parser.add_argument('--bpm',  action='store_true',help='Analyze the Bad Pixel Masks')
    parser.add_argument('--gain', action='store_true',help='Analyze the CCD gains')
    parser.add_argument('--pixflat','--flat',  action='store_true', help='Analyze the Pixel Flats')
    parser.add_argument('--percentile', type=float, required=False, nargs=2, help='Percentile Cuts for Pixel Flat Analysis (ex: 45 80)', default=None)
    parser.add_argument('--vrange', type=float, required=False, nargs=2, help='scale range for Pixel Flat Analysis (ex: .992 1.008)', default=None)
    parser.add_argument('--mjd', required=False, nargs='+', default=None, help='List of MJDs to include')
    args = parser.parse_args()
    sav_base = ptt.join(getenv('SPECFLAT_WORK_DIR'), 'Analysis')
    if args.savedir is not None:
        sav_base = args.savedir
    if args.tagged_loaded:
        args.tagged_loc = getenv('SPECFLAT_DIR')
    if args.gain:
        gainfiles = find_files(ptt.join(args.staged_dir,'gains','*','boss_gain.log'))
        Gain_analysis(gainfiles, show=False, term=True, savdir=sav_base)

    if args.tagged_loc:
        flag = '_tagged'
    else:
        flag = ''

    obs = 'APO' if not args.lco else 'LCO'
    if args.bpm:   
        run_analysis(savdir=ptt.join(sav_base,'bhm'+flag), OBS=obs, caltype='bpm', specflat_product=args.tagged_loc, 
                    clobber=args.clobber, show=False, term=True, stage_dir = args.staged_dir)
    if args.bias:
        run_analysis(savdir=ptt.join(sav_base,'bias'+flag), OBS=obs, caltype='bias', specflat_product=args.tagged_loc, 
                    clobber=args.clobber, show=False, term=True, stage_dir = args.staged_dir)
    if args.pixflat:
        run_analysis(savdir=ptt.join(sav_base,'Flat'+flag), OBS=obs, caltype='flat', specflat_product=args.tagged_loc, 
                    clobber=args.clobber, show=False, term=True, stage_dir = args.staged_dir,percentile=args.percentile,
                    vrange = args.vrange, mjd_only = args.mjd)

            
