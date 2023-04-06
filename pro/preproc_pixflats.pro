; NAME:
;   preproc_pixflats
;
; PURPOSE:
;   preprocess lossy-fiber pixel flat exposures with sdssproc and
;   group into multi-HDU files per camera.
;
; CALLING SEQUENCE:
;   preproc_pixflats, mjd, camera, expstart, expstop
;
; INPUTS:
;   mjd      - Modified Julian Day of data to pre-process
;   camera   - r1, r2, b1, or b2
;   expstart - first exposure to process
;   expstop  - last exposure to process
;
; OPTIONAL INPUTS:
;   outdir   - output directory; defaults to current directory
;
; OUTPUTS:
;   Creates files outdir/imgflat-CAMERA-EXPTIME.fit
;   where EXPTIME is the integer EXPTIME from the input FITS headers.
;
; PROCEDURES CALLED:
;   sdssproc()
;   mrdfits()
;
; TO DO:
;   Add override for $BOSS_SPECTRO_REDUX indir
;   Add override for prefix (default imgflat)
;   Add /clobber option
;
; REVISION HISTORY:
;   2013-05-21  Created by Stephen Bailey, LBL
;-------------------------------------------------------------------------

;- Fill in masked (ivar=0) regions with local median
;- This doesn't actually work so well, but leaving code here for reference
;- While investigating other possibilities.
function maskfilter, image, ivar
    d = 15  ;- Hard coded median filter half-width
    
    ; Create image with masked regions=NaN so that median will ignore them
    nanimg = image
    nanimg[where(ivar eq 0)] = !values.f_nan

    result = image
    nx = (size(image, /dim))[0]
    ny = (size(image, /dim))[1]
    for row=0,ny-1 do begin
        ;- If the entire row is masked, don't bother (e.g. at CCD edges)
        iz = where(ivar[*,row] eq 0)
        if n_elements(iz) lt nx then begin        
            for i=0,n_elements(iz)-1 do begin
                col = iz[i]
                ymin = max([0, row-d])
                ymax = min([row+d, ny-1])
                xmin = max([0, col-d])
                xmax = min([col+d, nx-1])
                ;- Median will ignore the NaNs
                result[col, row] = median(nanimg[xmin:xmax, ymin:ymax])
            endfor
        endif
    endfor
    
    ; Fix regions that were masked larger than the median filter
    ; Set them back to input image rather than leaving NaN values
    ii = where(finite(result, /nan) gt 0, nbad)
    if (nbad gt 0) then result[ii] = image[ii]
    
    return, result
end

pro preproc_pixflats, mjd, camera, expstart, expstop, outdir=outdir,indir=indir
   
    RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR


   if not keyword_set(outdir) then outdir = '.'
   outdir += '/'   ; for later outdir+outfile
   if not keyword_set(indir) then indir = getenv('BOSS_SPECTRO_DATA') + '/'
   index=0
   for expid=expstart, expstop do begin
       camexpid = camera + '-' + string(expid, format="(I08)")
       infile = strcompress(string(mjd), /remove_all)+'/sdR-'+camexpid+'.fit.gz'
       print,indir+infile
       if not file_test(indir + infile) then continue
       x = mrdfits(indir + infile, 0, hdr, /silent)
       exptime = uint(sxpar(hdr, 'EXPTIME'))
       
       outfile = 'imgflat-'+camera+'-'
       outfile += strcompress(string(exptime), /remove_all)
       outfile += '.fits'
       
       print, infile, exptime
       ;- Create infile list 
       sdssproc, indir+infile, image, ivar, /silent, /applybias, /nopixflat

       ;- Write mask image to HDU 0
       if not file_test(outdir+outfile) then $
         mwrfits, (ivar eq 0), outdir+outfile, /create
       ;- Make file list   
       sxaddpar, hdr_list,'NFILE', strcompress(index+1,/remove)
       sxaddpar, hdr_list, 'FILE'+strcompress(index,/remove), infile 
       ;;; print,hdr_list
       mwrfits,image,outdir+outfile,hdr_list,/silent             

       ;- Write out problematic region for study
       ;; mwrfits, image[2550:2700, 1150:1300], outdir+outfile, /silent


       ;; mkhdr, hdr, image, /image
       ;; sxaddpar, hdr, 'EXTNAME', 'IMG-'+strupcase(camexpid)       
       ;; mwrfits, image, outdir+outfile, hdr, /silent
       ;; sxaddpar, hdr, 'EXTNAME', 'IVAR-'+strupcase(camexpid)       
       ;; mwrfits, ivar, outdir+outfile, hdr, /silent       
   index += 1 
   endfor
   print, "wrote output to", outdir
end

