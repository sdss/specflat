How to generate lossy-fiber pixflats

Currently this takes 3 steps, primarily because the code was developed
in 3 steps with intermediate files to save debugging iteration time.
If we start doing this more than once or twice a year it could be combined
into a single script.

e.g. to generate pixflats using the data from MJD 56055 exps
143494-143553 (sp1) and MJD 56056 exps 143560-143619 (sp2):

1. Within IDL, run preproc_pixflats (in $SPECFLAT_DIR/pro/) to run
   sdssproc on multiple exposures and group them my exposure time in
   multi-HDU fits files:

outdir = getenv('BOSS_SPECTRO_REDUX')+'/test/sbailey/pixflats/2012a/'
preproc_pixflats, 56055, 'r1', 143494, 143553, outdir=outdir
preproc_pixflats, 56055, 'b1', 143494, 143553, outdir=outdir
preproc_pixflats, 56056, 'r2', 143560, 143619, outdir=outdir
preproc_pixflats, 56056, 'b2', 143560, 143619, outdir=outdir

That generates imgflat-b1-150.fits, imgflat-b1-500.fits, etc. for b2, r1, r2.
You can discard the imgflat-r?-500.fits files since they are saturated.


2. Median flatten those files

export DATADIR=$BOSS_SPECTRO_REDUX/test/sbailey/pixflats/2012a/
python $SPECFLAT_DIR/bin/median_flatten.py \
    -i $DATADIR/imgflat-b1-150.fits -o $DATADIR/medflat-b1-150.fits

etc. for other cameras and exposure times


3. Median stack those images, e.g.

export CAMERA=b1
export DATADIR=$BOSS_SPECTRO_REDUX/test/sbailey/pixflats/2012a/
python $SPECFLAT_DIR/bin/median_stack.py \
    -o $DATADIR/pixflat-$CAMERA.fits \
    $DATADIR/medflat-$CAMERA-*.fits

