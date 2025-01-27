# GNS-scripts
Pipeline for GNS 1
GNS Data Reduction Pipeline


Scripts

All scripts are located in XXX 

In principle, all the scripts are identical, but there may be small differences depending on the filter.
Therefore we keep the scripts associated to the filters. 
We may want to simplify this.

Parameters need to be set/checked in the scripts before running them, for example: band, NDIT, field number, paths to data, etc.

Many scripts will write control images to the tmp directory. Do not forget to check those.
Do not forget to clean up raw and intermediate data when you are finished with a field.

IDL usage

Log into teacup with a terminal.
Change to the directory containing the scripts.
Call IDL.
IDL> script_name (do not include the ‘.pro’ suffix when calling a script)

When you open additional terminals on the remote machine, call ‘xterm’ within the terminal that is already open on teacup. Do not log in externally again because this will be costly on IDL licenses.

Raw data

Location in 
/home/data/raw/2015/[Field-Nr]/

Within this directory there are sub-directories for each filter (J H Ks) and one for the darks (Dark)

Within these subdirectories there are directories Flat (flat field data) Sky (sky observations for each field) Field (science observations for each target)

The flat and dark data are located in sub-directories corresponding to the science pointings (fields).

To get raw data you need to download them from the ESO archive. You can find them with date/instrument/filter combination on the archive or, more easily, with the OB identification number.

The observing blocks are called FAST_GC-[Field-Nr]-[JHKs] and FAST_GC-[Field-Nr]-[JHKs]-Sky
for object and sky observations, respectively.

De-select the ACQUISITION observation.
Copy the provided download shell script (downloadRequestXXXscript.sh) into the raw directory on 
the server. change permissions: chmod u+x downloadRequestXXXscript.sh Then run the shell script: ./downloadRequestXXXscript.sh
Do this for the object and sky observations. 

You also need to download darks and flat fields.
Use the ESO gasgano tool to move the different types of files to their respective locations.
Decompress the data with uncompress.

Raw SKY data go into /home/data/raw/2015/H/Sky/[Field-Nr]
IMPORTANT: PLEASE SAVE list.txt from previous data release: cp list.txt lis_DR1.txt
Then uncompress sky observations and create the list of files (ls HAWKI*fits > list.txt).

1) Calibration: Flat and Dark


dark.pro, [field], [band], [date]

    Raw dark calibration files can be found on teatime in 
    /home/data/working/GNS_2/Dark/[field_nr]
  
    Create list of raw files in the same directory (ls HAWK* > list.txt)
    
    Set the correct field number, path and band in dark.pro

    dark.pro

    Choose correct field, filter and date for dark

2)  Create flat

    Download flat field data from ESO Science Archive (must be correct band, twilight flat).
    /home/data/working/GNS_2/H/Flat/[field_nr]
    Create list of raw files in this directory with ls *fits > list.txt (for documentation).

    gunzip the data, open gasgano tool and load flat field data into gasgano
    Run recipe  hawki_cal_flat
   
     joinflats.pro  —> create flat field , bad pixel map (bpm) and mask for continuing data reduction
     
     Each quadrant is individually normalised to a median value of 1.0.  

    Delete raw files in “Flat" directory. Do not delete list.txt in order to know which raw data    
    were used to create the flat.

1) Data reduction

3) sky.pro

    Create normalized sky image (“master sky”). The last frame in the sky cubes is typically a mean    
    image of the sub-integration and is not used. The first frame is often of lower quality in the 
    FastPhot mode and is therefore also discarded. Make sure to set the correct NDIT variable.
 
    The dark is subtracted form the sky. Stars (or most of their area) are filled in with random values 
    around the mean sky. Finally, each quadrant is normalised by its sky value. 

    In the end the master sky will look almost like the flat field. There are some clear differences, 
    though. They can be mostly seen in the form of vertical bright and dark stripes. Since the sky is  
    dark subtracted, the stripes cannot be caused by the dark. I suppose that we are dealing with 
    some bias that may depend on illumination or readout mode. Perhaps this could be remedied if     
    we were using flat field data acquired in FASTPHOT mode.

4) fullbpm,pro

    Create final bpm (bad pixel map) from master sky and dark. 
    I do not remember why we do not use the bpms already created from the flat and the dark in 
    steps 1 and 2, but we (or Paco) probably had good reasons for proceeding as we do. The 
    differences should be negligible in any way.

5) makemask_python.py

    Make a final mask from the bpm and from the preliminary mask (that masks the chip borders). It 
    will be used to mask regions on the chips that are too large to be interpolated (individual pixels 
    will be interpolated, but this is hard to do for clusters).

6) reduce.pro

   Create list of raw files 
   cd /home/data/GNS/[Bnad]/Field/[field_nr]
   ls HAWKI*fits.Z > list.txt

   The sky is fitted to each individual exposure because the sky values vary consistently with the  exposure series. Fitting the sky to the long exposures can also be done and will accelerate the code somewhat. I have also run a test (reduce_test.pro) on fitting the sky to each exposure and each quadrant separately, but this is complex and even slower. The difference between the quadrants is only of the order of 1%. So I have decided that this is not necessary.

Careful/To Do: Make sure that the right sky is selected. It may be that for some fields we have no good sky. We may even fare better with creating a single sky for all fields.
Sometimes, it is possible to use sky observations taken before AND after the object observations. 


7) cleancubes,pro
    This program cleans the exposures from the few hot/cold pixels that may remain in the data.
    This step is important because individual hot pixels can mess up the holography procedure.


8) makecubes.pro
    Create data cubes for each individual chip.
 
    After this step and before the next one you can delete the cubes in the ‘../cubes' directory:
    cd ../cubes/
    rm cube*fits.gz 

9) Dejitter: Align the individual exposures relative to the first one.
    
    findstars_dejitter.pro (needs extractpsf_fast.pro) + align_dejitter.py + dejitter.pro

    ASTROALIGN is used for this purpose.

    Call: 
    findstars_dejitter.pro  — to find stars in each exposure

    python align_dejitter.py [band] [field_nr] [n_exp]  — to find stars common between     
        exposure 1 (another one can be chosen as reference) and the other ones
        It is sufficient to run this for one chip (by default chip 1)
    
    dejitter.pro, field,  band, chip — use lists of common stars and correct for median offset

    The stars in chip 1 are used by default to determine relative alignment.

    Dejitter will remove the last exposure (the long exposure) from the cubes. 
    The long exposures will then be treated separately.

   NOTE: After dejitter, the astrometric information in the image headers will be removed.




11) longexposure.pro

      Create long exposure images for each chip.

12) findstars_lxp,pro, [field] , [chip]

      Find stars on long exposure image of a given field and chip.
     
       This star list is needed for astrometric alignment of the imaging from the different pointings 
      and filters with the VVV survey.

      Alternativey, if findstars_lxp.pro does not work well (it is automatic), you can use 
      findstars_lxp_manual.pro
      You need to select some (10-20) isolated, unsaturated reference stars to create the PSF.
    
13) Preparation of star list for VVV reference field

      This step has to be performed only once for a given field (not necessary for each filter)

      We use the VVV survey for aligning the images astrometrically with each other and with the 
      absolute astrometric reference frame.

     We do this in /home/data/VVV:

      We use the J-band because it has least saturation issues. Note that J-band images can look     
      quite different from K-band images in the GC (because of the extreme extinction).

      Here you first extract a small reference image for the region of interest from a larger VVV tile 
      (b333 or b334). Subsequently, the PSF is extracted and the stars are found on this VVV image.
   
      **prep_refim.pro**, [field] [band] (use  “J" for band -unless you have a very 
      good reason for a different choice) 
      y_size is set for a DIT of 1.26s and needs 
      to be adjusted for data with different windowing (e.g. y_size=512 for DIT = 0.85s)
      You may need the printed output of the size of the reference image to edit these parameters in 
      the next script (alignquadrants.pro). 

     python vvv_stars.py [field] [filter] [yy-mm-dd]
     The epoch and filter need to be provided so that the proper motions can be applied.


 **Outputs**: Image Field[number].fits.gz and star list Field[numer]_stars.txt amd Field[numer]_stars_fine.txt  in VVV directory. “fine”is then used by default. It takes into account the proper motions of the sources and therefore their change of position between the VVV epoch and the epoch of HAWK-I observation.

14) Now align the data with the VVV reference frame

     align_VVV.py +  alignquadrants.pro
 
  Align the HAWK-I data with the VVV data and with each other. The code uses the star lists 
  from steps 12 (HAWK-I) and 13 (VVV) for this purpose.  This code has to be run for each     
  chip.

align_VVV.py will use the VVV list without proper motion corrections.  This approach is necessary so that astroalign finds a solution easily. However, alignquadrants.pro will then identify the correct stars from the list with proper motion corrections.

NOTE:
  1) Please check the result of the alignment: Open the aligned long exposures    
      (lnx_jitter_X_aligned,fits.gz) of the four chips in tile mode in ds9 and then let them blink. 
      They should be precisely aligned.
  2) This program saves the parameters for second order polynomial image alignment in the
      IDL structures AlignmentPars_chipX in the data directory. These parameters will be needed 
      by the program alignframes.pro (see below).

3) Alignquadrants uses the auxiliary function gridselect.pro. It makes sure that the stars that are    
    used for alignment are distributed evenly across the field, which may avoid bias in 
    higher order transforms. Ngrid is currently set to 20 (corresponding to a 20 x 20 grid). The 
     optimal choice has not been tested.
 
alignquadrants will use the proper motion corrected stars from VVV.

  Use the lnx_jitter_[?]_VVV.fits.gz images to check alignment. Alignment must be precise,    
  otherwise we can run into problems much later (calibration, at the end of it all).

 NOTES: 
  (1) alignquadrants.pro creates a file named 'xy_off_f'+strn(field_nr)+'.SAV'
 This file will be used by alignframes.pro and also, much later, by calibrate.pro
 This fiel contains xy offsets of the final images that crop the field slightly, to avoid unnecessary   
  zeros at the edges. Any loss to or change of these files can lead to problems later (in 
  particular, in calibrate.pro).

(2) In alignquadrants you can chose whether to use proper motion-corrected VVV positions (default … ‘_stars_fine,txt' in line 50) or not. 

16) mosaic.pro 

     Create a large mosaic image to check the alignment of the chips.
     Look at the mosaics carefully to check whether everything is ok (images look ok, not elongated 
     or doubled stars).

    There will also be a mosaic at the VVV scale, to check and compare with VVV.


17) alignframes.pro 
      Uses the distortion solution computed by alignquadrants.pro to align every individual frame. 
      Normally you do not need to adjust the parameters in this script (except band and field).
      If you need to adjust some of them (such as the offset inside the aligned frame or the size of 
      the reference image), then you can find the corresponding values in alignquadrants.pro

   

18) subcubes.pro
      Create cubes of small fields (about 60” x 60”) on which the holography algorithm will be run.
      Call: subcubes, field_nr, chip_nr

Note on size of the sub-fields: Possibly, we could use larger sub-fields, which would save 
time. However, the prime reason for using smaller sub-fields is to compensate for differential tip-tilt motions in the field. We have never tested how important this differential tip-tilt is, but it appears to be clearly visible by eye. Do we need to use smaller fields for high precision PMs?


19) findstars_holo.pro

     Find stars in the aligned longexposure images in order to identify reference stars and 
     secondary stars for the holography algorithm.
     I have changed the threshold from [3.,3.] to [5.,5.]. I think this will not change the quality of the 
     holographic reconstruction significantly, but will save some computing time.

     Note: simple shift-and-add images look hardly any better than the longexposure images.
               Therefore, we use the simpler longepxosures to prepare holography.

   findstars_holo.pro works automatically (needs extractpsf.pro).

   There is an alternative version, findstars_holo_manual.pro. It allows the user to select the PSF    
    reference stars manually if deemed necessary.

20) runholo.pro
     
     Create sub-directories tmp1, tmp2, tmp3, and tmp4 in tmp/ and chip1, chip2, chip3, chip 4 in /
     cubeims before running the code.
     
                   Weighting?
                Test: No iteration in PSF extraction: No significant deterioration, speed hardly different.
                        I stick with iterative PSF extraction.
                I have eliminated a problem with the fitting of a constant background to the 
                PSF reference images

 Changes: Normalization of PSFs implemented correctly; PSF-submasks were not shifted correctly in previous verison; weighting of PSFs before superposition implemented

Holography code is holo_full.pro  in the IDL library folder MyRoutines. 

Edge effects (ripples) will be mostly taken care of by this code, but I have not tested whether there may be remnant effects. The latter will probably be minimised by cutting off the edges of the sub-images in the reconstruct.pro code. 

holo_full.pro implements Jackknife sampling. provides as output the holography images (and sub-images) for each cube as well as the exposure map, i.e. a map with the number of exposures contributing to each pixel (for the full image only). Noise maps are also provided.

nsub is the number of jackknife samples created. In case of the GLAO data we use nsub = 6.














