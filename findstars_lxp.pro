PRO FINDSTARS_LXP, field, band, chip_nr


  basedir = '/home/data/GNS/2015/' + band + '/'
data_path = basedir + strn(field) + '/ims/'
out_path =basedir + strn(field) +  '/data/'
tmp_path = basedir + strn(field) +  '/tmp/'


; Find stars in images of each Chip
; that can serve for image alignment.


; General StarFinder settings
; ---------------------------
psf_size = 60.
maskrad = 20
threshold = [10.]
back_box = maskrad
deblend = 0
min_correlation = 0.7
weigh_med = 0
unweighted = 1
upper_lev = 3.0e4
n_fwhm_back = 9.0
n_fwhm_fit = 2.0
n_fwhm_match = 1.0
mag_fac = 2L
n_width = 3.0
norm_max = 1
correl_mag = 2.0
niter = 2
compbg = 1
rel_thresh = 1
guide_x = ""
guide_y = ""

;loop chips
; ----------------

for i_chip = chip_nr, chip_nr do begin
  chip_nr = strn(i_chip)

  im = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '.fits.gz')
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]
  noise = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '_sig.fits.gz')

 ; extract PSF
 ; --------------------------

   psf_fwhm = 4.0
   mask = replicate(1,n1,n2) ; in principle we do not need to mask anything
   psf = extractpsf(im, noise, mask, maskrad, psf_fwhm, tmp_path, /UNWEIGHTED, SF_THRESH = 10.)
   writefits, data_path + 'psf_lnx_'+ strn(i_chip) + '.fits' , psf
   print, 'Chip ' + chip_nr + '.txt, PSF FWHM: ' + strn(psf_fwhm*0.106) + ' "'

   
 ; Run StarFinder
 ; --------------------------------------------

   starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
   forprint, TEXTOUT= out_path + 'stars' + '_' + chip_nr + '.txt', x, y, f, COMMENT = '  x  y  f'
   writefits, out_path + 'found_stars' + '_' + chip_nr + '.fits', stars, /COMPRESS
   
endfor

END

