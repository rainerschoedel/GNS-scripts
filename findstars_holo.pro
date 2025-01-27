PRO FINDSTARS_HOLO, field, band, chip_nr

;field_nr = 1
;chip_nr = 1


data_path = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/data/'
tmp_path = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/tmp/'

; Find stars in images of each Chip
; that can serve for image alignment.


; General StarFinder settings
; ---------------------------
psf_size = 60.
maskrad = 20
threshold = [5.]
back_box = 0.
deblend = 0
min_correlation = 0.7
weigh_med = 0
unweighted = 1
upper_lev = 4.0e4
n_fwhm_back = 9.0
n_fwhm_fit = 2.0
n_fwhm_match = 1.0
mag_fac = 2L
n_width = 3.0
norm_max = 1
correl_mag = 2.0
niter = 2
compbg = 0
rel_thresh = 1
guide_x = ""
guide_y = ""

;loop chips
; ----------------

for i_chip = chip_nr, chip_nr do begin
  chip_nr = strn(i_chip)
  im = readfits(data_path + 'lnx_aligned_' + strn(i_chip)+ '.fits.gz')
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]
  noise = readfits(data_path + 'lnx_aligned_' + strn(i_chip)+ '_sig.fits.gz')

  
 ; create PSF
 ; --------------------------
  mmm, im, skymod, skysigma, skyskew
  psf_fwhm = 4.0
  mask = replicate(1,n1,n2) ; in principle we do not need to mask anything
  psf = extractpsf_fast(im, median(noise), mask, maskrad, psf_fwhm, tmp_path, UNWEIGHTED=1)
  psf = psf/total(psf)  ; normalization of PSF
  writefits, data_path + 'psf_chip' + chip_nr + '_holo.fits', psf
;   STOP


 ;    Extract PSF for field and re-run StarFinder
 ; --------------------------------------------

   threshold = [5.,5.]
   back_box = 0.
   deblend = 1
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
   forprint, TEXTOUT= data_path + 'stars' + '_' + chip_nr + '_holo.txt', /NOCOMMENT, x, y, f, c
; OPTIONAL OUTPUT FOR DEBUGGING
   writefits, data_path + 'im' + '_' + chip_nr + '.fits' , im
   writefits, data_path + 'resid' + '_' + chip_nr + '_holo.fits' , im - stars, /COMPRESS
   writefits, data_path + 'stars' + '_' + chip_nr + '_holo.fits' , stars, /COMPRESS
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
   map = image_model(x,y,f,n1,n2,'gaussian', dat)
   writefits, data_path + 'map' + '_' + chip_nr + '_holo.fits', map, /COMPRESS

endfor


END

