PRO FINDSTARS_dejitter, field, band, chip_nr

; PURPOSE:
; Find stars in long exposures images of each cube
; that can serve for image alignment.


;band = 'Ks'
;field = '20'

basedir =  '/home/data/GNS/2015/' + band + '/' + field + '/'
cube_path = basedir +   'ims/'
out_path = basedir +   'data/'
tmp_path = basedir +  'tmp/'

chip_nr = strn(chip_nr)

n_exp = 49 ; total number of exposures normally 49, but  B6 - 44 at Ks




; General StarFinder settings
; ---------------------------
psf_size = 60.
maskrad = 20
threshold = [10.]
back_box = 60.
deblend = 0
deblost = 0
min_correlation = 0.7 ; 0.7 is default
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
compbg = 0
rel_thresh = 1
guide_x = ""
guide_y = ""

for i_e = 1, n_exp do begin

  cube = readfits(cube_path + 'chip' + chip_nr + '_cube' + strn(i_e) + '.fits.gz')
  sz = size(cube)
  n1 = sz[1]
  n2 = sz[2]
  if (sz[0] eq 3) then begin
     n3 = sz[3]
;     im = cube[512:1535,512:1535,n3-1]
     im = cube[*,*,n3-1]
  endif else begin
     n3 = 1                     ; some cubes may only contain a long exposure 
;     im = cube[512:1535,512:1535]
     im = cube
  endelse

  
 ; create PSF
 ; --------------------------
  mmm, im, skymod, skysigma, skyskew
  psf_fwhm = 4.0
  mask = replicate(1,n1,n2) ; in principle we do not need to mask anything
  psf = extractpsf_fast(im, skysigma, mask, maskrad, psf_fwhm, tmp_path, UNWEIGHTED=1)
  writefits, tmp_path + 'psf_'+ strn(i_e) + '.fits' , psf
 

 ; run StarFinder
 ; --------------------------------------------

  noise = im
  noise[*,*] = skysigma
  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', DEBLOST = deblost, $
	ESTIMATE_BG = compbg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
   forprint, TEXTOUT= out_path + 'dejitter_stars_' + chip_nr + '_' + strn(i_e) + '.txt', COMMENT='  x  y  f', x, y, f
   writefits, tmp_path + 'dejitter_stars_'  + chip_nr + '_' +strn(i_e) + '.fits', stars, /COMPRESS

   
endfor


END

