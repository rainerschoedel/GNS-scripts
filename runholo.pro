PRO RUNHOLO, field, band, chip

; THESE PARAMETERS CAN BE EDITED
; USUALLY, THEY STAY UNCHANGED

debug = 0

  
iter = ''  ; if '', then results of SF on long exposure is used - default
           ; if '2', then  results of SF on holo rebin 1 is used

case band of

   'Ks':  BEGIN
      ZP = 25.5 ; rough mean ZP for all HAWK-I detectors
                ; estimated from ESO QC data base
                ; USE ZP FOR CORRECT FILTER!
      mag_min = 15  ; faintest mag for reference stars
      mag_max = 11  ; brightest mag for reference stars
    END
   'H':  BEGIN
      ZP = 26.3                 ; rough mean ZP for all HAWK-I detectors
                ; estimated from ESO QC data base
                ; USE ZP FOR CORRECT FILTER!
      mag_min = 17  ; faintest mag for reference stars
      mag_max = 11  ; brightest mag for reference stars
    END
   'J':  BEGIN
      ZP = 26.3                 ; rough mean ZP for all HAWK-I detectors
                ; estimated from ESO QC data base
                ; USE ZP FOR CORRECT FILTER!
      mag_min = 18  ; faintest mag for reference stars
      mag_max = 11  ; brightest mag for reference stars
    END
endcase
   
DIT = 1.26                      ; DIT [s] of observat1ions

psf_path = '/home/data/GNS/2015/'+band+'/' + field + '/data/'
in_path = '/home/data/GNS/2015/'+band+'/' + field + '/subcubes/'
out_path = '/home/data/GNS/2015/'+band+'/' + field + '/cubeims/'
outpsf_path = '/home/data/GNS/2015/'+band+'/' + field + '/psfs/'
data_path = '/home/data/GNS/2015/'+band+'/' + field + '/data/'
tmp_path = '/home/data/GNS/2015/'+band+'/' + field + '/tmp/'
tmp_path_s = tmp_path +  'tmp'+strn(chip)+'/'

x_cube = 600 ; xaxis length of sub-cube
y_cube = 600 ; yaxis length of sub-cube
x_large = 2400  ; xaxis length of large cube
y_large = 1200  ; yaxis length of large cube
      
SIGDEV = 3.0 ; standard deviations used in RESISTANT_MEAN

; Use LXP PSF as a reference for PSF size and seeing
 ; PSF should be roughly the same for all chips
 ; use PSF from long exposure, NOT holography (NOT holo2)
lnx_psf = readfits(psf_path + 'psf_chip'+strn(chip)+'_holo.fits')
psf_fwhm = fwhm(lnx_psf)

; Parameters for reference star selection
; ---------------------------------------
delta_r = 3.*psf_fwhm
delta_mag = 2.5 ; 2.5 mag correspond to a factor of 10
n_ref_max = 20  ; max number of reference stars to use in second iteration

; Parameters for holography
; --------------------------

rebfac = 2
nsub = 10         ; number of jackknife images (HOLO_MOSAIC.PRO)
nsigma = [2.,2.] ; noise thesholds for PSF 
satlevel = 4.0e4 ; saturation threshold of data
unweighted = 0
psfnoise = 1
smoothmask = 0 ; to suppress edge effects from mask
holo_iter = n_elements(nsigma) - 1
normrad = fwhm(lnx_psf)
starlist = tmp_path_s + 'stars.txt'
psf_size = rebfac * (x_cube < y_cube)
airy = psf_gaussian(NPIXEL=psf_size,FWHM=2.0*rebfac,/NORMALIZE,/DOUBLE)
bord = 0  ; width of cosine shaped transition region from 0 to 3 sigma noise threshold
out_iter = 10    ; output intermediate result after processing out_iter images
psfout = 0       ; psf output?
holoout = 0
subpix = 1       ; sub-pixel alignment of PSFs, Y/N
psfavg = 0      ; set > 0 for PSF mean superposition, else median
clip = 0
minsupp = 0.5
maxsupp = 1.5
weightframes = 0

maskrad = round(4*psf_fwhm)
boxhw = round(6*psf_fwhm)   ; half width of size of PSF image for holography
psf_border = 1*psf_fwhm
;maskrad = round(8*psf_fwhm)
;boxhw = round(10*psf_fwhm)   ; half width of size of PSF image for holography
;psf_border = 0*psf_fwhm

;print, maskrad, boxhw
;determine number of pixels in circular aperture within maskrad/2
; needed for reference star selection
dummy = fltarr(2*boxhw+1,2*boxhw+1)
dummy[*,*] = 1
circmask = circ_mask(dummy,boxhw,boxhw,maskrad)
n_inner = total(circmask)

; large circular mask to avoid reference sources near saturated stars
dummy = fltarr(4*boxhw+1,4*boxhw+1)
dummy[*,*] = 1
circbig = circ_mask(dummy,2*boxhw,2*boxhw,2*maskrad)


x_sub_shift = x_cube/2
y_sub_shift = y_cube/2
nx = x_large/x_sub_shift - 1
ny = y_large/y_sub_shift - 1


 chip_nr = strn(chip)
 ; read positions and fluxes of stars detected on SSA image of this chip
 readcol, data_path + 'stars' + '_' + chip_nr + '_holo' + iter + '.txt', x_chip, y_chip, f_chip, correl
; for i_x = 0, 2 do begin
for i_x = 0, nx -1 do begin    ;(original)
  for i_y = 0, ny -1 do begin

   ; Compute longexposure image (for debugging) and support region
   filenam = '_' + strn(i_x) + '_' + strn(i_y)
   readcol, in_path  + 'chip' + chip_nr + '/masklist' + filenam + '.txt', mnames, FORMAT='A'
   readcol, in_path  + 'chip' + chip_nr + '/list' + filenam + '.txt', cnames, FORMAT='A'
   lxp = fltarr(x_cube,y_cube)
   support = fltarr(x_cube,y_cube)
   nm = n_elements(mnames)
   for j = 0, nm-1 do begin
    maskcube = readfits(in_path + 'chip' +  chip_nr + '/' + mnames[j])
    support = support + total(maskcube,3)
    cube = readfits(in_path + 'chip' +  chip_nr + '/' + cnames[j])
    lxp = lxp + total(cube,3)
   endfor
   accept = where(support gt 0)
   lxp[accept] = lxp[accept]/support[accept]
   support = support/max(support)
   writefits, tmp_path_s + 'lxp.fits', lxp
   writefits, tmp_path_s + 'support.fits', support

   ; Select stars in this sub-field
   xlo = i_x * x_sub_shift
   ylo = i_y * y_sub_shift   
   x_sub = x_chip - xlo
   y_sub = y_chip - ylo
   subind = where(x_sub ge 0 and x_sub le (x_cube-1) and y_sub ge 0 and y_sub le (y_cube-1), n_sub) 
   x_sub = x_sub[subind]
   y_sub = y_sub[subind]
   f_sub = f_chip[subind]
   ; save list of stars to pass it on to holo_mosaic
   forprint, TEXTOUT=starlist, /NOCOMMENT, x_sub, y_sub, f_sub
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
   substars = image_model(x_sub,y_sub,f_sub,x_cube,y_cube,'gaussian', dat)
   writefits, tmp_path_s + 'sources.fits', substars

  ; Select PSF  reference stars among stars in selected brightness range
   ; Use a triple criterion: 
   ; (1) Minimum and maximum brightness.
   ; (2) No brighter star within delta_r.
   ; (3) Any star within delta_r  must be at least delta_mag fainter.
 ; ------------------------------------------------------------------- 

   m_sub = ZP - 2.5 * alog10(f_sub/DIT)

   ord = sort(m_sub)
   x_sub = x_sub[ord]
   y_sub = y_sub[ord]
   m_sub = m_sub[ord]

   ; select by magnitude 
   ref_ind = where(m_sub le mag_min and m_sub ge mag_max,n_ref) 
   x_psf = x_sub[ref_ind]
   y_psf = y_sub[ref_ind]
   m_psf = m_sub[ref_ind]

  ; PSF stars must be isolated
   isolated_stars, x_sub, y_sub, m_sub, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
   x_psf =x_psf[ind_iso]
   y_psf =y_psf[ind_iso]
   m_psf =m_psf[ind_iso]
   isolated_stars, x_sub, y_sub, m_sub, x_psf, y_psf, m_psf, 0, maskrad, ind_iso
   x_psf =x_psf[ind_iso]
   y_psf =y_psf[ind_iso]
   m_psf =m_psf[ind_iso]
   n_ref = n_elements(m_psf)
   print, 'Found '+ strn(n_ref) + ' isolated reference stars.'

  ; exclude potentially saturated reference stars and those close to saturated sources
   sat_pixels = where(lxp gt satlevel, complement=not_saturated,n_saturated)
   if (n_saturated gt 0) then begin
     sat_mask = lxp
     sat_mask[not_saturated] = 0
     sat_mask[sat_pixels] = 1
     sat_mask = CONVOLVE(sat_mask,circbig)
     goodpix = where(sat_mask lt 1,complement=maskpix)
     sat_mask[maskpix] = 0
     sat_mask[goodpix] = 1
     writefits, tmp_path_s + 'saturation_mask.fits', sat_mask
     accept = []
     for s = 0, n_ref-1 do begin
       xx = round(x_psf[s])
       yy = round(y_psf[s])
       xx = 0 > xx & xx = (x_cube-1) < xx
       yy = 0 > yy & yy = (y_cube-1) < yy
       if (sat_mask[xx,yy] gt 0) then accept = [accept,s]
     endfor
     x_psf = x_psf[accept]
     y_psf = y_psf[accept]
     m_psf = m_psf[accept]
   endif
   n_ref = n_elements(m_psf)
   print, 'Found '+ strn(n_ref) + ' isolated and unsaturated reference stars.'


   ; Sort from brightest to faintest: Important!
   ; ------------------------------------------
   ord = sort(m_psf)
   x_psf = x_psf[ord]
   y_psf = y_psf[ord]
   m_psf = m_psf[ord]
   n_ref = n_elements(m_psf)
   f_psf = 10^(0.4*(ZP-m_psf))
   print, 'Found '+ strn(n_ref) + ' reference stars.'
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
   refstars = image_model(x_psf,y_psf,f_psf,x_cube,y_cube,'gaussian', dat)
   writefits, tmp_path_s + 'refstars.fits', refstars

  ; Run holography
  ; -------------------- 
    filenam = '_' + strn(i_x) + '_' + strn(i_y)
    indir =  in_path  + 'chip' + chip_nr + '/'
    inlist = in_path  + 'chip' + chip_nr + '/list' + filenam + '.txt'
    maskdir = in_path  + 'chip' + chip_nr + '/'
    masklist = maskdir + 'masklist' + filenam + '.txt'

    outnam = out_path + 'chip' + chip_nr + '/' + 'holo_' + strn(i_x) + '_' +  strn(i_y)
    refsources = fltarr(3,n_ref)
    refsources[0,*] = x_psf
    refsources[1,*] = y_psf
    refsources[2,*] = f_psf
    mrad = maskrad
    bhw = boxhw
    nrad = normrad
    n_mask_secondary = 0
    out_psf_path_s =  outpsf_path+ 'chip'  + chip_nr + '/'

    holo_iter = n_elements(nsigma) - 1
    ; using the saturation mask will not
    ; necessarily improve first estimation
    ; of PSF
    ; ESTIM_BG usually not necessary and costs execution time

    mrad = maskrad
    nrad = normrad
    correct_sky = 0
    contrast_thresh = 1
    
    HOLO_FULL, indir, inlist, outnam,  mrad, nrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, DEBUG = debug, iter=holo_iter, AIRY=airy, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, tmpdir=tmp_path_s, BOXHW=boxhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, REBITER=0, RAWOUT = 0, N_MASK_SECONDARY=n_mask_secondary, PSF_FRAC = psf_frac, CORRECT_SKY = correct_sky, SMOOTHMASK = smoothmask, N_REF_MAX = n_ref_max, PR = pr, SATLEVEL=satlevel, PSF_BORDER=psf_border, ESTIM_BG = 0, SAT_MASK = 0, CONTRAST_THRESH=contrast_thresh
    
       
   ; Make noise map from jackknife images
    im = readfits(outnam  + '.fits.gz')
    wt = readfits(outnam + '_expmap.fits.gz')
    sz = size(im)
    n1 = sz[1] & n2 = sz[2]
    imcube = fltarr(n1,n2,nsub)
    meanim = fltarr(n1,n2)
    noise = fltarr(n1,n2)
    map = fltarr(n1,n2)
    valid = where(wt gt 0, complement=invalid)
    map[valid] = 1
    map[invalid] = 0
    for is = 0, nsub-1 do begin
      im = readfits(outnam  + '_s' + strn(is+1) + '.fits.gz')
      imcube[*,*,is] = im
    endfor
    nnjack = (nsub-1)/float(nsub)
    for ix =0, n1-1 do begin
      for iy = 0, n2-1 do begin
        vals = imcube[ix,iy,*]
        mm =  mean(vals)
        meanim[ix,iy] = mm
        noise[ix,iy] = sqrt(total((vals - mm)^2)*nnjack)
      endfor
    endfor
    writefits, outnam + '_ub.fits', meanim, /COMPRESS
    writefits, outnam + '_sigma.fits', noise, /COMPRESS
    writefits, outnam + '_map.fits', map, /COMPRESS

  
   ; Delete all temporary files
   print, 'Finished sub-field ' + strn(i_x) + ', ' + strn(i_y)

   spawn, 'rm ' + tmp_path_s + '*.fits.gz' 
   spawn, 'rm ' + tmp_path_s + '*.fits' 

  endfor
 endfor

 print, 'Finished chip ' + chip_nr

END
