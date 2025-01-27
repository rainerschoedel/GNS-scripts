PRO ALIGNQUADRANTS, field, band, chip

; PURPOSE: Use one (to a few, one is usually sufficient) reference stars per quadrant for a
; pre-alignment with VVV. Then search for more stars and refine the
; alignment.
; 
; This script loops over the four quadrants for each pointing.
; The pointing ('field') must be set by hand.

; ------------EDIT HERE---------
; size of VVV reference image
; This is the printed output  the end of prep_refim.pro

xsize_ref = 1595
ysize_ref = 1000

chip_nr = strn(chip)

ngrid = 20 ; Defines size of star grid for alignment 

; ----------------CAN BE EDITED, but USUALLY NOT NECESSARY ----------


; To create small aligned longexposure images
; So that we can cut away the unnecessary zeros around the images.

; size of dejittered HAWK-I long exposures
xsize_quad = 2400
ysize_quad = 1200
; Approximate offsets of fields within aligned HAWK-I frame
; These offsets can be seen inside the 
; images ;nx_jitter_?_aligned.fits.gz
; that are being produced by this script
; you can stimete aproximate these parametres with the x and y translation from astroalign
x_off = [200,2400,2400,250]
y_off = [650,650,1550,1550]
SAVE, x_off, y_off, FILENAME = 'xy_off_f'+strn(field)+'.SAV'


; VVV --> HAWKI
scale = 0.34/0.106

; define paths
VVV       ='/home/data/VVV/'
im_path   = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/ims/'
data_path = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/data/'
; Use next line for standard alignment
; vvv_stars =  VVV +'/Fields/J/Field' + strn(field) + '_stars.txt'
; Use next line for alignment with high quality stars and PM correction
vvv_stars =  VVV +'/Fields/J/Field' + strn(field) + '_' + band + '_stars_fine.txt'
tmp_path  = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/tmp/'

; Read list of reference stars in VVV
; and scale to HAWK-I pixel scale
; ------------------------------------
readcol, vvv_stars, x_vvv, y_vvv, a, d, J, Ks, FORMAT='F'
;readcol, vvv_stars, x_vvv, y_vvv, a, d, m, sx, sy, sm, c, FORMAT='F'
x_vvv_scaled = x_vvv * scale
y_vvv_scaled = y_vvv * scale

; Define size of HAWK-I mosaic
; -----------------------------
xsize_mosaic = round(xsize_ref * scale) + 50
ysize_mosaic = round(ysize_ref * scale) + 50
print, 'x, y size of HAWK-I mosaic: ' + strn(xsize_mosaic) + ', ' + strn(ysize_mosaic)

;Read list of stars in HAWK-I image
; ----------------------------------
 readcol, data_path + 'stars_' + chip_nr + '.txt', x, y, f
 n_stars = n_elements(f)

; read lists of common stars found by align_VVV.py
 ; -------------------------------------------------
 readcol, data_path + 'aa_stars_hawki_' + chip_nr + '.txt', x_ref, y_ref
 readcol, data_path + 'aa_stars_vvv_' + chip_nr + '.txt', x_ref_vvv, y_ref_vvv
 x_ref_vvv_scaled = x_ref_vvv * scale
 y_ref_vvv_scaled = y_ref_vvv * scale
print, 'Common stars from astroalign: ' + strn(n_elements(x_ref)) + ', ' + strn(n_elements(x_ref_vvv))

; Transform and find common stars
; -------------------------------

 print, 'Initial coarse alignment.'
 degree = 1
 polywarp, x_ref_vvv_scaled, y_ref_vvv_scaled, x_ref, y_ref, degree, Kx, Ky
 print, Kx
 print, Ky
 xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
 yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
 dmax = 2.0
 compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
 nc = n_elements(subc1)
 print, 'Found ' + strn(nc) + ' common stars.'

 
 ; iterative degree 1 alignment
 ; ------------------------------

  degree = 1
  dmax = 1.0
  lim_it = 1
  max_iter = 10
  count=0
  comm=[]
  it=0
  lim_it=2 ;cosider convergence when the number of common stars increases or stasy constant fot 'lim_it' times.
  this_iter = 0
	 
  while (count lt lim_it) and (this_iter lt max_iter) do begin
    it = it+1
    grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
    polywarp, x_vvv_scaled[subc1[grid_idx]], y_vvv_scaled[subc1[grid_idx]], x[subc2[grid_idx]], y[subc2[grid_idx]], degree, Kx, Ky
    print, Kx
    print, Ky
    xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
    yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
    compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
    nc = n_elements(subc1)
    print, 'Iteration ' + strn(it)
    print, 'Found ' + strn(nc) + ' common stars.'
    comm=[comm,nc]
    if (n_elements(comm) gt 2) then begin
      if comm[-2] ge comm[-1] then begin
        count=count+1
      endif else begin
        count = 0
      endelse
    endif
    this_iter = this_iter + 1
  endwhile

 ; iterative degree 2 alignment
 ; ------------------------------

  ; next line serves to switch degree 2 on/off
  if 1 then begin

    degree = 2
    dmax = 1.0
     
   print, 'Now Degree ' + strn(degree) + ' alignment.'
   lim_it = 1
   count=0
   comm=[]
   it=0
	 
  while count lt lim_it do begin
    it=it+1
    grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
    polywarp, x_vvv_scaled[subc1[grid_idx]], y_vvv_scaled[subc1[grid_idx]], x[subc2[grid_idx]], y[subc2[grid_idx]], degree, Kx, Ky

    print, Kx
    print, Ky
    xi = replicate(0.0,n_stars)
    yi = replicate(0.0,n_stars)
    for k = 0, degree do begin
      for m = 0, degree do begin
        xi = xi + Kx[m,k]*x^k*y^m
        yi = yi + Ky[m,k]*x^k*y^m
      endfor
    endfor
;    xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
;    yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x
    compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
    nc = n_elements(subc1)
    print, 'Iteration ' + strn(it)
    print, 'Found ' + strn(nc) + ' common stars.'
    comm=[comm,nc]
    if (n_elements(comm) gt 2) then begin
      if comm[-2] ge comm[-1] then begin
        count=count+1
      endif else begin
	count=0
	endelse
      endif
  endwhile
endif
  

 ; determine transformation parameters for image and save them
 ; for later use
 ; Kx and Ky describe second order polynomial transformation into
 ; a rescaled VVV (to HAWK-I pixel scale) reference frame
 ; ------------------------------------------------------------


 grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
 polywarp,  x[subc2[grid_idx]], y[subc2[grid_idx]], x_vvv_scaled[subc1[grid_idx]], y_vvv_scaled[subc1[grid_idx]], degree, Kx, Ky

  if (degree eq 2) then begin
    x_fin = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
    y_fin = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x 
  endif else begin
    x_fin = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
    y_fin = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
  endelse
  SAVE, Kx, Ky, FILENAME=data_path + 'AlignPars_chip' + chip_nr

 ; transform image and mask
 ; mask pixels that have coverage < 0.9 to avoid 
 ; later problems with division by small numbers
 ; ---------------------------------------------

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_wt.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned.fits.gz', transim, /COMPRESS

 transim_im = transim

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_sig.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0 )
 transim_noise = transim
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_sig.fits.gz', transim, /COMPRESS

 ; to check alignment with VVV, create
 ; a transformed image inside the VVV frame of reference
 ; ---------------------------------------------
 polywarp,   x[subc2[grid_idx]], y[subc2[grid_idx]], x_vvv[subc1[grid_idx]], y_vvv[subc1[grid_idx]], degree, Kx, Ky
 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV_wt.fits.gz', transim, /COMPRESS


; Now cut out the aligned HAWK-I images from the large frames to get
; rid of the edges
; --------------------------------------------------------------------

xlo = x_off[chip-1]
ylo = y_off[chip-1]
xhi = x_off[chip-1] + xsize_quad - 1
yhi = y_off[chip-1] + ysize_quad - 1

lnx = transim_im[xlo:xhi,ylo:yhi]
lnx_noise = transim_noise[xlo:xhi,ylo:yhi]
  
 
writefits, data_path + 'lnx_aligned_' + chip_nr + '.fits.gz', lnx, /COMPRESS
writefits, data_path + 'lnx_aligned_' + chip_nr + '_sig.fits.gz', lnx_noise, /COMPRESS


;To have an idea of the displacement that we have, we take the initial
;list and the corrected one to see the difference in positions.
; These plots can be used to detect systematic residual distortions.
; ------------------------------------------------------------------------------

; median global offset
x_dif = median(x_fin - x)
y_dif = median(y_fin - y) 
 
; origin of arrows in displaced image
x0_new = x_dif + x
y0_new = y_dif + y
 
; max size of image
x_dis = max(x_fin-x0_new)
y_dis = max(y_fin-y0_new)
 

print,  strn(x_dis*0.106) + 'arcsec    ' + strn(y_dis*0.106) + 'arcsec'

dx = x_dis*0.106
dy = y_dis*0.106

;forprint, TEXTOUT= data_path + 'displ_arcsec_distortion_chip' + chip_nr + '.txt', dx, dy, /NOCOMMENT 

;Drawing the detector


set_plot,'PS', /interpolate

device, XOFFSET=0, YOFFSET=0, $
   FILENAME= data_path + 'Alignment_resid_chip' + chip_nr +'.eps', XSIZE=30., YSIZE=30., $
   /portrait, /color, BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,1,0]
!P.CHARSIZE=1.8
!X.THICK=4
!Y.THICK=4
!P.THICK=4.0
!P.CHARTHICK=4

!P.COLOR=0

;cgoplot, x1c, y1c, Color='orange', PSYM=3

; compute shift between common re-scaled VVV and HAWK-I positions
; The arrow will point to the VVV position
; -----------------------------------------------------------------
xf = (-xi[subc2]+x1c)*60 
yf = (-yi[subc2]+y1c)*60
; subtract global shift
xf = xf - median(xf)
yf = yf - median(yf)
x0 = xi[subc2] - xlo
y0 = yi[subc2] - ylo
xf = xf+xi[subc2] - xlo
yf = yf+yi[subc2] - ylo

ran = RANDOMU(Seed, n_elements(xf))

num = sort(ran)
max_num = n_elements(ran) < 80
num = num[0:max_num-1]
a = xi[subc2[num]]
b = yi[subc2[num]]

cgplot, a, b, Color='black', XRANGE = [0,xsize_quad], YRANGE = [0,ysize_quad], XTITLE='X-Axis [pixels]', YTITLE='Y-Axis [pixels]', XSTYLE=1, PSYM=3, YSTYLE=1


ARROW, x0[num], y0[num], xf[num], yf[num], /DATA, COLOR=4, HSIZE=100


device, /close




END
