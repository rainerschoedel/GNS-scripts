PRO MOSAIC, field, band


refim = '/home/data/VVV/Fields/J/Field' + field + '.fits.gz'
im_path = '/home/data/GNS/2015/'+band+'/' + field + '/ims/'

; mosaic at VVV pixel scale
; see xsize_ref and ysize_ref in alignquadrants.pro
x_VVV =  1595
y_VVV = 1000
mosaic = fltarr(x_VVV,y_VVV)
weight = fltarr(x_VVV,y_VVV)

for chip = 1, 4 do begin

 chip_nr = strn(chip)
 im = readfits(im_path + 'lnx_jitter_'+strn(chip_nr)+ '_VVV.fits.gz')
 wt = readfits(im_path + 'lnx_jitter_'+strn(chip_nr) + '_VVV_wt.fits.gz')
 good = where(wt gt 3,complement=bad)

; wt[good] = 1
 wt[bad] = 0
 im = im * wt
 mosaic = mosaic + im
 weight = weight + wt

endfor
good = where(weight gt 0,complement=bad)
mosaic[good] = mosaic[good]/weight[good]
writefits, 'mosaic.fits.gz', mosaic, /COMPRESS
find, mosaic, xx, yy, flux, sharp, roundness, 1000.0, 4.0, [-1.0,1.0], [0.2,1.0]
print, 'Found ' + strn(n_elements(xx)) + ' stars for fine alignment.'

; shift the images, which have a small offset compared to the VVV
; reference image
ref = readfits(refim)
find, ref, xx_ref, yy_ref, flux, sharp, roundness, 1000.0, 4.0, [-1.0,1.0], [0.2,1.0]
print, 'Found ' + strn(n_elements(xx_ref)) + ' stars for fine alignment.'

; determine offset from lists of detected sources
dmax = 4.0
compare_lists, xx_ref, yy_ref, xx, yy, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2
nc = n_elements(x1c)
print, 'Found ' + strn(nc) + ' common stars.'
x_off = median(x1c - x2c)
y_off = median(y1c - y2c)
print, 'Offsets: ' + strn(x_off) + ', ' + strn(y_off)


;ref[bad] = 0
;correl_optimize, ref, smooth(mosaic,1), x_off, y_off, MAGNIFICATION=4, XOFF_INIT=0, YOFF_INIT=0
;print, x_off, y_off

writefits, im_path + 'lnx_mosaic_VVV.fits.gz', image_shift(mosaic,x_off,y_off), /COMPRESS
writefits, im_path + 'lnx_weight_VVV.fits.gz', image_shift(weight,x_off, y_off), /COMPRESS

; mosaic at HAWK-I pixel scale
x_hawki = 5166
y_hawki = 3258
mosaic = fltarr(x_hawki, y_hawki)
mosaic_sigma = fltarr(x_hawki, y_hawki)
weight = fltarr(x_hawki, y_hawki)

for chip = 1, 4 do begin

 chip_nr = strn(chip)
 im = readfits(im_path + 'lnx_jitter_'+strn(chip_nr)+ '_aligned.fits.gz')
 sigma = readfits(im_path + 'lnx_jitter_'+strn(chip_nr)+ '_aligned_sig.fits.gz')
 wt = readfits(im_path + 'lnx_jitter_'+strn(chip_nr) + '_aligned_wt.fits.gz')
 not_nan =  where(FINITE(sigma), complement=isnan)
 sigma[isnan] = 0

; exclude pixels with low weights or
; inexistent sigma (nan value in sigma array)
 good = where(wt gt 3 and sigma gt 0,complement=bad)
; wt[good] = 1
 wt[bad] = 0
 im = im * wt
 sigma2 = sigma^2 * wt
 mosaic = mosaic + im
 mosaic_sigma = mosaic_sigma + sigma2
 weight = weight + wt

endfor
good = where(weight gt 0)
mosaic[good] = mosaic[good]/weight[good]
mosaic_sigma[good] = sqrt(mosaic_sigma[good]/weight[good])
writefits, im_path + 'lnx_mosaic.fits.gz', mosaic, /COMPRESS
writefits, im_path + 'lnx_mosaic_sigma.fits.gz', mosaic_sigma, /COMPRESS
writefits, im_path + 'lnx_weight.fits.gz', weight, /COMPRESS

END
