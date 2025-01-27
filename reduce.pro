PRO REDUCE, field, band

debug = 0
nam = ''

;field = 'B1'
;band = 'Ks'

raw_path = '/home/data/raw/2015/' + band +'/Field/' + field + '/'
common_path = '/home/data/GNS/2015/' + band + '/' + field + '/ims/'
tmp_path = '/home/data/GNS/2015/' + band + '/' + field + '/tmp/'
out_path = '/home/data/GNS/2015/' + band + '/' + field + '/cubes/'
mask_name = 'mask.fits'
bpm_name = 'bpm.fits'
sky_name = 'sky.fits'
dark_name = 'dark.fits'
flat_name = 'flat_' + band + '.fits'
list = 'list.txt'
outlist = 'cubelist.txt'

readcol, common_path + 'gains.txt', gains
sky_frac = 0.1     ; fraction of pixels of lowest value
                   ; that will be used to create the matched sky
                   ; The MMM algorithm is relatively insensitive 
                   ; to this value

dark = readfits(common_path + dark_name)
sky = readfits(common_path + sky_name)
bpm = readfits(common_path + bpm_name)
mask = readfits(common_path + mask_name)
flat = readfits(common_path + flat_name)

; work on the four chips (quadrants) individually
;sky_qs = fltarr(2048,768,4)
;sky_qs[*,*,0] = sky[0:2047,0:767]
;sky_qs[*,*,1] = sky[2048:4095,0:767]
;sky_qs[*,*,2] = sky[2048:4095,768:1535]
;sky_qs[*,*,3] = sky[0:2047,768:1535]


; identify bad pixels for interpolation
bad = where(bpm gt 0 and mask gt 0)
xy = array_indices(bpm,bad)
xbad = xy[0,*]
ybad = xy[1,*]

; identify good pixels for entire image and for each quadrant
good = where(bpm lt 1 and mask gt 0)
;good_q1 = where(bpm[0:2047,0:767] lt 1 and mask[0:2047,0:767] gt 0)
;good_q2 = where(bpm[2048:4095,0:767] lt 1 and mask[2048:4095,0:767] gt 0)
;good_q3 = where(bpm[2048:4095,768:1535] lt 1 and mask[2048:4095,768:1535] gt 0)
;good_q4 = where(bpm[0:2047,768:1535] lt 1 and mask[0:2047,768:1535] gt 0)
;good_qs = [Ptr_New(good_q1, /No_Copy), Ptr_New(good_q2, /No_Copy), Ptr_New(good_q3, /No_Copy), Ptr_New(good_q4, /No_Copy)]

openw, lun, (out_path + outlist), /get_lun, lun ; open output file to write

openr, inp, (raw_path + list), /get_lun  ; open input file for reading


cn = 0L
while (not (EOF(inp))) do begin
   readf, inp, nam
   cube = readfits(raw_path + nam, header) 
   sz = size(cube)
   n3 = sz[3]
   ; subtract dark
   for i_c = 0, n3-1 do begin
     cube[*,*,i_c] = cube[*,*,i_c] - dark
   endfor
;   writefits, tmp_path + 'raw.fits', cube

   ; apply gain factors to the four quadrants
   ; gain factors are computed in sky.pro
   cube[0:2047,0:767,*] = cube[0:2047,0:767,*] / gains[0]
   cube[2048:4095,0:767,*] = cube[2048:4095,0:767,*] / gains[1]
   cube[2048:4095,768:1535,*] = cube[2048:4095,768:1535,*] / gains[2]
   cube[0:2047,768:1535,*] = cube[0:2047,768:1535,*] / gains[3]
;   writefits, tmp_path + 'rawxgains.fits', cube

   ; Match the master sky to the sky level in this cube
   ; With the gain differences between the chips
   ; already having been taken care of in sky.pro and 
   ; above in this script, I can use the entire field 
   ; to determine the sky offset
   ; I can thus avoid scaling problems due to 
   ; the presence of crowded fields on some chips
;   lng = cube[*,*,n3-1] ; last image in cube is long exposure
;   lng[good] = lng[good]/flat[good]
;   lng = replace_pix(lng,xbad,ybad)
;   skyvals = sky[good]
;   vals = lng[good]
;   ord = sort(vals)
;   vals = vals[ord]
;   skyvals = skyvals[ord]
;   skythresh = vals[round(sky_frac*n_elements(vals)) - 1]
;   sky_ind = where(vals lt skythresh,skycount)
;   mmm, skyvals[sky_ind], skymod, skysigma, skyskew 
;   mmm, vals[sky_ind], imskymod, imskysigma, imskyskew 
;   factor = imskymod/skymod
;   this_sky = sky * factor
;   print, 'Sky factor:' 
;   print, factor
;   print, 'Relative uncertainty: '
;   print, strn(sqrt((skysigma/skymod)^2  +(imskysigma/imskymod)^2))
;   print

   ; Now reduce the frames
   ; dark has already been subtracted above
   ; gain has already been calibrated above
   for j = 0, n3-1 do begin
      im = cube[*,*,j]
 ;     writefits, tmp_path + 'raw.fits', im
      
      ; scale sky
      im[good] = im[good]/flat[good]
      im = replace_pix(im,xbad,ybad)
      skyvals = sky[good]
      vals = im[good]
      ord = sort(vals)
      vals = vals[ord]
      skyvals = skyvals[ord]
      skythresh = vals[round(sky_frac*n_elements(vals)) - 1]
      sky_ind = where(vals lt skythresh,skycount)
      mmm, skyvals[sky_ind], skymod, skysigma, skyskew 
      mmm, vals[sky_ind], imskymod, imskysigma, imskyskew 
      factor = imskymod/skymod
      this_sky = sky * factor
      print, 'Sky factor:' 
      print, factor
      print, 'Relative uncertainty: '
      print, strn(sqrt((skysigma/skymod)^2  +(imskysigma/imskymod)^2))
      print

;      writefits, tmp_path + 'this_sky.fits', this_sky
      im = cube[*,*,j]
      im = im - this_sky
      im[good] = im[good]/flat[good]
      im = replace_pix(im,xbad,ybad)
  ;    writefits, tmp_path + 'reduced.fits', im * mask
      cube[*,*,j] = im * mask

   endfor
   cn = cn + 1

   sxaddpar, header, 'NAXIS3', n3
   sxdelpar, header, 'CHECKSUM'
   sxdelpar, header, 'CDATASUM'
   writefits, out_path + 'cube' + strn(cn) + '.fits', cube, header, /COMPRESS

;   writefits, tmp_path + 'reduced.fits', cube
 
;   writefits,tmp_path + 'long_sigma.fits', stddev(cube,DIMENSION=3)

   printf, lun, 'cube' + strn(cn) + '.fits.gz'

   print, nam
endwhile

free_lun, inp
free_lun, lun

;stop

print, "reduce.pro ended"

END
