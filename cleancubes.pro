PRO CLEANCUBES, field, band

;field = '20'
;band = 'Ks'

common_path = '/home/data/GNS/2015/' + band + '/' + field + '/ims/'
in_path = '/home/data/GNS/2015/' + band + '/' + field + '/cubes/'
out_path = '/home/data/GNS/2015/' + band + '/' + field + '/cleaned/'
tmp_path = '/home/data/GNS/2015/' + band + '/' + field + '/tmp/'
mask_name = 'mask.fits'
n_sigma = 7. ; value must be high, otherwise valid pixels of bright stars will be corrected (PSF varies between frames)!
debug = 0
filt_box = 5 ; width of box for sigma filtering



if not(KEYWORD_SET(n_sigma)) then n_sigma = 5.
if not(KEYWORD_SET(box_width)) then box_width = 5
if not(KEYWORD_SET(debug)) then debug = 0
if not(KEYWORD_SET(sigma_dev)) then sigma_dev = 5.



nam = ''

inlist = 'cubelist.txt'
outlist = 'cubelist.txt'


mask = readfits(common_path + mask_name)

openw, lun, (out_path + outlist), /get_lun, lun ; open output file to write
openr, inp, (in_path + inlist), /get_lun  ; open input file for reading
cn = 0L
while (not (EOF(inp))) do begin
   readf, inp, nam
   cube = readfits(in_path + nam, header)
   sz = size(cube)
   n3 = sz[3]

   ;  Fid and interpolate bad pixels
   for j = 0, n3 -1 do begin
      im = cube[*,*,j]

      if debug then writefits, tmp_path + 'im_raw.fits', im
      filtered_im = sigma_filter(im, filt_box, N_SIGMA=n_sigma, /ITERATE) * mask ; multiply with mask to avoid non-zero values near mask edges
      bad = where((im - filtered_im) ne 0, n_bad)
      if (n_bad gt 0) then begin
         print, 'Found ' + strn(n_bad) + ' bad pixels in image ' + strn(j+1) + ' of cube ' + nam
         xy = array_indices(im,bad)
         xbad = xy[0,*]
         ybad = xy[1,*]
         im = replace_pix(im,xbad,ybad)
         nok = im
         nok[*,*] = 0
         nok[bad] = 1
         if debug then begin
            writefits, tmp_path + 'clean.fits', im
            writefits, tmp_path + 'nok.fits', nok
            STOP
         endif
      endif
      cube[*,*,j] = im * mask
   endfor
   cn = cn + 1
   writefits, out_path + 'cube'+ strn(cn) + '.fits', cube, header, /COMPRESS
 
   printf, lun, 'cube' + strn(cn) + '.fits'
   print, nam
endwhile

free_lun, inp
free_lun, lun

print, "cleancubes.pro ended"

END
