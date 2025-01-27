pro dejitter, field, band, chip

basedir = '/home/data/GNS/2015/'+band+'/' + strn(field)
indir = basedir + '/ims/'
outdir = basedir + '/cubes/'
datadir = basedir + '/data/'

; Choose final image size large enough to accommodate jittered
; images for 2048 x 768 pixels with 60" jitter box
; in x: (2048 * 0.106 + 60.)/0.106 = 2614.04
; in y: (768 * 0.106 + 60.)/0.106 = 1334.04
; Later on, the images will be chopped up in so-called sub-cubes.
; The sub-cubes are 600 x 600 pixels and are shifted by 300 pixels
; To make this fit nicely, we now round the above numbers to
; 2700 and 1500
size_old_x = 2048
size_old_y = 768
sizex_new = 2700
sizey_new = 1500

; For jitter box 30" use these numbers: comment for jitter box 60" 
;sizex_new = 2400
;sizey_new = 1200
;
; The jitter_offset paraemter makes sure that an exposure
; with zero offset falls into the middle of the frame after alignment
jitter_offset = round((60./2.)/0.106) ; maximal offset of a jittered exposure from initial pointing
; For jitter box 30 use this number
;jitter_offset = round((30./2.)/0.106) ; comment for jitter box 60"
;jitter_offset = round((30./2.)/0.106)+1.0 ; use only for F9 reduction, comment otherwise
;jitter offset = 330 ; This value has been chosen empirically (systematics in FITS Header cumoffset?)

n_offsets = 48
openw, out1, outdir + 'list.txt', /get_lun
openw, out2, outdir + 'masklist.txt', /get_lun
openw, out3, outdir + 'jitter_offsets.txt', /get_lun


lnx_cube = fltarr(sizex_new,sizey_new,n_offsets)
mask_cube = fltarr(sizex_new,sizey_new,n_offsets)

for i = 1, n_offsets do begin

  print, 'Pointing : ' + strn(i)   

  printf, out1, 'cube_jitter_' + strn(i) + '.fits.gz'
  printf, out2, 'mask_jitter_' + strn(i) + '.fits.gz'

  ; Determine the initial offset fro exposure 1
  ; Then determine jitters from offsets of stellar positions in chip1 and correct al chips
  
  if i eq 1 then begin
     
    cube = readfits(indir + 'chip1_cube' + strn(i) + '.fits.gz', header)
    mask = readfits(indir + 'mask_chip1.fits.gz')
    sz = size(cube)
    elements_cube = sz[3]

    ; Read offset from header
    x_ini = strsplit(header[427],'HIERARCH ESO SEQ CUMOFFSETX = ', ESCAPE = '/', /extract)
    y_ini = strsplit(header[428],'HIERARCH ESO SEQ CUMOFFSETY = ', ESCAPE = '/', /extract) 
; use header[608-609] only if reducing Ks band from field 7 (2018
; observations), comment otherwise and use 427 and 428
    x_ini = fix(x_ini[0])
    y_ini = fix(y_ini[0])
    print, 'Offsets from Fits Header: ' + strn(x_ini) + ', ' + strn(y_ini)
    ref_header = header
  
    pos_x = jitter_offset - x_ini
    pos_y = jitter_offset - y_ini

    x_off = 0
    y_off = 0
    printf, out3, x_off, y_off, FORMAT='(2F8.2)'

  endif else begin

   ; Determine offset from common stars
   ; identified by align_dejitter.py
   readcol, datadir + 'dejitter_common_ref_1' + '_' + strn(i) + '.txt', xref, yref, /SILENT
   readcol, datadir + 'dejitter_common_1' + '_' + strn(i) + '.txt', xim, yim, /SILENT

   x_off = long(median(xref - xim))
   y_off = long(median(yref - yim))
   print, 'Median offsets with respect to exposure 1 in x and y: ' + strn(x_off) + ', ' + strn(y_off)
   printf, out3, x_off, y_off, FORMAT='(2F8.2)'

 endelse
     
 ; Now dejitter chip

  cube = readfits(indir + 'chip' + strn(chip) + '_cube' + strn(i) + '.fits.gz', header)
  mask = readfits(indir + 'mask_chip' + strn(chip) + '.fits.gz')
  sz = size(cube)
  elements_cube = sz[3]
  new_cube = fltarr(sizex_new,sizey_new, elements_cube-1)

  new_mask = fltarr(sizex_new,sizey_new)
  lnx_tmp = fltarr(sizex_new,sizey_new)
    
  lxp = MEAN(cube[*,*,0:elements_cube-2],DIM=3) ; Do not use long exposure from ESO cube, which may have errors (in case of unwanted telescpoce offsets)
  
  lnx_tmp[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1] = lxp
  new_mask[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1] = mask
  mask_cube[*,*,i-1] = new_mask
  
  new_cube[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1,*] = cube[*,*,0:elements_cube-2]
  
  lnx_cube[*,*,i-1] = lnx_tmp[*,*]
  mask_cube[*,*,i-1] = new_mask[*,*]   

  writefits, outdir + 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_cube, /COMPRESS
  writefits, outdir + 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_mask, /COMPRESS
  print, 'Dejittered chip ' + strn(chip) + '.'

endfor

; output: loongexposure cubes for chip1 1 to 4
writefits, outdir + 'lnx_jitter_cube_' + strn(chip) + '.fits', lnx_cube, ref_header, /COMPRESS
writefits, outdir + 'lnx_jitter_mask_' + strn(chip) + '.fits', mask_cube, /COMPRESS

free_lun, out1, out2, out3

;stop
 

end
